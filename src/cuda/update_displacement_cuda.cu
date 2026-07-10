/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                              CNRS, France
 !                       and Princeton University, USA
 !                 (there are currently many more authors!)
 !                           (c) October 2017
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 3 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#include "mesh_constants_cuda.h"
#include <cublas_v2.h>

/* ----------------------------------------------------------------------------------------------- */

// elastic wavefield

/* ----------------------------------------------------------------------------------------------- */


// __global__ void UpdateDispVeloc_kernel(realw_p displ,
//                                        realw_p veloc,
//                                        realw_p accel,
//                                        int size,
//                                        realw deltat,
//                                        realw deltatsqover2,
//                                        realw deltatover2) {

//   // two dimensional array of blocks on grid where each block has one dimensional array of threads
//   int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  
//   // because of block and grid sizing problems, there is a small
//   // amount of buffer at the end of the calculation
//   if (id < size) {
//     realw acc = accel[id];
//     displ[id] = displ[id] + deltat*veloc[id] + deltatsqover2*acc;
//     veloc[id] = veloc[id] + deltatover2*acc;
//     accel[id] = 0.0f; // can do this using memset...not sure if faster,probably not
//   }

// // -----------------
// // total of: 6 FLOP per thread (without int id calculation at beginning)
// //
// //           8 * 4 BYTE = 32 DRAM accesses per thread
// //
// // arithmetic intensity: 6 FLOP / 32 BYTES ~ 0.19 FLOP/BYTE
// // -----------------
// // nvprof: 24599250 flops for 4099875 threads -> 6 FLOP per thread
// }

__global__ void UpdateDispVeloc_kernel(
  realw_p displ,
  realw_p veloc,
  realw_p accel,
  int size,
  realw deltat,
  realw deltatsqover2,
  realw deltatover2
) {

    // 1. Calculate ID for the VECTOR (chunk of 4), not the individual float
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    // 2. Main Vectorized Loop (Processes 4 items at a time)
    // We treat 'size' as the number of floats, so we check if the *entire vector* fits.
    int idx = tid * 4;
    
    if (idx + 4 <= size) {
        // Load data as 128-bit vectors
        float4 d_vec = reinterpret_cast<float4*>(displ)[tid];
        float4 v_vec = reinterpret_cast<float4*>(veloc)[tid];
        float4 a_vec = reinterpret_cast<float4*>(accel)[tid];

        // Perform computation on all 4 components (.x, .y, .z, .w)
        
        // Update Displacement
        d_vec.x += deltat * v_vec.x + deltatsqover2 * a_vec.x;
        d_vec.y += deltat * v_vec.y + deltatsqover2 * a_vec.y;
        d_vec.z += deltat * v_vec.z + deltatsqover2 * a_vec.z;
        d_vec.w += deltat * v_vec.w + deltatsqover2 * a_vec.w;

        // Update Velocity
        v_vec.x += deltatover2 * a_vec.x;
        v_vec.y += deltatover2 * a_vec.y;
        v_vec.z += deltatover2 * a_vec.z;
        v_vec.w += deltatover2 * a_vec.w;

        // Write back results
        reinterpret_cast<float4*>(displ)[tid] = d_vec;
        reinterpret_cast<float4*>(veloc)[tid] = v_vec;
        
        // Zero out acceleration (Write 0.0f to all 4 components)
        // This generates a single ST.128 instruction of zeros.
        float4 zero_vec = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        reinterpret_cast<float4*>(accel)[tid] = zero_vec;
    }
    // 3. Tail Handling (Edge case for non-multiples of 4)
    // Only the very last active thread will ever enter this block.
    else if (idx < size) {
        for (int k = 0; k < 4; ++k) {
            int current_idx = idx + k;
            if (current_idx < size) {
                realw acc = accel[current_idx];
                displ[current_idx] += deltat * veloc[current_idx] + deltatsqover2 * acc;
                veloc[current_idx] += deltatover2 * acc;
                accel[current_idx] = 0.0f;
            }
        }
    }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(update_displacement_cuda,
              UPDATE_DISPLACMENT_CUDA)(long* Mesh_pointer,
                                          realw* deltat_F,
                                          realw* deltatsqover2_F,
                                          realw* deltatover2_F,
                                          realw* b_deltat_F,
                                          realw* b_deltatsqover2_F,
                                          realw* b_deltatover2_F) {

  TRACE("\tupdate_displacement_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  int size = (NDIM * mp->NGLOB_AB + 3) / 4; // number of float4 elements
  int nblocks = (size + BLOCKSIZE_KERNEL1 - 1) / BLOCKSIZE_KERNEL1;
  dim3 grid(nblocks,1,1);
  dim3 threads(BLOCKSIZE_KERNEL1,1,1);

  // Cuda timing
  cudaEvent_t start,stop;
  if (CUDA_TIMING_UPDATE ){
    start_timing_cuda(&start,&stop);
  }

  // debug
  //realw max_d,max_v,max_a;
  //max_d = get_device_array_maximum_value(mp->d_displ, size);
  //max_v = get_device_array_maximum_value(mp->d_veloc, size);
  //max_a = get_device_array_maximum_value(mp->d_accel, size);
  //printf("rank %d - max displ: %f veloc: %f accel: %f\n",mp->myrank,max_d,max_v,max_a);

  //launch kernel
  UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                NDIM*mp->NGLOB_AB,deltat,deltatsqover2,deltatover2);

  // kernel for backward fields
  if (mp->simulation_type == 3) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;

    UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                                                  NDIM*mp->NGLOB_AB,b_deltat,b_deltatsqover2,b_deltatover2);
  }

  // Cuda timing
  if (CUDA_TIMING_UPDATE ){
    realw flops,time;
    stop_timing_cuda(&start,&stop,"UpdateDispVeloc_kernel",&time);
    // time in seconds
    time = time / 1000.;
    // performance: 6 FLOPS per thread
    flops = 6.0 * size;
    //printf("  performance: %f GFlop/s num_blocks x/y: %d %d threads: %d\n", flops/time * 1.e-9,num_blocks_x,num_blocks_y,size);
    printf("  performance: %f GFlop/s\n", flops/time * 1.e-9);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("update_displacement_cuda");
#endif 
} 

extern "C" void 
FC_FUNC_(update_displacement_cuda_ade,
        UPDATE_DISPLACEMENT_CUDA_ADE) (
                                          long* Mesh_pointer,
                                          realw* deltat_F,
                                          realw* deltatsqover2_F,
                                          realw* deltatover2_F,
                                          realw* b_deltat_F,
                                          realw* b_deltatsqover2_F,
                                          realw* b_deltatover2_F) 
{

  TRACE("\tupdate_displacement_cuda_ade");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  int size = (NDIM * mp->NGLOB_AB + 3) / 4; // number of float4 elements
  int nblocks = (size + BLOCKSIZE_KERNEL1 - 1) / BLOCKSIZE_KERNEL1;
  dim3 grid(nblocks,1,1);
  dim3 threads(BLOCKSIZE_KERNEL1,1,1);

  UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                NDIM*mp->NGLOB_AB,deltat,deltatsqover2,deltatover2);

  if(mp->simulation_type == 3 && (!mp->SUBSAMPLE_FWD_WAVEFIELD)) {
    UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(
      mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
      NDIM*mp->NGLOB_AB,*b_deltat_F,*b_deltatsqover2_F,*b_deltatover2_F);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// acoustic wavefield

// KERNEL 1
/* ----------------------------------------------------------------------------------------------- */

__global__ void UpdatePotential_kernel(field* potential_acoustic,
                                       field* potential_dot_acoustic,
                                       field* potential_dot_dot_acoustic,
                                       int size,
                                       realw deltat,
                                       realw deltatsqover2,
                                       realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    field p_dot_dot = potential_dot_dot_acoustic[id];
    field p_dot = potential_dot_acoustic[id];
    potential_acoustic[id] +=   deltat*p_dot
                              + deltatsqover2*p_dot_dot;

    potential_dot_acoustic[id] = p_dot + deltatover2*p_dot_dot;

    potential_dot_dot_acoustic[id] = Make_field(0.f);
  }

// -----------------
// total of: 6 FLOP per thread (without id calculation)
//
//           8 * 4 BYTE = 32 DRAM accesses per thread
//
// arithmetic intensity: 6 FLOP / 32 BYTES ~ 0.19 FLOP/BYTE
// -----------------
//
// nvprof: nvprof --metrics flops_sp ./xspecfem3D
//          -> 8199750 FLOPS (Single) floating-point operations for 1366625 threads
//                                    1366625 (NGLOB) -> 10677 * 128 active threads- 31 ghost threads
//          -> 6 FLOP per thread


}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(it_update_displacement_ac_cuda,
              it_update_displacement_ac_cuda)(long* Mesh_pointer,
                                               realw* deltat_F,
                                               realw* deltatsqover2_F,
                                               realw* deltatover2_F,
                                               realw* b_deltat_F,
                                               realw* b_deltatsqover2_F,
                                               realw* b_deltatover2_F) {
  TRACE("\tit_update_displacement_ac_cuda");
  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //launch kernel
  // forward wavefields
  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  // Cuda timing
  cudaEvent_t start,stop;
  if (CUDA_TIMING_UPDATE ){
    start_timing_cuda(&start,&stop);
  }

  UpdatePotential_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                                                 mp->d_potential_dot_acoustic,
                                                                 mp->d_potential_dot_dot_acoustic,
                                                                 size,deltat,deltatsqover2,deltatover2);

  // backward/reconstructed wavefields
  if (mp->simulation_type == 3) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;

    UpdatePotential_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,
                                                                  mp->d_b_potential_dot_acoustic,
                                                                  mp->d_b_potential_dot_dot_acoustic,
                                                                  size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

  // Cuda timing
  if (CUDA_TIMING_UPDATE ){
    realw flops,time;
    stop_timing_cuda(&start,&stop,"UpdatePotential_kernel",&time);
    // time in seconds
    time = time / 1000.;
    // performance
    // see with: nvprof --metrics flops_sp ./xspecfem3D
    //           -> using 8199750 FLOPS (Single) floating-point operations for 1366625 threads
    //              = 6 FLOPS per thread
    flops = 6.0 * size;
    //printf("  performance: %f GFlop/s num_blocks x/y: %d %d threads: %d\n", flops/time * 1.e-9,num_blocks_x,num_blocks_y,size);
    printf("  performance: %f GFlop/s\n", flops/time * 1.e-9);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("it_update_displacement_ac_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// elastic domains

// KERNEL 3

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_cuda_device(realw* veloc,
                                     realw* accel,
                                     realw* b_veloc,
                                     realw* b_accel,
                                     int size,
                                     int simulation_type,
                                     realw deltatover2,
                                     realw b_deltatover2,
                                     realw* rmassx,
                                     realw* rmassy,
                                     realw* rmassz) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  realw rx,ry,rz;
  realw ax,ay,az;
  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    rx = rmassx[id];
    ry = rmassy[id];
    rz = rmassz[id];
    ax = accel[3*id  ]*rx;
    ay = accel[3*id+1]*ry;
    az = accel[3*id+2]*rz;

    accel[3*id]   = ax;
    accel[3*id+1] = ay;
    accel[3*id+2] = az;

    veloc[3*id]   += deltatover2*ax;
    veloc[3*id+1] += deltatover2*ay;
    veloc[3*id+2] += deltatover2*az;

    if (simulation_type==3){
      ax = b_accel[3*id  ]*rx;
      ay = b_accel[3*id+1]*ry;
      az = b_accel[3*id+2]*rz;

      b_accel[3*id]   = ax;
      b_accel[3*id+1] = ay;
      b_accel[3*id+2] = az;

      b_veloc[3*id]   += b_deltatover2*ax;
      b_veloc[3*id+1] += b_deltatover2*ay;
      b_veloc[3*id+2] += b_deltatover2*az;

    }

  }

}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_accel_cuda_device(realw* accel,
                                           realw* b_accel,
                                           int size,
                                           int simulation_type,
                                           realw* rmassx,
                                           realw* rmassy,
                                           realw* rmassz) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  realw rx,ry,rz;
  realw ax,ay,az;
  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    rx = rmassx[id];
    ry = rmassy[id];
    rz = rmassz[id];
    ax = accel[3*id  ]*rx;
    ay = accel[3*id+1]*ry;
    az = accel[3*id+2]*rz;
    accel[3*id  ] = ax;
    accel[3*id+1] = ay;
    accel[3*id+2] = az;

    if (simulation_type==3){
      ax = b_accel[3*id  ]*rx;
      ay = b_accel[3*id+1]*ry;
      az = b_accel[3*id+2]*rz;

      b_accel[3*id]   = ax;
      b_accel[3*id+1] = ay;
      b_accel[3*id+2] = az;

    }

  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_veloc_cuda_device(realw* veloc,
                                           realw* accel,
                                           int size,
                                           realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id+1] = veloc[3*id+1] + deltatover2*accel[3*id+1];
    veloc[3*id+2] = veloc[3*id+2] + deltatover2*accel[3*id+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_a_cuda,
              KERNEL_3_A_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F,
                               int* APPROXIMATE_OCEAN_LOAD) {

  TRACE("\tkernel_3_a_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // check whether we can update accel and veloc, or only accel at this point
  if (*APPROXIMATE_OCEAN_LOAD == 0){
   realw deltatover2 = *deltatover2_F;
   realw b_deltatover2 = *b_deltatover2_F;
   // updates both, accel and veloc
   kernel_3_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                                 mp->d_accel,
                                                                 mp->d_b_veloc,
                                                                 mp->d_b_accel,
                                                                 size,mp->simulation_type,deltatover2,b_deltatover2,
                                                                 mp->d_rmassx,mp->d_rmassy,mp->d_rmassz);
  }else{
   // updates only accel
   kernel_3_accel_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                       mp->d_b_accel,
                                                                       size,
                                                                       mp->simulation_type,
                                                                       mp->d_rmassx,
                                                                       mp->d_rmassy,
                                                                       mp->d_rmassz);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel 3 a");
#endif
}



extern "C"
void FC_FUNC_(kernel_3_c_cuda,
              KERNEL_3_C_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F,
                               int* APPROXIMATE_OCEAN_LOAD) {

  TRACE("\tkernel_3_c_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // check whether we can update accel and veloc, or only accel at this point
  if (*APPROXIMATE_OCEAN_LOAD == 0){
   realw deltatover2 = *deltatover2_F;
   realw b_deltatover2 = *b_deltatover2_F;
   // updates both, accel and veloc
   kernel_3_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                                 mp->d_accel,
                                                                 mp->d_b_veloc,
                                                                 mp->d_b_accel,
                                                                 size,1,deltatover2,b_deltatover2,
                                                                 mp->d_rmassx,mp->d_rmassy,mp->d_rmassz);
  }else{
   // updates only accel
   kernel_3_accel_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                       mp->d_b_accel,
                                                                       size,
                                                                       1,
                                                                       mp->d_rmassx,
                                                                       mp->d_rmassy,
                                                                       mp->d_rmassz);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel 3 a");
#endif
}

// __global__ void kernel_apply_mass(realw_p accel, realw_const_p rmassx,
//                                   realw_const_p rmassy,realw_const_p rmassz,
//                                   int nglob)
// {
//   int idx = threadIdx.x + blockIdx.x * blockDim.x;
//   if (idx < nglob) {
//     // Using registers instead of repeated memory access
//     realw rx = rmassx[idx];
//     realw ry = rmassy[idx];
//     realw rz = rmassz[idx];

//     int base_idx = idx * NDIM;
    
//     accel[base_idx]     *= rx;
//     accel[base_idx + 1] *= ry;
//     accel[base_idx + 2] *= rz;
//   }
// }

__global__ void kernel_apply_mass(
  realw_p accel, realw_const_p rmassx,
  realw_const_p rmassy,realw_const_p rmassz,
  int nglob
) 
{
    // Process 4 particles per thread
    int idx = (threadIdx.x + blockIdx.x * blockDim.x);
    int p_idx = idx * 4; // Particle index start

    // --- SAFETY CHECK ---
    // Ensure we don't read past the end. 
    // Optimization applies only if we have a full chunk of 4 particles.
    if (p_idx + 4 <= nglob) {
        
        // 1. Load Masses (Structure of Arrays - Easy)
        // We load 4 masses at once using float4
        float4 mx = reinterpret_cast<const float4*>(rmassx)[idx];
        float4 my = reinterpret_cast<const float4*>(rmassy)[idx];
        float4 mz = reinterpret_cast<const float4*>(rmassz)[idx];

        // 2. Load Accel (Array of Structures - Tricky)
        // 4 particles * 3 components = 12 floats.
        // We load them as 3 contiguous float4 vectors.
        // The accel pointer is treated as float4*.
        // Base index in float4 space is 3 * idx.
        float4* accel_vec_ptr = reinterpret_cast<float4*>(accel);
        
        float4 v1 = accel_vec_ptr[3 * idx + 0]; // x0, y0, z0, x1
        float4 v2 = accel_vec_ptr[3 * idx + 1]; // y1, z1, x2, y2
        float4 v3 = accel_vec_ptr[3 * idx + 2]; // z2, x3, y3, z3

        // 3. Math (Apply Mass) - "Shuffle" the data in registers
        // Particle 0
        v1.x *= mx.x; // x0
        v1.y *= my.x; // y0
        v1.z *= mz.x; // z0
        
        // Particle 1
        v1.w *= mx.y; // x1
        v2.x *= my.y; // y1
        v2.y *= mz.y; // z1
        
        // Particle 2
        v2.z *= mx.z; // x2
        v2.w *= my.z; // y2
        v3.x *= mz.z; // z2
        
        // Particle 3
        v3.y *= mx.w; // x3
        v3.z *= my.w; // y3
        v3.w *= mz.w; // z3

        // 4. Store Back
        accel_vec_ptr[3 * idx + 0] = v1;
        accel_vec_ptr[3 * idx + 1] = v2;
        accel_vec_ptr[3 * idx + 2] = v3;
    }
    else {
        // --- TAIL HANDLING ---
        // Fallback for the last 1-3 particles if nglob is not div by 4
        for (int k = 0; k < 4; ++k) {
            int curr = p_idx + k;
            if (curr < nglob) {
                realw rx = rmassx[curr];
                realw ry = rmassy[curr];
                realw rz = rmassz[curr];

                int base = curr * 3; // NDIM = 3
                accel[base]     *= rx;
                accel[base + 1] *= ry;
                accel[base + 2] *= rz;
            }
        }
    }
}

extern "C"
void apply_massmat_device_(long* Mesh_pointer,const int *backward)
{
  TRACE("\tapply_massmat_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = (mp->NGLOB_AB + 3) / 4; // number of float4 elements
  //size = mp->NGLOB_AB;

  int blocksize = 128;
  int nb = (size + blocksize - 1) / blocksize;


  if(*backward) {
    kernel_apply_mass <<<nb,blocksize,0,mp->compute_stream>>> (
      mp->d_b_accel,mp->d_rmassx,mp->d_rmassy,mp->d_rmassz,mp->NGLOB_AB
    );
  }
  else {
    kernel_apply_mass <<<nb,blocksize,0,mp->compute_stream>>> (
      mp->d_accel,mp->d_rmassx,mp->d_rmassy,mp->d_rmassz,mp->NGLOB_AB
    );
  }

}

__global__ void kernel_UpdateVeloc(realw_p veloc, const realw *accel,
                                   realw dtover2,int npts)
{
  // int idx = threadIdx.x + blockIdx.x * blockDim.x;
  // if(idx < npts) {
  //   veloc[idx] += accel[idx] * dtover2;
  // }

// 1. Calculate thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    // 2. Calculate the data index for the vector (chunk of 4)
    int idx = tid * 4;

    // 3. Main Vectorized Path
    // Checks if a full chunk of 4 floats fits within bounds
    if (idx + 4 <= npts) {
        // Load data as 128-bit vectors
        // Note: We cast using 'tid' because we are accessing the array of float4s
        float4 v_vec = reinterpret_cast<float4*>(veloc)[tid];
        float4 a_vec = reinterpret_cast<const float4*>(accel)[tid];

        // Perform Math on all 4 components
        v_vec.x += a_vec.x * dtover2;
        v_vec.y += a_vec.y * dtover2;
        v_vec.z += a_vec.z * dtover2;
        v_vec.w += a_vec.w * dtover2;

        // Store back result
        reinterpret_cast<float4*>(veloc)[tid] = v_vec;
    }
    // 4. Tail Handling (Process remaining 1-3 floats)
    else if (idx < npts) {
        for (int k = 0; k < 4; ++k) {
            int current_idx = idx + k;
            if (current_idx < npts) {
                veloc[current_idx] += accel[current_idx] * dtover2;
            }
        }
    }
}

extern "C"
void update_velocity_device_(long* Mesh_pointer,realw *delta2ov2_f,const int *backward)
{
  TRACE("\tupdate_velocity_newmark");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = (mp->NGLOB_AB*NDIM + 3) / 4; // number of float4 elements

  int blocksize = BLOCKSIZE_KERNEL3;
  int nb = (size + blocksize - 1) / blocksize;
  realw dtover2 = *delta2ov2_f;

  if(!*backward) {
    kernel_UpdateVeloc<<<nb,blocksize,0,mp->compute_stream >>> (
      mp->d_veloc,mp->d_accel,dtover2,NDIM*mp->NGLOB_AB
    );
  }
  else {
    kernel_UpdateVeloc<<<nb,blocksize,0,mp->compute_stream >>> (
      mp->d_b_veloc,mp->d_b_accel,dtover2,NDIM*mp->NGLOB_AB
    );
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_b_cuda,
              KERNEL_3_B_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F) {
  TRACE("\tkernel_3_b_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltatover2 = *deltatover2_F;
  // updates only veloc at this point
  kernel_3_veloc_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                                      mp->d_accel,
                                                                      size,deltatover2);

  if (mp->simulation_type == 3) {
    realw b_deltatover2 = *b_deltatover2_F;
    kernel_3_veloc_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_b_veloc,
                                                                        mp->d_b_accel,
                                                                        size,b_deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel 3 b");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// acoustic domains

// KERNEL 3

/* ----------------------------------------------------------------------------------------------- */


__global__ void kernel_3_acoustic_cuda_device(field* potential_dot_acoustic,
                                                field* potential_dot_dot_acoustic,
                                                field* b_potential_dot_acoustic,
                                                field* b_potential_dot_dot_acoustic,
                                                int simulation_type,
                                                int size,
                                                realw deltatover2,
                                                realw b_deltatover2,
                                                realw* rmass_acoustic) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  realw rmass;
  field p_dot_dot;
  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    rmass = rmass_acoustic[id];
    // multiplies pressure with the inverse of the mass matrix
    p_dot_dot = rmass*potential_dot_dot_acoustic[id];
    potential_dot_dot_acoustic[id] = p_dot_dot;
    potential_dot_acoustic[id] += deltatover2*p_dot_dot;
    if (simulation_type==3) {
      p_dot_dot = rmass*b_potential_dot_dot_acoustic[id];
      b_potential_dot_dot_acoustic[id] = p_dot_dot;
      b_potential_dot_acoustic[id] += b_deltatover2*p_dot_dot;
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_acoustic_cuda,
              KERNEL_3_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                      realw* deltatover2_F,
                                      realw* b_deltatover2_F) {

TRACE("kernel_3_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltaover2 = *deltatover2_F;
  realw b_deltaover2 = *b_deltatover2_F;

  kernel_3_acoustic_cuda_device<<< grid, threads>>>(mp->d_potential_dot_acoustic,
                                                    mp->d_potential_dot_dot_acoustic,
                                                    mp->d_b_potential_dot_acoustic,
                                                    mp->d_b_potential_dot_dot_acoustic,
                                                    mp->simulation_type,
                                                    size,
                                                    deltaover2,
                                                    b_deltaover2,
                                                    mp->d_rmass_acoustic);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel 3 ");
#endif
} 

static __global__ void 
kernel_apply_dirichlet_mask_to_accel(realw* accel, const realw* dirichlet_mask, int nglob)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < nglob) {
    // Using registers instead of repeated memory access
    realw r = dirichlet_mask[idx];

    int base_idx = idx * NDIM;
    
    accel[base_idx]     *= r;
    accel[base_idx + 1] *= r;
    accel[base_idx + 2] *= r;
  }
}

extern "C"
void FC_FUNC_(apply_dirichlet_mask_to_accel,
  APPLY_DIRICHLET_MASK_TO_ACCEL)(long* Mesh_pointer)
{
  TRACE("\tapply_dirichlet_mask_to_accel");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int nb = (size + blocksize - 1) / blocksize;

  kernel_apply_dirichlet_mask_to_accel<<<nb,blocksize,0,mp->compute_stream>>>(
    mp->d_accel,mp->d_mask_dirichlet,size
  );
}