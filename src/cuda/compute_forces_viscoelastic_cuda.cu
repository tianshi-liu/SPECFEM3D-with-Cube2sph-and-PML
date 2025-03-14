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


#include "utils.cuh"

/* ----------------------------------------------------------------------------------------------- */

#ifdef USE_TEXTURES_FIELDS
realw_texture d_displ_tex;
realw_texture d_veloc_tex;
realw_texture d_accel_tex;
//backward/reconstructed
realw_texture d_b_displ_tex;
realw_texture d_b_veloc_tex;
realw_texture d_b_accel_tex;

//note: texture variables are implicitly static, and cannot be passed as arguments to cuda kernels;
//      thus, 1) we thus use if-statements (FORWARD_OR_ADJOINT) to determine from which texture to fetch from
//            2) we use templates
//      since if-statements are a bit slower as the variable is only known at runtime, we use option 2)

// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_displ(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_veloc(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_accel(int x);

// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_displ<1>(int x) { return tex1Dfetch(d_displ_tex, x); }
template<> __device__ float texfetch_veloc<1>(int x) { return tex1Dfetch(d_veloc_tex, x); }
template<> __device__ float texfetch_accel<1>(int x) { return tex1Dfetch(d_accel_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_displ<3>(int x) { return tex1Dfetch(d_b_displ_tex, x); }
template<> __device__ float texfetch_veloc<3>(int x) { return tex1Dfetch(d_b_veloc_tex, x); }
template<> __device__ float texfetch_accel<3>(int x) { return tex1Dfetch(d_b_accel_tex, x); }

#endif

#ifdef USE_TEXTURES_CONSTANTS
realw_texture d_hprime_xx_tex;
#endif



/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for elastic domains

/* ----------------------------------------------------------------------------------------------- */

// note:
// kernel_2 is split into several kernels:
//  - a kernel without attenuation and for isotropic media: Kernel_2_noatt_iso_impl()
//  - a kernel without attenuation and for isotropic media with coloring: Kernel_2_noatt_iso_col_impl()
//  - a kernel without attenuation and for isotropic media with gravity: Kernel_2_noatt_iso_grav_impl()
//  - a kernel without attenuation and for anisotropic media: Kernel_2_noatt_ani_impl()
//  - a kernel including attenuation: Kernel_2_att_impl()
//
// this should help with performance:
// the high number of registers needed for our kernels limits the occupancy; separation tries to reduce this.


// kernel without attenuation
//
// we use templates to distinguish between calls with forward or adjoint texture fields

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_impl(const int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        const int* d_irregular_element_number,
                        realw_p d_displ,
                        realw_p d_accel,
                        realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                        realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                        realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                        const realw xix_regular,const realw jacobian_regular,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw_const_p d_kappav,realw_const_p d_muv){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .false.
//  COMPUTE_AND_STORE_STRAIN  = .false.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;

  int I,J,K;
  int iglob,offset;
  int working_element, ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

// arithmetic intensity: ratio of number-of-arithmetic-operations / number-of-bytes-accessed-on-DRAM
//
// hand-counts on floating-point operations: counts addition/subtraction/multiplication/division
//                                           no counts for operations on indices in for-loops (compiler will likely unrool loops)
//
//                                           counts accesses to global memory, but no shared memory or register loads/stores
//                                           float has 4 bytes

// counts:
// 2 FLOP

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // limits thread ids to range [0,125-1]
  if (tx >= NGLL3) tx = NGLL3 - 1;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  ispec_irreg = d_irregular_element_number[working_element] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

// counts:
// + 8 FLOP
//
// ( 1 int + 2 float) * 128 threads = 1024 BYTE

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
  }

// counts:
// + 5 FLOP
//
// + 3 float * 125 threads = 1500 BYTE

  kappal = d_kappav[offset];
  mul = d_muv[offset];

// counts:
// + 0 FLOP
//
// 2 * 1 float * 128 threads = 1024 BYTE

  // local index
  K = (tx/NGLL2);
  J = ((tx-K*NGLL2)/NGLLX);
  I = (tx-K*NGLL2-J*NGLLX);

// counts:
// + 8 FLOP
//
// + 0 BYTE

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
// counts:
// + 0 FLOP
//
// 2 * 1 float * 25 threads = 200 BYTE

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes the spatial derivatives duxdxl ... depending on the regularity of the element
  get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                          &gammaxl,&gammayl,&gammazl,&jacobianl,I,J,K,tx,
                          &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,
                          &tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx,
                          &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                          d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,ispec_irreg,xix_regular,0);

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

// counts:
// + 22 FLOP
//
// + 0 BYTE

  // form dot product with test vector, symmetric form

  // 1. cut-plane xi
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  xixl,xiyl,xizl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,1);
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  etaxl,etayl,etazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,2);
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  gammaxl,gammayl,gammazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,3);
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

// counts:
// + 3 * 3 * 6 FLOP = 54 FLOP
// + 3 * 100 FLOP = 300 FLOP
//
// + 0 BYTE

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

// counts:
// + 3 * 2 FLOP = 6 FLOP
//
// + 3 float * 128 threads = 1536 BYTE

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

// counts:
// + 3 * 6 FLOP = 18 FLOP
//
// + 0 BYTE

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {
    atomicAdd(&d_accel[iglob*3], sum_terms1);
    atomicAdd(&d_accel[iglob*3+1], sum_terms2);
    atomicAdd(&d_accel[iglob*3+2], sum_terms3);
  }

// counts:
// + 8 FLOP
//
// + 3 float * 125 threads = 1500 BYTE


// counts:
// -----------------
// total of: 791 FLOP per thread
//           ~ 128 * 791 = 101248 FLOP per block
//
//           11392 BYTE DRAM accesses per block
//
// arithmetic intensity: 101120 FLOP / 11392 BYTES ~ 8.9 FLOP/BYTE
// -----------------
//
// nvprof: nvprof --metrics flops_sp ./xspecfem3D
//          -> 883146240 FLOPS (Single) floating-point operations for 20736 elements
//          -> 42590 FLOP per block
// arithmetic intensity: 42590 FLOP / 11392 BYTES ~ 3.74 FLOP/BYTE
//
// roofline model: Kepler K20x
// ---------------------------
//   for a Kepler K20x card, the peak single-precision performance is about 3.95 TFlop/s.
//   global memory access has a bandwidth of ~ 250 GB/s.
//
//   memory bandwidth: 250 GB/s
//   single-precision peak performance: 3.95 TFlop/s -> corner arithmetic intensity = 3950./250. ~ 15.8 flop/byte
//
//   elastic kernel has an arithmetic intensity of: hand-counts   ~ 8.9 flop/byte
//                                                  nvprof-counts ~ 42590./11392. flop/byte = 3.74 flop/byte
//
//   -> we can only achieve about: (hand-counts)   56% of the peak performance
//                                 (nvprof-counts) 24% of the peak performance -> 935.0 GFlop/s
//
// roofline model: Tesla K20c (Kepler architecture: http://www.nvidia.com/content/tesla/pdf/Tesla-KSeries-Overview-LR.pdf)
// ---------------------------
//   memory bandwidth: 208 GB/s
//   single-precision peak performance: 3.52 TFlop/s -> corner arithmetic intensity = 3520 / 208 ~ 16.9 flop/byte
//
//   we can only achieve about: (hand-counts)   52% of the peak performance
//                              (nvprof-counts) 22% of the peak performance -> 779.0 GFlop/s - measured: 647.3 GFlop/s


} // kernel_2_noatt_iso_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_strain_impl(int nb_blocks_to_compute,
                              const int* d_ibool,
                              const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                              const int d_iphase,
                              const int* d_irregular_element_number,
                              realw_p d_displ,
                              realw_p d_accel,
                              realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                              realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                              realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                              const realw xix_regular,const realw jacobian_regular,
                              realw_const_p d_hprime_xx,
                              realw_const_p d_hprimewgll_xx,
                              realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                              realw_const_p d_kappav,realw_const_p d_muv,
                              const int COMPUTE_AND_STORE_STRAIN,
                              realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                              realw_p epsilondev_xz,realw_p epsilondev_yz,
                              realw_p epsilon_trace_over_3,
                              const int SIMULATION_TYPE){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3 - 1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element,ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  ispec_irreg = d_irregular_element_number[working_element] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
  }



  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes the spatial derivatives duxdxl ... depending on the regularity of the element
  get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                          &gammaxl,&gammayl,&gammazl,&jacobianl,I,J,K,tx,
                          &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,
                          &tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx,
                          &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                          d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,ispec_irreg,xix_regular,0);

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // computes deviatoric strain for kernel calculations
  if (COMPUTE_AND_STORE_STRAIN) {
    // save deviatoric strain for Runge-Kutta scheme
    if (threadIdx.x < NGLL3) {
      realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
      // local storage: stresses at this current time step
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = duxdxl - templ; // epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = duydyl - templ; // epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = 0.5f * duxdyl_plus_duydxl; // epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = 0.5f * duzdxl_plus_duxdzl; // epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = 0.5f * duzdyl_plus_duydzl; //epsilondev_yz_loc;
      // kernel simulations
      if (SIMULATION_TYPE == 3){
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    } // threadIdx.x
  }

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology
  kappal = d_kappav[offset];
  mul = d_muv[offset];

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

  // form dot product with test vector, symmetric form

  // 1. cut-plane xi
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  xixl,xiyl,xizl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,1);
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  etaxl,etayl,etazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,2);
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  gammaxl,gammayl,gammazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,3);
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {
    atomicAdd(&d_accel[iglob*3], sum_terms1);
    atomicAdd(&d_accel[iglob*3+1], sum_terms2);
    atomicAdd(&d_accel[iglob*3+2], sum_terms3);
  } // threadIdx.x

} // kernel_2_noatt_iso_strain_impl()

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_col_impl(int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        const int use_mesh_coloring_gpu,
                        realw_p d_displ,
                        realw_p d_accel,
                        realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                        realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                        realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw_const_p d_kappav,realw_const_p d_muv,
                        const int COMPUTE_AND_STORE_STRAIN,
                        realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                        realw_p epsilondev_xz,realw_p epsilondev_yz,
                        realw_p epsilon_trace_over_3,
                        const int SIMULATION_TYPE){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .true.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3-1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if (use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  }
#endif
  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
  }

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xiyl = d_xiy[offset]; // all subsequent without to avoid over-use of texture for coalescent access
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)
                    -xiyl*(etaxl*gammazl-etazl*gammaxl)
                    +xizl*(etaxl*gammayl-etayl*gammaxl));

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // computes deviatoric strain for kernel calculations
  if (COMPUTE_AND_STORE_STRAIN) {
    // save deviatoric strain for Runge-Kutta scheme
    if (threadIdx.x < NGLL3) {
      realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
      // local storage: stresses at this current time step
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = duxdxl - templ; // epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = duydyl - templ; // epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = 0.5f * duxdyl_plus_duydxl; // epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = 0.5f * duzdxl_plus_duxdzl; // epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = 0.5f * duzdyl_plus_duydzl; //epsilondev_yz_loc;
      // kernel simulations
      if (SIMULATION_TYPE == 3){
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    } // threadIdx.x
  }

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology
  kappal = d_kappav[offset];
  mul = d_muv[offset];

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

  // form dot product with test vector, symmetric form
  // 1. cut-plane xi
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl); // sh_tempx1
    sh_tempy[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl); // sh_tempy1
    sh_tempz[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl); // sh_tempx2
    sh_tempy[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl); // sh_tempy2
    sh_tempz[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl); // sh_tempz2
  }
  __syncthreads();
  // 2. cut-plane eta
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {
    sh_tempx[tx] = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl); // sh_tempx3
    sh_tempy[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl); // sh_tempy3
    sh_tempz[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();
  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

  } // threadIdx.x

} // kernel_2_noatt_iso_col_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_grav_impl(int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        const int* d_irregular_element_number,
                        const int use_mesh_coloring_gpu,
                        realw_p d_displ,
                        realw_p d_accel,
                        realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                        realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                        realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                        const realw xix_regular,const realw jacobian_regular,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw_const_p d_kappav,realw_const_p d_muv,
                        const int COMPUTE_AND_STORE_STRAIN,
                        realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                        realw_p epsilondev_xz,realw_p epsilondev_yz,
                        realw_p epsilon_trace_over_3,
                        const int SIMULATION_TYPE,
                        const int gravity,
                        realw_const_p d_minus_g,
                        realw_const_p d_minus_deriv_gravity,
                        realw_const_p d_rhostore,
                        realw_const_p wgll_cube ){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .true.
//

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3-1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element,ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1 = 0.f;
  realw rho_s_H2 = 0.f;
  realw rho_s_H3 = 0.f;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if (use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  }
#endif
  ispec_irreg = d_irregular_element_number[working_element] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes the spatial derivatives duxdxl ... depending on the regularity of the element
  get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                          &gammaxl,&gammayl,&gammazl,&jacobianl,I,J,K,tx,
                          &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,
                          &tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx,
                          &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                          d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,ispec_irreg,xix_regular,0);

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // computes deviatoric strain for kernel calculations
  if (COMPUTE_AND_STORE_STRAIN) {
    // save deviatoric strain for Runge-Kutta scheme
    if (threadIdx.x < NGLL3) {
      realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
      // local storage: stresses at this current time step
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = duxdxl - templ; // epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = duydyl - templ; // epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = 0.5f * duxdyl_plus_duydxl; // epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = 0.5f * duzdxl_plus_duxdzl; // epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = 0.5f * duzdyl_plus_duydzl; //epsilondev_yz_loc;
      // kernel simulations
      if (SIMULATION_TYPE == 3){
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    } // threadIdx.x
  }

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology
  kappal = d_kappav[offset];
  mul = d_muv[offset];

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

  // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
  sigma_yx = sigma_xy;
  sigma_zx = sigma_xz;
  sigma_zy = sigma_yz;

  if (gravity ){
    //  computes non-symmetric terms for gravity
    compute_element_gravity(tx,working_element,&iglob,d_minus_g,d_minus_deriv_gravity,
                            d_rhostore,wgll_cube,jacobianl,
                            sh_tempx,sh_tempy,sh_tempz,
                            &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                            &rho_s_H1,&rho_s_H2,&rho_s_H3);
  }

  // form dot product with test vector, symmetric form

  // 1. cut-plane xi
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  xixl,xiyl,xizl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,1);
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  etaxl,etayl,etazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,2);
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  gammaxl,gammayl,gammazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,3);
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);


  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

  // adds gravity term
  sum_terms1 += rho_s_H1;
  sum_terms2 += rho_s_H2;
  sum_terms3 += rho_s_H3;

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

  } // threadIdx.x

} // kernel_2_noatt_iso_grav_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_ani_impl(int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        const int* d_irregular_element_number,
                        const int use_mesh_coloring_gpu,
                        realw_p d_displ,
                        realw_p d_accel,
                        realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                        realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                        realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                        const realw xix_regular,const realw jacobian_regular,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw_const_p d_kappav,realw_const_p d_muv,
                        const int COMPUTE_AND_STORE_STRAIN,
                        realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                        realw_p epsilondev_xz,realw_p epsilondev_yz,
                        realw_p epsilon_trace_over_3,
                        const int SIMULATION_TYPE,
                        const int ANISOTROPY,
                        realw* d_c11store,realw* d_c12store,realw* d_c13store,
                        realw* d_c14store,realw* d_c15store,realw* d_c16store,
                        realw* d_c22store,realw* d_c23store,realw* d_c24store,
                        realw* d_c25store,realw* d_c26store,realw* d_c33store,
                        realw* d_c34store,realw* d_c35store,realw* d_c36store,
                        realw* d_c44store,realw* d_c45store,realw* d_c46store,
                        realw* d_c55store,realw* d_c56store,realw* d_c66store,
                        const int gravity,
                        realw_const_p d_minus_g,
                        realw_const_p d_minus_deriv_gravity,
                        realw_const_p d_rhostore,
                        realw_const_p wgll_cube ){

// elastic compute kernel without attenuation for anisotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .true.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)


  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3-1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element,ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if (use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  }
#endif
  ispec_irreg = d_irregular_element_number[working_element] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes the spatial derivatives duxdxl ... depending on the regularity of the element
  get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                          &gammaxl,&gammayl,&gammazl,&jacobianl,I,J,K,tx,
                          &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,
                          &tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx,
                          &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                          d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,ispec_irreg,xix_regular,0);

  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // computes deviatoric strain for kernel calculations
  if (COMPUTE_AND_STORE_STRAIN) {
    realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
    // local storage: stresses at this current time step
    epsilondev_xx_loc = duxdxl - templ;
    epsilondev_yy_loc = duydyl - templ;
    epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl;
    epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl;
    epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl;

    if (SIMULATION_TYPE == 3){
      if (threadIdx.x < NGLL3) {
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    }
  }

  // full anisotropic case, stress calculations
  if (ANISOTROPY){
    c11 = d_c11store[offset];
    c12 = d_c12store[offset];
    c13 = d_c13store[offset];
    c14 = d_c14store[offset];
    c15 = d_c15store[offset];
    c16 = d_c16store[offset];
    c22 = d_c22store[offset];
    c23 = d_c23store[offset];
    c24 = d_c24store[offset];
    c25 = d_c25store[offset];
    c26 = d_c26store[offset];
    c33 = d_c33store[offset];
    c34 = d_c34store[offset];
    c35 = d_c35store[offset];
    c36 = d_c36store[offset];
    c44 = d_c44store[offset];
    c45 = d_c45store[offset];
    c46 = d_c46store[offset];
    c55 = d_c55store[offset];
    c56 = d_c56store[offset];
    c66 = d_c66store[offset];

    sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
               c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
    sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
               c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
    sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
               c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
    sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
               c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
    sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
               c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
    sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
               c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;

  }else{

    // isotropic case

    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
    lambdal = lambdalplus2mul - 2.0f * mul;

    // compute the six components of the stress tensor sigma
    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
    sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

    sigma_xy = mul*duxdyl_plus_duydxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
    sigma_yz = mul*duzdyl_plus_duydzl;
  }

  // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
  sigma_yx = sigma_xy;
  sigma_zx = sigma_xz;
  sigma_zy = sigma_yz;

  if (gravity ){
    //  computes non-symmetric terms for gravity
    compute_element_gravity(tx,working_element,&iglob,d_minus_g,d_minus_deriv_gravity,
                            d_rhostore,wgll_cube,jacobianl,
                            sh_tempx,sh_tempy,sh_tempz,
                            &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                            &rho_s_H1,&rho_s_H2,&rho_s_H3);
  }

  // form dot product with test vector, symmetric form

  // 1. cut-plane xi
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  xixl,xiyl,xizl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,1);
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  etaxl,etayl,etazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,2);
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                  gammaxl,gammayl,gammazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,3);
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

  // adds gravity term
  if (gravity ){
    sum_terms1 += rho_s_H1;
    sum_terms2 += rho_s_H2;
    sum_terms3 += rho_s_H3;
  }

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

    // save deviatoric strain for Runge-Kutta scheme
    if (COMPUTE_AND_STORE_STRAIN ){
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = epsilondev_yz_loc;
    }

  } // threadIdx.x

} // kernel_2_noatt_ani_impl()


/* ----------------------------------------------------------------------------------------------- */

// kernel with attenuation
//
// we use templates to distinguish between calls with forward or adjoint texture fields

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_att_impl(int nb_blocks_to_compute,
                  const int* d_ibool,
                  const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                  const int d_iphase,
                  const int* d_irregular_element_number,
                  const int use_mesh_coloring_gpu,
                  const realw d_deltat,
                  realw_p d_displ,
                  realw_const_p d_veloc,
                  realw_p d_accel,
                  realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                  realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                  realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                  const realw xix_regular,const realw jacobian_regular,
                  realw_const_p d_hprime_xx,
                  realw_const_p d_hprimewgll_xx,
                  realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                  realw_const_p d_kappav,realw_const_p d_muv,
                  realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                  realw_p epsilondev_xz,realw_p epsilondev_yz,
                  realw_p epsilon_trace_over_3,
                  const int SIMULATION_TYPE,
                  const int NSPEC,
                  realw_const_p factor_common,
                  realw_p R_xx,realw_p R_yy,realw_p R_xy,realw_p R_xz,realw_p R_yz,
                  realw_const_p alphaval,realw_const_p betaval,realw_const_p gammaval,
                  const int ANISOTROPY,
                  realw_const_p d_c11store,realw_const_p d_c12store,realw_const_p d_c13store,
                  realw_const_p d_c14store,realw_const_p d_c15store,realw_const_p d_c16store,
                  realw_const_p d_c22store,realw_const_p d_c23store,realw_const_p d_c24store,
                  realw_const_p d_c25store,realw_const_p d_c26store,realw_const_p d_c33store,
                  realw_const_p d_c34store,realw_const_p d_c35store,realw_const_p d_c36store,
                  realw_const_p d_c44store,realw_const_p d_c45store,realw_const_p d_c46store,
                  realw_const_p d_c55store,realw_const_p d_c56store,realw_const_p d_c66store,
                  const int gravity,
                  realw_const_p d_minus_g,
                  realw_const_p d_minus_deriv_gravity,
                  realw_const_p d_rhostore,
                  realw_const_p wgll_cube){


// elastic compute kernel with attenuation
// holds for: ATTENUATION = .true.
//            COMPUTE_AND_STORE_STRAIN = .true. (always true for attenuation)

  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  unsigned short int active;
  int iglob,offset;
  int working_element,ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw templ;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;
  realw Rxx_loc[N_SLS],Ryy_loc[N_SLS],Rxy_loc[N_SLS],Rxz_loc[N_SLS],Ryz_loc[N_SLS];


  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  // re-assigns shared array to decrease shared memory usage
  // note: this will re-use s_temp arrays from above to save shared memory
  realw* s_dummyx_loc_att = (realw*) sh_tempx;
  realw* s_dummyy_loc_att = (realw*) sh_tempy;
  realw* s_dummyz_loc_att = (realw*) sh_tempz;


  if (bx >= nb_blocks_to_compute) return;
  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  // because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active = (tx < NGLL3) ? 1:0;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (active ){

#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    //mesh coloring
    if (use_mesh_coloring_gpu ){
      working_element = bx;
    }else{
      // iphase-1 and working_element-1 for Fortran->C array conventions
      working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)]-1;
    }
#endif
    ispec_irreg = d_irregular_element_number[working_element] - 1;

    // local padded index
    offset = working_element*NGLL3_PADDED + tx;

    // global index
    iglob = d_ibool[offset] - 1;

    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,s_dummyx_loc,s_dummyy_loc,s_dummyz_loc);

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // use first order Taylor expansion of displacement for local storage of stresses
    // at this current time step, to fix attenuation in a consistent way
#ifdef USE_TEXTURES_FIELDS
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3);
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 1);
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 2);
#else
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * d_veloc[iglob*3];
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * d_veloc[iglob*3 + 1];
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * d_veloc[iglob*3 + 2];
#endif
  }// active

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active ){
   // attenuation
  // computes the spatial derivatives duxdxl ... depending on the regularity of the element
    get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                            &gammaxl,&gammayl,&gammazl,&jacobianl,I,J,K,tx,
                            &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,
                            &tempx3l,&tempy3l,&tempz3l,s_dummyx_loc_att,s_dummyy_loc_att,s_dummyz_loc_att,sh_hprime_xx,
                            &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,ispec_irreg,xix_regular,0);

    // attenuation
    // computes deviatoric strain attenuation and/or for kernel calculations
    templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333
    // local storage: stresses at this current time step
    epsilondev_xx_loc = duxdxl - templ;
    epsilondev_yy_loc = duydyl - templ;
    epsilondev_xy_loc = 0.5f * (duxdyl + duydxl);
    epsilondev_xz_loc = 0.5f * (duzdxl + duxdzl);
    epsilondev_yz_loc = 0.5f * (duzdyl + duydzl);

    if (SIMULATION_TYPE == 3) epsilon_trace_over_3[tx + working_element*NGLL3] = templ;

    // computes the spatial derivatives duxdxl ... depending on the regularity of the element
    get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                            &gammaxl,&gammayl,&gammazl,&jacobianl,I,J,K,tx,
                            &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,
                            &tempx3l,&tempy3l,&tempz3l,s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,sh_hprime_xx,
                            &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,ispec_irreg,xix_regular,1);

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

    // full anisotropic case, stress calculations
    if (ANISOTROPY){
      c11 = d_c11store[offset];
      c12 = d_c12store[offset];
      c13 = d_c13store[offset];
      c14 = d_c14store[offset];
      c15 = d_c15store[offset];
      c16 = d_c16store[offset];
      c22 = d_c22store[offset];
      c23 = d_c23store[offset];
      c24 = d_c24store[offset];
      c25 = d_c25store[offset];
      c26 = d_c26store[offset];
      c33 = d_c33store[offset];
      c34 = d_c34store[offset];
      c35 = d_c35store[offset];
      c36 = d_c36store[offset];
      c44 = d_c44store[offset];
      c45 = d_c45store[offset];
      c46 = d_c46store[offset];
      c55 = d_c55store[offset];
      c56 = d_c56store[offset];
      c66 = d_c66store[offset];

      sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
                 c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
      sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
                 c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
      sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
                 c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
      sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
                 c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
      sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
                 c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
      sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
                 c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;
    }else{

      // isotropic case

      // compute elements with an elastic isotropic rheology
      kappal = get_global_cr( &d_kappav[offset] );
      mul = get_global_cr( &d_muv[offset] );

      // attenuation

      lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
      lambdal = lambdalplus2mul - 2.0f * mul;

      // compute the six components of the stress tensor sigma
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
      sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

      sigma_xy = mul*duxdyl_plus_duydxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_yz = mul*duzdyl_plus_duydzl;
    }

    // attenuation
    // subtracts memory variables if attenuation
    compute_element_att_stress(tx,working_element,NSPEC,
                               R_xx,R_yy,R_xy,R_xz,R_yz,Rxx_loc,Ryy_loc,Rxy_loc,Rxz_loc,Ryz_loc,
                               &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz);

    // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;

    if (gravity ){
      //  computes non-symmetric terms for gravity
      compute_element_gravity(tx,working_element,&iglob,d_minus_g,d_minus_deriv_gravity,
                              d_rhostore,wgll_cube,jacobianl,
                              s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                              &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                              &rho_s_H1,&rho_s_H2,&rho_s_H3);
    }

  } // active

  // re-assigns sh_hprime_xx to load hprimewgll
  // note: the sync seems to be necessary, otherwise there is more jitter, not sure why...
  __syncthreads();
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprime_xx);
  }

// JC JC here we will need to add GPU support for the new C-PML routines
  __syncthreads();

  if (active ){
    // form dot product with test vector, symmetric form
    // 1. cut-plane xi
    __syncthreads();
    get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                    xixl,xiyl,xizl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,1);
    sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

    // 2. cut-plane eta
    __syncthreads();
    get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                    etaxl,etayl,etazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,2);
    sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

    // 3. cut-plane gamma
    __syncthreads();
    get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                    gammaxl,gammayl,gammazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,3);
    sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

    fac1 = d_wgllwgll_yz[K*NGLLX+J];
    fac2 = d_wgllwgll_xz[K*NGLLX+I];
    fac3 = d_wgllwgll_xy[J*NGLLX+I];

    sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
    sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
    sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

    // adds gravity term
    if (gravity ){
      sum_terms1 += rho_s_H1;
      sum_terms2 += rho_s_H2;
      sum_terms3 += rho_s_H3;
    }

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

// JC JC here we will need to add GPU support for the new C-PML routines

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }
    else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

    // attenuation
    // update memory variables based upon the Runge-Kutta scheme
    compute_element_att_memory(tx,working_element,NSPEC,
                               mul,
                               factor_common,alphaval,betaval,gammaval,
                               R_xx,R_yy,R_xy,R_xz,R_yz,Rxx_loc,Ryy_loc,Rxy_loc,Rxz_loc,Ryz_loc,
                               epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,
                               epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc);

    // save deviatoric strain for Runge-Kutta scheme
    // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
    epsilondev_xx[tx + working_element*NGLL3] = epsilondev_xx_loc;
    epsilondev_yy[tx + working_element*NGLL3] = epsilondev_yy_loc;
    epsilondev_xy[tx + working_element*NGLL3] = epsilondev_xy_loc;
    epsilondev_xz[tx + working_element*NGLL3] = epsilondev_xz_loc;
    epsilondev_yz[tx + working_element*NGLL3] = epsilondev_yz_loc;
  } // if (active)

// JC JC here we will need to add GPU support for the new C-PML routines

} // kernel_2_att_impl()


/* ----------------------------------------------------------------------------------------------- */

/*

// original kernel
// please leave it here for reference ...

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_att_org_impl(int nb_blocks_to_compute,
                  const int* d_ibool,
                  const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                  const int d_iphase,
                  const int use_mesh_coloring_gpu,
                  realw d_deltat,
                  realw_const_p d_displ,
                  realw_const_p d_veloc,
                  realw_p d_accel,
                  realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                  realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                  realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                  realw_const_p d_hprime_xx,
                  realw_const_p d_hprimewgll_xx,
                  realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                  realw_const_p d_kappav,realw_const_p d_muv,
                  realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                  realw_p epsilondev_xz,realw_p epsilondev_yz,
                  realw_p epsilon_trace_over_3,
                  const int SIMULATION_TYPE,
                  const int NSPEC,
                  realw_const_p one_minus_sum_beta,realw_const_p factor_common,
                  realw_p R_xx,realw_p R_yy,realw_p R_xy,realw_p R_xz,realw_p R_yz,
                  realw_const_p alphaval,realw_const_p betaval,realw_const_p gammaval,
                  const int ANISOTROPY,
                  realw_const_p d_c11store,realw_const_p d_c12store,realw_const_p d_c13store,
                  realw_const_p d_c14store,realw_const_p d_c15store,realw_const_p d_c16store,
                  realw_const_p d_c22store,realw_const_p d_c23store,realw_const_p d_c24store,
                  realw_const_p d_c25store,realw_const_p d_c26store,realw_const_p d_c33store,
                  realw_const_p d_c34store,realw_const_p d_c35store,realw_const_p d_c36store,
                  realw_const_p d_c44store,realw_const_p d_c45store,realw_const_p d_c46store,
                  realw_const_p d_c55store,realw_const_p d_c56store,realw_const_p d_c66store,
                  const int gravity,
                  realw_const_p d_minus_g,
                  realw_const_p d_minus_deriv_gravity,
                  realw_const_p d_rhostore,
                  realw_const_p wgll_cube){


// elastic compute kernel with attenuation
// holds for: ATTENUATION = .true.
//            COMPUTE_AND_STORE_STRAIN = .true. (always true for attenuation)

  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  unsigned short int active;
  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw templ;

  realw tempx1l_att,tempx2l_att,tempx3l_att,tempy1l_att,tempy2l_att,tempy3l_att,tempz1l_att,tempz2l_att,tempz3l_att;
  realw duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att,duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att;
  realw duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
#endif

  __shared__ realw sh_tempx1[NGLL3];
  __shared__ realw sh_tempx2[NGLL3];
  __shared__ realw sh_tempx3[NGLL3];

  __shared__ realw sh_tempy1[NGLL3];
  __shared__ realw sh_tempy2[NGLL3];
  __shared__ realw sh_tempy3[NGLL3];

  __shared__ realw sh_tempz1[NGLL3];
  __shared__ realw sh_tempz2[NGLL3];
  __shared__ realw sh_tempz3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  // re-assigns shared array to decrease shared memory usage
  // note: this will re-use s_temp arrays from above to save shared memory
  realw* s_dummyx_loc_att = (realw*) sh_tempx1;
  realw* s_dummyy_loc_att = (realw*) sh_tempx2;
  realw* s_dummyz_loc_att = (realw*) sh_tempx3;

  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  // because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (active ){

#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    //mesh coloring
    if (use_mesh_coloring_gpu ){
      working_element = bx;
    }else{
      // iphase-1 and working_element-1 for Fortran->C array conventions
      working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)]-1;
    }
#endif
    // local padded index
    offset = working_element*NGLL3_PADDED + tx;

    // global index
    iglob = d_ibool[offset] - 1;

    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,s_dummyx_loc,s_dummyy_loc,s_dummyz_loc);

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // use first order Taylor expansion of displacement for local storage of stresses
    // at this current time step, to fix attenuation in a consistent way
#ifdef USE_TEXTURES_FIELDS
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3);
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 1);
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * texfetch_veloc<FORWARD_OR_ADJOINT>(iglob*3 + 2);
#else
    s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * d_veloc[iglob*3];
    s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * d_veloc[iglob*3 + 1];
    s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * d_veloc[iglob*3 + 2];
#endif
  }// active

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active ){

#ifndef MANUALLY_UNROLLED_LOOPS

    tempx1l = 0.f;
    tempx2l = 0.f;
    tempx3l = 0.f;

    tempy1l = 0.f;
    tempy2l = 0.f;
    tempy3l = 0.f;

    tempz1l = 0.f;
    tempz2l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      //assumes that hprime_xx = hprime_yy = hprime_zz
      fac1 = sh_hprime_xx[l*NGLLX+I];
      tempx1l += s_dummyx_loc[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l += s_dummyy_loc[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l += s_dummyz_loc[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[l*NGLLX+J];
      tempx2l += s_dummyx_loc[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l += s_dummyy_loc[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l += s_dummyz_loc[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[l*NGLLX+K];
      tempx3l += s_dummyx_loc[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l += s_dummyy_loc[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l += s_dummyz_loc[l*NGLL2+J*NGLLX+I]*fac3;
    }

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    tempx1l_att = 0.f;
    tempx2l_att = 0.f;
    tempx3l_att = 0.f;

    tempy1l_att = 0.f;
    tempy2l_att = 0.f;
    tempy3l_att = 0.f;

    tempz1l_att = 0.f;
    tempz2l_att = 0.f;
    tempz3l_att = 0.f;

    for (l=0;l<NGLLX;l++) {
      fac1 = sh_hprime_xx[l*NGLLX+I];
      tempx1l_att += s_dummyx_loc_att[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l_att += s_dummyy_loc_att[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l_att += s_dummyz_loc_att[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[l*NGLLX+J];
      tempx2l_att += s_dummyx_loc_att[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l_att += s_dummyy_loc_att[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l_att += s_dummyz_loc_att[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[l*NGLLX+K];
      tempx3l_att += s_dummyx_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l_att += s_dummyy_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l_att += s_dummyz_loc_att[l*NGLL2+J*NGLLX+I]*fac3;
    }

#else
    tempx1l = s_dummyx_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempy1l = s_dummyy_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempz1l = s_dummyz_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempx2l = s_dummyx_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyx_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyx_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempy2l = s_dummyy_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyy_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyy_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempz2l = s_dummyz_loc[K*NGLL2+I]*d_hprime_xx[J]
            + s_dummyz_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
            + s_dummyz_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempx3l = s_dummyx_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyx_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyx_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyx_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyx_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempy3l = s_dummyy_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyy_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyy_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyy_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyy_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempz3l = s_dummyz_loc[J*NGLLX+I]*d_hprime_xx[K]
            + s_dummyz_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
            + s_dummyz_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
            + s_dummyz_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
            + s_dummyz_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    tempx1l_att = s_dummyx_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempy1l_att = s_dummyy_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempz1l_att = s_dummyz_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
                  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

    tempx2l_att = s_dummyx_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyx_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyx_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempy2l_att = s_dummyy_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyy_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyy_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempz2l_att = s_dummyz_loc_att[K*NGLL2+I]*d_hprime_xx[J]
                  + s_dummyz_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
                  + s_dummyz_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

    tempx3l_att = s_dummyx_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyx_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyx_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyx_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyx_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempy3l_att = s_dummyy_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyy_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyy_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyy_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyy_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

    tempz3l_att = s_dummyz_loc_att[J*NGLLX+I]*d_hprime_xx[K]
                  + s_dummyz_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
                  + s_dummyz_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
                  + s_dummyz_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
                  + s_dummyz_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];
#endif

    // compute derivatives of ux, uy and uz with respect to x, y and z
    xixl = get_global_cr( &d_xix[offset] );
    xiyl = get_global_cr( &d_xiy[offset] );
    xizl = get_global_cr( &d_xiz[offset] );
    etaxl = get_global_cr( &d_etax[offset] );
    etayl = get_global_cr( &d_etay[offset] );
    etazl = get_global_cr( &d_etaz[offset] );
    gammaxl = get_global_cr( &d_gammax[offset] );
    gammayl = get_global_cr( &d_gammay[offset] );
    gammazl = get_global_cr( &d_gammaz[offset] );

    jacobianl = 1.0f / (xixl*(etayl*gammazl-etazl*gammayl)
                       -xiyl*(etaxl*gammazl-etazl*gammaxl)
                       +xizl*(etaxl*gammayl-etayl*gammaxl));

    duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
    duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
    duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

    duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
    duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
    duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

    duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
    duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
    duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

// JC JC here we will need to add GPU support for the new C-PML routines

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

// JC JC here we will need to add GPU support for the new C-PML routines

    // attenuation
    // temporary variables used for fixing attenuation in a consistent way
    duxdxl_att = xixl*tempx1l_att + etaxl*tempx2l_att + gammaxl*tempx3l_att;
    duxdyl_att = xiyl*tempx1l_att + etayl*tempx2l_att + gammayl*tempx3l_att;
    duxdzl_att = xizl*tempx1l_att + etazl*tempx2l_att + gammazl*tempx3l_att;

    duydxl_att = xixl*tempy1l_att + etaxl*tempy2l_att + gammaxl*tempy3l_att;
    duydyl_att = xiyl*tempy1l_att + etayl*tempy2l_att + gammayl*tempy3l_att;
    duydzl_att = xizl*tempy1l_att + etazl*tempy2l_att + gammazl*tempy3l_att;

    duzdxl_att = xixl*tempz1l_att + etaxl*tempz2l_att + gammaxl*tempz3l_att;
    duzdyl_att = xiyl*tempz1l_att + etayl*tempz2l_att + gammayl*tempz3l_att;
    duzdzl_att = xizl*tempz1l_att + etazl*tempz2l_att + gammazl*tempz3l_att;

    // precompute some sums to save CPU time
    duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att;
    duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att;
    duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att;

    // attenuation
    // computes deviatoric strain attenuation and/or for kernel calculations
    templ = 0.33333333333333333333f * (duxdxl_att + duydyl_att + duzdzl_att); // 1./3. = 0.33333
    // local storage: stresses at this current time step
    epsilondev_xx_loc = duxdxl_att - templ;
    epsilondev_yy_loc = duydyl_att - templ;
    epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl_att;
    epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl_att;
    epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl_att;

    if (SIMULATION_TYPE == 3) {
      epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
    }

    // full anisotropic case, stress calculations
    if (ANISOTROPY){
      c11 = d_c11store[offset];
      c12 = d_c12store[offset];
      c13 = d_c13store[offset];
      c14 = d_c14store[offset];
      c15 = d_c15store[offset];
      c16 = d_c16store[offset];
      c22 = d_c22store[offset];
      c23 = d_c23store[offset];
      c24 = d_c24store[offset];
      c25 = d_c25store[offset];
      c26 = d_c26store[offset];
      c33 = d_c33store[offset];
      c34 = d_c34store[offset];
      c35 = d_c35store[offset];
      c36 = d_c36store[offset];
      c44 = d_c44store[offset];
      c45 = d_c45store[offset];
      c46 = d_c46store[offset];
      c55 = d_c55store[offset];
      c56 = d_c56store[offset];
      c66 = d_c66store[offset];

      sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
                 c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
      sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
                 c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
      sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
                 c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
      sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
                 c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
      sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
                 c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
      sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
                 c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;
    }else{

      // isotropic case

      // compute elements with an elastic isotropic rheology
      kappal = get_global_cr( &d_kappav[offset] );
      mul = get_global_cr( &d_muv[offset] );

      // attenuation
      // use unrelaxed parameters if attenuation
      mul  = mul * get_global_cr( &one_minus_sum_beta[tx+working_element*NGLL3] ); // (i,j,k,ispec)

      lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
      lambdal = lambdalplus2mul - 2.0f * mul;

      // compute the six components of the stress tensor sigma
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
      sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

      sigma_xy = mul*duxdyl_plus_duydxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_yz = mul*duzdyl_plus_duydzl;
    }

    // attenuation
    // subtracts memory variables if attenuation
    compute_element_att_stress(tx,working_element,NSPEC,
                               R_xx,R_yy,R_xy,R_xz,R_yz,
                               &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz);

    // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;

    if (gravity ){
      //  computes non-symmetric terms for gravity
      compute_element_gravity(tx,working_element,&iglob,d_minus_g,d_minus_deriv_gravity,
                              d_rhostore,wgll_cube,jacobianl,
                              s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                              &sigma_xx,&sigma_yy,&sigma_xz,&sigma_yz,
                              &rho_s_H1,&rho_s_H2,&rho_s_H3);
    }
  } // active

  //note: due to re-assignement of s_dummyx_loc_att,..,we need to sync before updating sh_tempx1...
  __syncthreads();

  if (active ){
    // form dot product with test vector, non-symmetric form
    sh_tempx1[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl);
    sh_tempy1[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl);
    sh_tempz1[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

    sh_tempx2[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl);
    sh_tempy2[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl);
    sh_tempz2[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

    sh_tempx3[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl);
    sh_tempy3[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl);
    sh_tempz3[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);
  }

  // re-assigns sh_hprime_xx to load hprimewgll
  // note: the sync seems to be necessary, otherwise there is more jitter, not sure why...
  __syncthreads();
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprime_xx);
  }

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

// JC JC here we will need to add GPU support for the new C-PML routines

  if (active ){

#ifndef MANUALLY_UNROLLED_LOOPS

    tempx1l = 0.f;
    tempy1l = 0.f;
    tempz1l = 0.f;

    tempx2l = 0.f;
    tempy2l = 0.f;
    tempz2l = 0.f;

    tempx3l = 0.f;
    tempy3l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      // assumes hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      fac1 = sh_hprime_xx[I*NGLLX+l]; //  d_hprimewgll_xx[I*NGLLX+l];
      tempx1l += sh_tempx1[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l += sh_tempy1[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l += sh_tempz1[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[J*NGLLX+l]; // d_hprimewgll_xx[J*NGLLX+l];
      tempx2l += sh_tempx2[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l += sh_tempy2[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l += sh_tempz2[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[K*NGLLX+l]; // d_hprimewgll_xx[K*NGLLX+l];
      tempx3l += sh_tempx3[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l += sh_tempy3[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l += sh_tempz3[l*NGLL2+J*NGLLX+I]*fac3;
    }
#else
    tempx1l = sh_tempx1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempx1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempx1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempx1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempx1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempy1l = sh_tempy1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempy1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempy1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempy1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempy1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempz1l = sh_tempz1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
            + sh_tempz1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
            + sh_tempz1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
            + sh_tempz1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
            + sh_tempz1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

    tempx2l = sh_tempx2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempx2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempx2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempx2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempx2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempy2l = sh_tempy2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempy2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempy2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempy2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempy2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempz2l = sh_tempz2[K*NGLL2+I]*d_hprimewgll_xx[J*NGLLX]
            + sh_tempz2[K*NGLL2+NGLLX+I]*d_hprimewgll_xx[J*NGLLX+1]
            + sh_tempz2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+2]
            + sh_tempz2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+3]
            + sh_tempz2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_xx[J*NGLLX+4];

    tempx3l = sh_tempx3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempx3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempx3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempx3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempx3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempy3l = sh_tempy3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempy3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempy3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempy3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempy3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];

    tempz3l = sh_tempz3[J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX]
            + sh_tempz3[NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+1]
            + sh_tempz3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+2]
            + sh_tempz3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+3]
            + sh_tempz3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_xx[K*NGLLX+4];
#endif

    fac1 = d_wgllwgll_yz[K*NGLLX+J];
    fac2 = d_wgllwgll_xz[K*NGLLX+I];
    fac3 = d_wgllwgll_xy[J*NGLLX+I];

    sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
    sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
    sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

    // adds gravity term
    if (gravity ){
      sum_terms1 += rho_s_H1;
      sum_terms2 += rho_s_H2;
      sum_terms3 += rho_s_H3;
    }

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

// JC JC here we will need to add GPU support for the new C-PML routines

#else // MESH_COLORING

    //mesh coloring
    if (use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }
    else {
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
    } // if (use_mesh_coloring_gpu)

#endif // MESH_COLORING

    // attenuation
    // update memory variables based upon the Runge-Kutta scheme
    compute_element_att_memory(tx,working_element,NSPEC,
                               d_muv,
                               factor_common,alphaval,betaval,gammaval,
                               R_xx,R_yy,R_xy,R_xz,R_yz,
                               epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,
                               epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc);

    // save deviatoric strain for Runge-Kutta scheme
    // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
    epsilondev_xx[tx + working_element*NGLL3] = epsilondev_xx_loc;
    epsilondev_yy[tx + working_element*NGLL3] = epsilondev_yy_loc;
    epsilondev_xy[tx + working_element*NGLL3] = epsilondev_xy_loc;
    epsilondev_xz[tx + working_element*NGLL3] = epsilondev_xz_loc;
    epsilondev_yz[tx + working_element*NGLL3] = epsilondev_yz_loc;
  } // if (active)

// JC JC here we will need to add GPU support for the new C-PML routines

} // kernel_2_att_impl()

*/

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_kelvinvoigt_impl(const int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        const int* d_irregular_element_number,
                        realw* d_kelvin_voigt_eta,
                        realw_p d_displ,
                        realw_const_p d_veloc,
                        realw_p d_accel,
                        realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                        realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                        realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                        const realw xix_regular,const realw jacobian_regular,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                        realw_const_p d_kappav,realw_const_p d_muv){

// elastic compute kernel without attenuation for isotropic elements with kelvin voigt damping aroung the fault
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .false.
//  COMPUTE_AND_STORE_STRAIN  = .false.
//  mp -> Kelvin_Voigt_damping = .true.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;

  int I,J,K;
  int iglob,offset;
  int working_element,ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw kelvin_voigt_eta;
  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

// arithmetic intensity: ratio of number-of-arithmetic-operations / number-of-bytes-accessed-on-DRAM
//
// hand-counts on floating-point operations: counts addition/subtraction/multiplication/division
//                                           no counts for operations on indices in for-loops (compiler will likely unrool loops)
//
//                                           counts accesses to global memory, but no shared memory or register loads/stores
//                                           float has 4 bytes

// counts:
// 2 FLOP
  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // limits thread ids to range [0,125-1]
  if (tx >= NGLL3) tx = NGLL3 - 1;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);

    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

// counts:
// + 0 FLOP
//
// 2 * 1 float * 25 threads = 200 BYTE

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  ispec_irreg = d_irregular_element_number[working_element] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;
  // fetch the value of kelvin_voigt eta
  kelvin_voigt_eta = d_kelvin_voigt_eta[bx] ;


// counts:
// + 7 FLOP
//
// + 2 float * 128 threads = 1024 BYTE

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL3 ){
    // copy displacement from global memory to shared memory
    load_shared_memory_displ_visco<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,d_veloc,kelvin_voigt_eta,sh_tempx,sh_tempy,sh_tempz);
  }

// counts:
// + 5 FLOP
//
// + 3 float * 125 threads = 1500 BYTE

  kappal = d_kappav[offset];
  mul = d_muv[offset];

// counts:
// + 0 FLOP
//
// + 2 * 1 float * 128 threads = 1024 BYTE


  // local index
  K = (tx/NGLL2);
  J = ((tx-K*NGLL2)/NGLLX);
  I = (tx-K*NGLL2-J*NGLLX);

// counts:
// + 8 FLOP
//
// + 0 BYTE

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes the spatial derivatives duxdxl ... depending on the regularity of the element
  get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
                          &gammaxl,&gammayl,&gammazl,&jacobianl,I,J,K,tx,
                          &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,
                          &tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx,
                          &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                          d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,ispec_irreg,xix_regular,0);


  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  sigma_xy = mul*duxdyl_plus_duydxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;
  sigma_yz = mul*duzdyl_plus_duydzl;

// counts:
// + 22 FLOP
//
// + 0 BYTE
  // form dot product with test vector, symmetric form

  // 1. cut-plane xi
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  xixl,xiyl,xizl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,1);
  sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 2. cut-plane eta
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  etaxl,etayl,etazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,2);
  sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                  gammaxl,gammayl,gammazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,3);
  sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

// counts:
// + 3 * 3 * 6 FLOP = 54 FLOP
// + 3 * 100 FLOP = 300 FLOP
//
// + 0 BYTE

  // gets double weights
  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

// counts:
// + 3 * 2 FLOP = 6 FLOP
//
// + 3 float * 128 threads = 1536 BYTE

  sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
  sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
  sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

// counts:
// + 3 * 6 FLOP = 18 FLOP
//
// + 0 BYTE

  // assembles acceleration array
  if (threadIdx.x < NGLL3) {
    atomicAdd(&d_accel[iglob*3], sum_terms1);
    atomicAdd(&d_accel[iglob*3+1], sum_terms2);
    atomicAdd(&d_accel[iglob*3+2], sum_terms3);
  }
} // kernel_2_noatt_iso_kelvinvoigt_impl()

/* ----------------------------------------------------------------------------------------------- */






void Kernel_2(int nb_blocks_to_compute,Mesh* mp,int d_iphase,realw d_deltat,
              int COMPUTE_AND_STORE_STRAIN,
              int ATTENUATION,
              int ANISOTROPY,
              int* d_ibool,
              realw* d_xix,realw* d_xiy,realw* d_xiz,
              realw* d_etax,realw* d_etay,realw* d_etaz,
              realw* d_gammax,realw* d_gammay,realw* d_gammaz,
              realw* d_kappav,
              realw* d_muv,
              realw* d_epsilondev_xx,realw* d_epsilondev_yy,realw* d_epsilondev_xy,
              realw* d_epsilondev_xz,realw* d_epsilondev_yz,
              realw* d_epsilon_trace_over_3,
              realw* d_factor_common,
              realw* d_R_xx,realw* d_R_yy,realw* d_R_xy,
              realw* d_R_xz,realw* d_R_yz,
              realw* d_b_epsilondev_xx,realw* d_b_epsilondev_yy,realw* d_b_epsilondev_xy,
              realw* d_b_epsilondev_xz,realw* d_b_epsilondev_yz,
              realw* d_b_epsilon_trace_over_3,
              realw* d_b_R_xx,realw* d_b_R_yy,realw* d_b_R_xy,
              realw* d_b_R_xz,realw* d_b_R_yz,
              realw* d_c11store,realw* d_c12store,realw* d_c13store,
              realw* d_c14store,realw* d_c15store,realw* d_c16store,
              realw* d_c22store,realw* d_c23store,realw* d_c24store,
              realw* d_c25store,realw* d_c26store,realw* d_c33store,
              realw* d_c34store,realw* d_c35store,realw* d_c36store,
              realw* d_c44store,realw* d_c45store,realw* d_c46store,
              realw* d_c55store,realw* d_c56store,realw* d_c66store,
              realw* d_rhostore){

  TRACE("\tKernel_2");

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel Kernel 2");
#endif

  // if the grid can handle the number of blocks, we let it be 1D
  // grid_2_x = nb_elem_color;
  // nb_elem_color is just how many blocks we are computing now

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  cudaEvent_t start,stop;
  if (CUDA_TIMING ){
    start_timing_cuda(&start,&stop);
  }

  // cuda kernel call
  if (ATTENUATION ){
    // compute kernels with attenuation
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
    Kernel_2_att_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                d_ibool,
                                                                mp->d_phase_ispec_inner_elastic,
                                                                mp->num_phase_ispec_elastic,
                                                                d_iphase,
                                                                mp->d_irregular_element_number,
                                                                mp->use_mesh_coloring_gpu,
                                                                d_deltat,
                                                                mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                d_xix, d_xiy, d_xiz,
                                                                d_etax, d_etay, d_etaz,
                                                                d_gammax, d_gammay, d_gammaz,
                                                                mp->xix_regular,mp->jacobian_regular,
                                                                mp->d_hprime_xx,
                                                                mp->d_hprimewgll_xx,
                                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                d_kappav, d_muv,
                                                                d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                d_epsilondev_xz,d_epsilondev_yz,
                                                                d_epsilon_trace_over_3,
                                                                mp->simulation_type,
                                                                mp->NSPEC_AB,
                                                                d_factor_common,
                                                                d_R_xx,d_R_yy,d_R_xy,d_R_xz,d_R_yz,
                                                                mp->d_alphaval,mp->d_betaval,mp->d_gammaval,
                                                                ANISOTROPY,
                                                                d_c11store,d_c12store,d_c13store,
                                                                d_c14store,d_c15store,d_c16store,
                                                                d_c22store,d_c23store,d_c24store,
                                                                d_c25store,d_c26store,d_c33store,
                                                                d_c34store,d_c35store,d_c36store,
                                                                d_c44store,d_c45store,d_c46store,
                                                                d_c55store,d_c56store,d_c66store,
                                                                mp->gravity,
                                                                mp->d_minus_g,
                                                                mp->d_minus_deriv_gravity,
                                                                d_rhostore,
                                                                mp->d_wgll_cube);

    if (mp->simulation_type == 3) {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      Kernel_2_att_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                   d_ibool,
                                                                   mp->d_phase_ispec_inner_elastic,
                                                                   mp->num_phase_ispec_elastic,
                                                                   d_iphase,
                                                                   mp->d_irregular_element_number,
                                                                   mp->use_mesh_coloring_gpu,
                                                                   d_deltat,
                                                                   mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                                                   d_xix, d_xiy, d_xiz,
                                                                   d_etax, d_etay, d_etaz,
                                                                   d_gammax, d_gammay, d_gammaz,
                                                                   mp->xix_regular,mp->jacobian_regular,
                                                                   mp->d_hprime_xx,
                                                                   mp->d_hprimewgll_xx,
                                                                   mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                   d_kappav, d_muv,
                                                                   d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                   d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                   d_b_epsilon_trace_over_3,
                                                                   mp->simulation_type,
                                                                   mp->NSPEC_AB,
                                                                   d_factor_common,
                                                                   d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                                                   mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                                                   ANISOTROPY,
                                                                   d_c11store,d_c12store,d_c13store,
                                                                   d_c14store,d_c15store,d_c16store,
                                                                   d_c22store,d_c23store,d_c24store,
                                                                   d_c25store,d_c26store,d_c33store,
                                                                   d_c34store,d_c35store,d_c36store,
                                                                   d_c44store,d_c45store,d_c46store,
                                                                   d_c55store,d_c56store,d_c66store,
                                                                   mp->gravity,
                                                                   mp->d_minus_g,
                                                                   mp->d_minus_deriv_gravity,
                                                                   d_rhostore,
                                                                   mp->d_wgll_cube);
    }
  }else{
    // compute kernels without attenuation
    if (ANISOTROPY ){
      // full anisotropy
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      Kernel_2_noatt_ani_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->d_irregular_element_number,
                                                                        mp->use_mesh_coloring_gpu,
                                                                        mp->d_displ,
                                                                        mp->d_accel,
                                                                        d_xix, d_xiy, d_xiz,
                                                                        d_etax, d_etay, d_etaz,
                                                                        d_gammax, d_gammay, d_gammaz,
                                                                        mp->xix_regular,mp->jacobian_regular,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                        d_kappav, d_muv,
                                                                        COMPUTE_AND_STORE_STRAIN,
                                                                        d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                        d_epsilondev_xz,d_epsilondev_yz,
                                                                        d_epsilon_trace_over_3,
                                                                        mp->simulation_type,
                                                                        ANISOTROPY,
                                                                        d_c11store,d_c12store,d_c13store,
                                                                        d_c14store,d_c15store,d_c16store,
                                                                        d_c22store,d_c23store,d_c24store,
                                                                        d_c25store,d_c26store,d_c33store,
                                                                        d_c34store,d_c35store,d_c36store,
                                                                        d_c44store,d_c45store,d_c46store,
                                                                        d_c55store,d_c56store,d_c66store,
                                                                        mp->gravity,
                                                                        mp->d_minus_g,
                                                                        mp->d_minus_deriv_gravity,
                                                                        d_rhostore,
                                                                        mp->d_wgll_cube );
      // backward/reconstructed wavefield
      if (mp->simulation_type == 3) {
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        Kernel_2_noatt_ani_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                           d_ibool,
                                                                           mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                           d_iphase,
                                                                           mp->d_irregular_element_number,
                                                                           mp->use_mesh_coloring_gpu,
                                                                           mp->d_b_displ,
                                                                           mp->d_b_accel,
                                                                           d_xix, d_xiy, d_xiz,
                                                                           d_etax, d_etay, d_etaz,
                                                                           d_gammax, d_gammay, d_gammaz,
                                                                           mp->xix_regular,mp->jacobian_regular,
                                                                           mp->d_hprime_xx,
                                                                           mp->d_hprimewgll_xx,
                                                                           mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                           d_kappav, d_muv,
                                                                           COMPUTE_AND_STORE_STRAIN,
                                                                           d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                           d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                           d_b_epsilon_trace_over_3,
                                                                           mp->simulation_type,
                                                                           ANISOTROPY,
                                                                           d_c11store,d_c12store,d_c13store,
                                                                           d_c14store,d_c15store,d_c16store,
                                                                           d_c22store,d_c23store,d_c24store,
                                                                           d_c25store,d_c26store,d_c33store,
                                                                           d_c34store,d_c35store,d_c36store,
                                                                           d_c44store,d_c45store,d_c46store,
                                                                           d_c55store,d_c56store,d_c66store,
                                                                           mp->gravity,
                                                                           mp->d_minus_g,
                                                                           mp->d_minus_deriv_gravity,
                                                                           d_rhostore,
                                                                           mp->d_wgll_cube );
      }
    }else{
      // isotropic case
      if (mp->gravity){
        // with gravity
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
        Kernel_2_noatt_iso_grav_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                          d_ibool,
                                                                          mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                          d_iphase,
                                                                          mp->d_irregular_element_number,
                                                                          mp->use_mesh_coloring_gpu,
                                                                          mp->d_displ,
                                                                          mp->d_accel,
                                                                          d_xix, d_xiy, d_xiz,
                                                                          d_etax, d_etay, d_etaz,
                                                                          d_gammax, d_gammay, d_gammaz,
                                                                          mp->xix_regular,mp->jacobian_regular,
                                                                          mp->d_hprime_xx,
                                                                          mp->d_hprimewgll_xx,
                                                                          mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                          d_kappav, d_muv,
                                                                          COMPUTE_AND_STORE_STRAIN,
                                                                          d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                          d_epsilondev_xz,d_epsilondev_yz,
                                                                          d_epsilon_trace_over_3,
                                                                          mp->simulation_type,
                                                                          mp->gravity,
                                                                          mp->d_minus_g,
                                                                          mp->d_minus_deriv_gravity,
                                                                          d_rhostore,
                                                                          mp->d_wgll_cube );

        // backward/reconstructed wavefield
        if (mp->simulation_type == 3) {
          // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
          Kernel_2_noatt_iso_grav_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                             d_ibool,
                                                                             mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                             d_iphase,
                                                                             mp->d_irregular_element_number,
                                                                             mp->use_mesh_coloring_gpu,
                                                                             mp->d_b_displ,
                                                                             mp->d_b_accel,
                                                                             d_xix, d_xiy, d_xiz,
                                                                             d_etax, d_etay, d_etaz,
                                                                             d_gammax, d_gammay, d_gammaz,
                                                                             mp->xix_regular,mp->jacobian_regular,
                                                                             mp->d_hprime_xx,
                                                                             mp->d_hprimewgll_xx,
                                                                             mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                             d_kappav, d_muv,
                                                                             COMPUTE_AND_STORE_STRAIN,
                                                                             d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                             d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                             d_b_epsilon_trace_over_3,
                                                                             mp->simulation_type,
                                                                             mp->gravity,
                                                                             mp->d_minus_g,
                                                                             mp->d_minus_deriv_gravity,
                                                                             d_rhostore,
                                                                             mp->d_wgll_cube );
        }
      }else{
        // without gravity
        if (mp->use_mesh_coloring_gpu) {
          // with mesh coloring
          // forward wavefields -> FORWARD_OR_ADJOINT == 1
          Kernel_2_noatt_iso_col_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                            d_ibool,
                                                                            mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                            d_iphase,
                                                                            mp->use_mesh_coloring_gpu,
                                                                            mp->d_displ,
                                                                            mp->d_accel,
                                                                            d_xix, d_xiy, d_xiz,
                                                                            d_etax, d_etay, d_etaz,
                                                                            d_gammax, d_gammay, d_gammaz,
                                                                            mp->d_hprime_xx,
                                                                            mp->d_hprimewgll_xx,
                                                                            mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                            d_kappav, d_muv,
                                                                            COMPUTE_AND_STORE_STRAIN,
                                                                            d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                            d_epsilondev_xz,d_epsilondev_yz,
                                                                            d_epsilon_trace_over_3,
                                                                            mp->simulation_type);

          // backward/reconstructed wavefield
          if (mp->simulation_type == 3) {
            // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
            Kernel_2_noatt_iso_col_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                               d_ibool,
                                                                               mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                               d_iphase,
                                                                               mp->use_mesh_coloring_gpu,
                                                                               mp->d_b_displ,
                                                                               mp->d_b_accel,
                                                                               d_xix, d_xiy, d_xiz,
                                                                               d_etax, d_etay, d_etaz,
                                                                               d_gammax, d_gammay, d_gammaz,
                                                                               mp->d_hprime_xx,
                                                                               mp->d_hprimewgll_xx,
                                                                               mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                               d_kappav, d_muv,
                                                                               COMPUTE_AND_STORE_STRAIN,
                                                                               d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                               d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                               d_b_epsilon_trace_over_3,
                                                                               mp->simulation_type);
          }
        }else{
          // without mesh coloring
          if (COMPUTE_AND_STORE_STRAIN ){
            // stores strains
            // forward wavefields -> FORWARD_OR_ADJOINT == 1
            Kernel_2_noatt_iso_strain_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                              d_ibool,
                                                                              mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                              d_iphase,
                                                                              mp->d_irregular_element_number,
                                                                              mp->d_displ,
                                                                              mp->d_accel,
                                                                              d_xix, d_xiy, d_xiz,
                                                                              d_etax, d_etay, d_etaz,
                                                                              d_gammax, d_gammay, d_gammaz,
                                                                              mp->xix_regular,mp->jacobian_regular,
                                                                              mp->d_hprime_xx,
                                                                              mp->d_hprimewgll_xx,
                                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                              d_kappav, d_muv,
                                                                              COMPUTE_AND_STORE_STRAIN,
                                                                              d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                                                              d_epsilondev_xz,d_epsilondev_yz,
                                                                              d_epsilon_trace_over_3,
                                                                              mp->simulation_type);

            // backward/reconstructed wavefield
            if (mp->simulation_type == 3) {
              // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
              Kernel_2_noatt_iso_strain_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                 d_ibool,
                                                                                 mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                 d_iphase,
                                                                                 mp->d_irregular_element_number,
                                                                                 mp->d_b_displ,
                                                                                 mp->d_b_accel,
                                                                                 d_xix, d_xiy, d_xiz,
                                                                                 d_etax, d_etay, d_etaz,
                                                                                 d_gammax, d_gammay, d_gammaz,
                                                                                 mp->xix_regular,mp->jacobian_regular,
                                                                                 mp->d_hprime_xx,
                                                                                 mp->d_hprimewgll_xx,
                                                                                 mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                 d_kappav, d_muv,
                                                                                 COMPUTE_AND_STORE_STRAIN,
                                                                                 d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                                                                 d_b_epsilondev_xz,d_b_epsilondev_yz,
                                                                                 d_b_epsilon_trace_over_3,
                                                                                 mp->simulation_type);
            }
          }else{
            if (mp->Kelvin_Voigt_damping) {
                         // Kelvin_Voigt_damping == true means there is fault in this partition
              Kernel_2_noatt_iso_kelvinvoigt_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                           d_ibool,
                                                                                           mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                           d_iphase,
                                                                                           mp->d_irregular_element_number,
                                                                                           mp->d_Kelvin_Voigt_eta,
                                                                                           mp->d_displ,
                                                                                           mp->d_veloc,
                                                                                           mp->d_accel,
                                                                                           d_xix, d_xiy, d_xiz,
                                                                                           d_etax, d_etay, d_etaz,
                                                                                           d_gammax, d_gammay, d_gammaz,
                                                                                           mp->xix_regular,mp->jacobian_regular,
                                                                                           mp->d_hprime_xx,
                                                                                           mp->d_hprimewgll_xx,
                                                                                           mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                           d_kappav, d_muv);
                       }
            else{
            // without storing strains
            // forward wavefields -> FORWARD_OR_ADJOINT == 1
            Kernel_2_noatt_iso_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                              d_ibool,
                                                                              mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                              d_iphase,
                                                                              mp->d_irregular_element_number,
                                                                              mp->d_displ,
                                                                              mp->d_accel,
                                                                              d_xix, d_xiy, d_xiz,
                                                                              d_etax, d_etay, d_etaz,
                                                                              d_gammax, d_gammay, d_gammaz,
                                                                              mp->xix_regular,mp->jacobian_regular,
                                                                              mp->d_hprime_xx,
                                                                              mp->d_hprimewgll_xx,
                                                                              mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                              d_kappav, d_muv);

            // backward/reconstructed wavefield
            if (mp->simulation_type == 3) {
              // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
              Kernel_2_noatt_iso_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                                 d_ibool,
                                                                                 mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                                 d_iphase,
                                                                                 mp->d_irregular_element_number,
                                                                                 mp->d_b_displ,
                                                                                 mp->d_b_accel,
                                                                                 d_xix, d_xiy, d_xiz,
                                                                                 d_etax, d_etay, d_etaz,
                                                                                 d_gammax, d_gammay, d_gammaz,
                                                                                 mp->xix_regular,mp->jacobian_regular,
                                                                                 mp->d_hprime_xx,
                                                                                 mp->d_hprimewgll_xx,
                                                                                 mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                 d_kappav, d_muv);
            }
            }
          } // COMPUTE_AND_STORE_STRAIN
        } // use_mesh_coloring_gpu
      } // gravity
    } // ANISOTROPY
  } // ATTENUATION

  // Cuda timing
  if (CUDA_TIMING ){
    if (ATTENUATION ){
      stop_timing_cuda(&start,&stop,"Kernel_2_att_impl");
    }else{
      if (ANISOTROPY ){
        stop_timing_cuda(&start,&stop,"Kernel_2_noatt_ani_impl");
      }else{
        if (mp->gravity ){
          stop_timing_cuda(&start,&stop,"Kernel_2_noatt_iso_grav_impl");
        }else{
          if (COMPUTE_AND_STORE_STRAIN ){
            stop_timing_cuda(&start,&stop,"Kernel_2_noatt_iso_strain_impl");
          }else{
            realw time;
            stop_timing_cuda(&start,&stop,"Kernel_2_noatt_iso_impl",&time);
            // time in seconds
            time = time / 1000.;
            // performance
            // see with: nvprof --metrics flops_sp ./xspecfem3D -> using 883146240 FLOPS (Single) floating-point operations
            // hand-counts: 89344 * number-of-blocks
            realw flops = 89344 * nb_blocks_to_compute;
            printf("  performance: %f GFlops/s\n", flops/time *(1./1000./1000./1000.));
          }
        }
      }
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("Kernel_2_impl");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_forces_viscoelastic_cuda,
              COMPUTE_FORCES_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphase,
                                                realw* deltat,
                                                int* nspec_outer_elastic,
                                                int* nspec_inner_elastic,
                                                int* COMPUTE_AND_STORE_STRAIN,
                                                int* ATTENUATION,
                                                int* ANISOTROPY) {

  TRACE("\tcompute_forces_viscoelastic_cuda");
  // EPIK_TRACER("compute_forces_viscoelastic_cuda");
  //printf("Running compute_forces\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int num_elements;

  if (*iphase == 1)
    num_elements = *nspec_outer_elastic;
  else
    num_elements = *nspec_inner_elastic;

  // checks if anything to do
  if (num_elements == 0) return;

  // mesh coloring
  if (mp->use_mesh_coloring_gpu ){
    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         elastic elements also start with outer than inner element ordering
    int nb_colors,nb_blocks_to_compute;
    int istart;
    int offset,offset_nonpadded,offset_nonpadded_att2;

    // sets up color loop
    if (*iphase == 1){
      // outer elements
      nb_colors = mp->num_colors_outer_elastic;
      istart = 0;

      // array offsets
      offset = 0;
      offset_nonpadded = 0;
      offset_nonpadded_att2 = 0;
    }else{
      // inner elements (start after outer elements)
      nb_colors = mp->num_colors_outer_elastic + mp->num_colors_inner_elastic;
      istart = mp->num_colors_outer_elastic;

      // array offsets
      offset = (*nspec_outer_elastic) * NGLL3_PADDED;
      offset_nonpadded = (*nspec_outer_elastic) * NGLL3;
      offset_nonpadded_att2 = (*nspec_outer_elastic) * NGLL3 * N_SLS;
    }

    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      nb_blocks_to_compute = mp->h_num_elem_colors_elastic[icolor];

      // checks
      //if (nb_blocks_to_compute <= 0){
      //  printf("error number of elastic color blocks: %d -- color = %d \n",nb_blocks_to_compute,icolor);
      //  exit(EXIT_FAILURE);
      //}

      Kernel_2(nb_blocks_to_compute,mp,*iphase,*deltat,
               *COMPUTE_AND_STORE_STRAIN,
               *ATTENUATION,*ANISOTROPY,
               mp->d_ibool + offset,
               mp->d_xix + offset,mp->d_xiy + offset,mp->d_xiz + offset,
               mp->d_etax + offset,mp->d_etay + offset,mp->d_etaz + offset,
               mp->d_gammax + offset,mp->d_gammay + offset,mp->d_gammaz + offset,
               mp->d_kappav + offset,
               mp->d_muv + offset,
               mp->d_epsilondev_xx + offset_nonpadded,mp->d_epsilondev_yy + offset_nonpadded,mp->d_epsilondev_xy + offset_nonpadded,
               mp->d_epsilondev_xz + offset_nonpadded,mp->d_epsilondev_yz + offset_nonpadded,
               mp->d_epsilon_trace_over_3 + offset_nonpadded,
               mp->d_factor_common + offset_nonpadded_att2,
               mp->d_R_xx + offset_nonpadded,mp->d_R_yy + offset_nonpadded,mp->d_R_xy + offset_nonpadded,
               mp->d_R_xz + offset_nonpadded,mp->d_R_yz + offset_nonpadded,
               mp->d_b_epsilondev_xx + offset_nonpadded,mp->d_b_epsilondev_yy + offset_nonpadded,mp->d_b_epsilondev_xy + offset_nonpadded,
               mp->d_b_epsilondev_xz + offset_nonpadded,mp->d_b_epsilondev_yz + offset_nonpadded,
               mp->d_b_epsilon_trace_over_3 + offset_nonpadded,
               mp->d_b_R_xx + offset_nonpadded,mp->d_b_R_yy + offset_nonpadded,mp->d_b_R_xy + offset_nonpadded,
               mp->d_b_R_xz + offset_nonpadded,mp->d_b_R_yz + offset_nonpadded,
               mp->d_c11store + offset,mp->d_c12store + offset,mp->d_c13store + offset,
               mp->d_c14store + offset,mp->d_c15store + offset,mp->d_c16store + offset,
               mp->d_c22store + offset,mp->d_c23store + offset,mp->d_c24store + offset,
               mp->d_c25store + offset,mp->d_c26store + offset,mp->d_c33store + offset,
               mp->d_c34store + offset,mp->d_c35store + offset,mp->d_c36store + offset,
               mp->d_c44store + offset,mp->d_c45store + offset,mp->d_c46store + offset,
               mp->d_c55store + offset,mp->d_c56store + offset,mp->d_c66store + offset,
               mp->d_rhostore + offset);

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL3;
      // for factor_common array
      offset_nonpadded_att2 += nb_blocks_to_compute * NGLL3 * N_SLS;

      //note: we use the same stream, so kernels are executed one after the other
      //      thus, there should be no need to synchronize in case we run on only 1 process to avoid race-conditions

    }

  }else{
    // no mesh coloring: uses atomic updates
    Kernel_2(num_elements,mp,*iphase,*deltat,
             *COMPUTE_AND_STORE_STRAIN,
             *ATTENUATION,*ANISOTROPY,
             mp->d_ibool,
             mp->d_xix,mp->d_xiy,mp->d_xiz,
             mp->d_etax,mp->d_etay,mp->d_etaz,
             mp->d_gammax,mp->d_gammay,mp->d_gammaz,
             mp->d_kappav,
             mp->d_muv,
             mp->d_epsilondev_xx,mp->d_epsilondev_yy,mp->d_epsilondev_xy,
             mp->d_epsilondev_xz,mp->d_epsilondev_yz,
             mp->d_epsilon_trace_over_3,
             mp->d_factor_common,
             mp->d_R_xx,mp->d_R_yy,mp->d_R_xy,
             mp->d_R_xz,mp->d_R_yz,
             mp->d_b_epsilondev_xx,mp->d_b_epsilondev_yy,mp->d_b_epsilondev_xy,
             mp->d_b_epsilondev_xz,mp->d_b_epsilondev_yz,
             mp->d_b_epsilon_trace_over_3,
             mp->d_b_R_xx,mp->d_b_R_yy,mp->d_b_R_xy,
             mp->d_b_R_xz,mp->d_b_R_yz,
             mp->d_c11store,mp->d_c12store,mp->d_c13store,
             mp->d_c14store,mp->d_c15store,mp->d_c16store,
             mp->d_c22store,mp->d_c23store,mp->d_c24store,
             mp->d_c25store,mp->d_c26store,mp->d_c33store,
             mp->d_c34store,mp->d_c35store,mp->d_c36store,
             mp->d_c44store,mp->d_c45store,mp->d_c46store,
             mp->d_c55store,mp->d_c56store,mp->d_c66store,
             mp->d_rhostore);
  }
}

#ifdef TEST

__device__ __forceinline__ void 
add_displ_discontinuity(int ispec,int tx,realw_const_p displ_wd,const int* ibool_wd,
                       const int *ispec_to_elem_wd,realw *ux,realw *uy,
                      realw* uz)
{
  int ispec_wd = ispec_to_elem_wd[ispec] - 1;
  if(ispec_wd > -1) {
    int iglob_wd = ibool_wd[ispec_wd * NGLL3_PADDED + tx] - 1;
    if(iglob_wd > -1) {
      *ux += displ_wd[iglob_wd*NDIM + 0];
      *uy += displ_wd[iglob_wd*NDIM + 1];
      *uz += displ_wd[iglob_wd*NDIM + 2];
    }
  }
}


__device__ __forceinline__ void
compute_stress(int ANISOTROPY,const realw* d_c11store,const realw* d_c12store,const realw* d_c13store,
              const realw* d_c14store,const realw* d_c15store,const realw* d_c16store,
              const realw* d_c22store,const realw* d_c23store,const realw* d_c24store,
              const realw* d_c25store,const realw* d_c26store,const realw* d_c33store,
              const realw* d_c34store,const realw* d_c35store,const realw* d_c36store,
              const realw* d_c44store,const realw* d_c45store,const realw* d_c46store,
              const realw* d_c55store,const realw* d_c56store,const realw* d_c66store,
              realw_const_p d_kappav,realw_const_p d_muv,
              realw duxdxl,realw duxdyl,realw duxdzl,realw duydxl,realw duydyl,
              realw duydzl,realw duzdxl,realw duzdyl,realw duzdzl,
              realw *sigma_xx,realw *sigma_yy,
              realw *sigma_zz, realw *sigma_xy,realw *sigma_xz,realw *sigma_yz,int offset
            )
{
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  // precompute some sums to save CPU time
  duxdxl_plus_duydyl = duxdxl + duydyl;
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duydyl_plus_duzdzl = duydyl + duzdzl;
  duxdyl_plus_duydxl = duxdyl + duydxl;
  duzdxl_plus_duxdzl = duzdxl + duxdzl;
  duzdyl_plus_duydzl = duzdyl + duydzl;
  

  // full anisotropic case, stress calculations
  if (ANISOTROPY){
    realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
    c11 = d_c11store[offset];
    c12 = d_c12store[offset];
    c13 = d_c13store[offset];
    c14 = d_c14store[offset];
    c15 = d_c15store[offset];
    c16 = d_c16store[offset];
    c22 = d_c22store[offset];
    c23 = d_c23store[offset];
    c24 = d_c24store[offset];
    c25 = d_c25store[offset];
    c26 = d_c26store[offset];
    c33 = d_c33store[offset];
    c34 = d_c34store[offset];
    c35 = d_c35store[offset];
    c36 = d_c36store[offset];
    c44 = d_c44store[offset];
    c45 = d_c45store[offset];
    c46 = d_c46store[offset];
    c55 = d_c55store[offset];
    c56 = d_c56store[offset];
    c66 = d_c66store[offset];

    *sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
               c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
    *sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
               c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
    *sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
               c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
    *sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
               c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
    *sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
               c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
    *sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
               c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;

  }else{

    // isotropic case
    realw lambdal,mul,lambdalplus2mul,kappal;
    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
    lambdal = lambdalplus2mul - 2.0f * mul;

    // compute the six components of the stress tensor sigma
    *sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
    *sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
    *sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

    *sigma_xy = mul*duxdyl_plus_duydxl;
    *sigma_xz = mul*duzdxl_plus_duxdzl;
    *sigma_yz = mul*duzdyl_plus_duydzl;
  }
}

__device__ __forceinline__ void 
compute_sigma_adepml(int ispec,int i,int j,int k,int ANISO,
                    const realw* d_c11store,const realw* d_c12store,const realw* d_c13store,
                    const realw* d_c14store,const realw* d_c15store,const realw* d_c16store,
                    const realw* d_c22store,const realw* d_c23store,const realw* d_c24store,
                    const realw* d_c25store,const realw* d_c26store,const realw* d_c33store,
                    const realw* d_c34store,const realw* d_c35store,const realw* d_c36store,
                    const realw* d_c44store,const realw* d_c45store,const realw* d_c46store,
                    const realw* d_c55store,const realw* d_c56store,const realw* d_c66store,
                    realw_const_p d_kappav,realw_const_p d_muv,int offset,
                    realw duxdxl, realw duydxl, realw duzdxl, 
                    realw duxdyl, realw duydyl, realw duzdyl,
                    realw duxdzl, realw duydzl, realw duzdzl,
                    realw_const_p rtrans, realw_const_p rtrans_inv,
                    realw_const_p pml_kappa, realw_p dsigma, const realw* Qu)
{
  const int idx = ispec * NGLL3 + k * NGLL2 + j * NGLLX + i;
  #define Rmat(i,j)  rtrans[(idx * NDIM + j-1) * NDIM + i-1]
  #define Rmatinv(i,j) rtrans_inv[(idx * NDIM + j-1) * NDIM + i-1]
  #define Qmat(i,j) Qu[(idx * NDIM + j-1) * NDIM + i-1]
  realw k1l,k2l,k3l,
      r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,
      ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,
      pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl,
      Qu1xl,Qu1yl,Qu1zl,Qu2xl,Qu2yl,Qu2zl,Qu3xl,Qu3yl,Qu3zl,
      e11,e12,e13,e22,e23,e33;
  r1xl = Rmat(1,1);
  r2xl = Rmat(2,1);
  r3xl = Rmat(3,1);
  r1yl = Rmat(1,2);
  r2yl = Rmat(2,2);
  r3yl = Rmat(3,2);
  r1zl = Rmat(1,3);
  r2zl = Rmat(2,3);
  r3zl = Rmat(3,3);

  ri1xl = Rmatinv(1,1);
  ri2xl = Rmatinv(2,1);
  ri3xl = Rmatinv(3,1);
  ri1yl = Rmatinv(1,2);
  ri2yl = Rmatinv(2,2);
  ri3yl = Rmatinv(3,2);
  ri1zl = Rmatinv(1,3);
  ri2zl = Rmatinv(2,3);
  ri3zl = Rmatinv(3,3);

  k1l = pml_kappa[idx * NDIM + 1 - 1];
  k2l = pml_kappa[idx * NDIM + 2 - 1];
  k3l = pml_kappa[idx * NDIM + 3 - 1];

  pmlxxl = r1xl*ri1xl/k1l + r2xl*ri2xl/k2l + r3xl*ri3xl/k3l;
  pmlxyl = r1xl*ri1yl/k1l + r2xl*ri2yl/k2l + r3xl*ri3yl/k3l;
  pmlxzl = r1xl*ri1zl/k1l + r2xl*ri2zl/k2l + r3xl*ri3zl/k3l;
  pmlyxl = r1yl*ri1xl/k1l + r2yl*ri2xl/k2l + r3yl*ri3xl/k3l;
  pmlyyl = r1yl*ri1yl/k1l + r2yl*ri2yl/k2l + r3yl*ri3yl/k3l;
  pmlyzl = r1yl*ri1zl/k1l + r2yl*ri2zl/k2l + r3yl*ri3zl/k3l;
  pmlzxl = r1zl*ri1xl/k1l + r2zl*ri2xl/k2l + r3zl*ri3xl/k3l;
  pmlzyl = r1zl*ri1yl/k1l + r2zl*ri2yl/k2l + r3zl*ri3yl/k3l;
  pmlzzl = r1zl*ri1zl/k1l + r2zl*ri2zl/k2l + r3zl*ri3zl/k3l;

  Qu1xl = Qmat(1,1);
  Qu2xl = Qmat(2,1);
  Qu3xl = Qmat(3,1);
  Qu1yl = Qmat(1,2);
  Qu2yl = Qmat(2,2);
  Qu3yl = Qmat(3,2);
  Qu1zl = Qmat(1,3);
  Qu2zl = Qmat(2,3);
  Qu3zl = Qmat(3,3);

  e11 = pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ 
               r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl;
  e22 = pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ 
               r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl;
  e33 = pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ 
               r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl;
  e12 = pmlxxl*duydxl+pmlxyl*duydyl+pmlxzl*duydzl+ 
               r1xl*Qu1yl+r2xl*Qu2yl+r3xl*Qu3yl+ 
        pmlyxl*duxdxl+pmlyyl*duxdyl+pmlyzl*duxdzl+ 
               r1yl*Qu1xl+r2yl*Qu2xl+r3yl*Qu3xl;
  e13 = pmlxxl*duzdxl+pmlxyl*duzdyl+pmlxzl*duzdzl+ 
               r1xl*Qu1zl+r2xl*Qu2zl+r3xl*Qu3zl+ 
        pmlzxl*duxdxl+pmlzyl*duxdyl+pmlzzl*duxdzl+ 
               r1zl*Qu1xl+r2zl*Qu2xl+r3zl*Qu3xl;
  e23 = pmlyxl*duzdxl+pmlyyl*duzdyl+pmlyzl*duzdzl+ 
               r1yl*Qu1zl+r2yl*Qu2zl+r3yl*Qu3zl+ 
        pmlzxl*duydxl+pmlzyl*duydyl+pmlzzl*duydzl+ 
               r1zl*Qu1yl+r2zl*Qu2yl+r3zl*Qu3yl;
  
  if(! ANISO) {
    realw kappal = d_kappav[offset], mul = d_muv[offset];
    realw lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
    realw lambdal = lambdalplus2mul - 2.0f * mul;
    dsigma[1-1] = lambdalplus2mul*e11+lambdal*(e22+e33); // xx
    dsigma[2-1] = mul*e12; // xy
    dsigma[3-1] = mul*e13;; // xz
    dsigma[4-1] = lambdalplus2mul*e22+lambdal*(e11+e33); // yy
    dsigma[5-1] = mul*e23;; // yz
    dsigma[6-1] = lambdalplus2mul*e33+lambdal*(e11+e22); //zz
  }
  else {
    realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
    c11 = d_c11store[offset];
    c12 = d_c12store[offset];
    c13 = d_c13store[offset];
    c14 = d_c14store[offset];
    c15 = d_c15store[offset];
    c16 = d_c16store[offset];
    c22 = d_c22store[offset];
    c23 = d_c23store[offset];
    c24 = d_c24store[offset];
    c25 = d_c25store[offset];
    c26 = d_c26store[offset];
    c33 = d_c33store[offset];
    c34 = d_c34store[offset];
    c35 = d_c35store[offset];
    c36 = d_c36store[offset];
    c44 = d_c44store[offset];
    c45 = d_c45store[offset];
    c46 = d_c46store[offset];
    c55 = d_c55store[offset];
    c56 = d_c56store[offset];
    c66 = d_c66store[offset];

    dsigma[0] = c11 * e11 + c16 * e12 + c12 * e22 +
        c15 * e13 + c14 * e23 + c13 * e33;
    dsigma[3] = c12 * e11 + c26 * e12 + c22 * e22 +
        c25 * e13 + c24 * e23 + c23 * e33;
    dsigma[5] = c13 * e11 + c36 * e12 + c23 * e22 +
        c35 * e13 + c34 * e23 + c33 * e33;
    dsigma[1] = c16 * e11 + c66 * e12 + c26 * e22 +
        c56 * e13 + c46 * e23 + c36 * e33;
    dsigma[2] = c15 * e11 + c56 * e12 + c25 * e22 +
        c55 * e13 + c45 * e23 + c35 * e33;
    dsigma[4] = c14 * e11 + c46 * e12 + c24 * e22 +
        c45 * e13 + c44 * e23 + c34 * e33;
  }

  #undef Qmat
  #undef Rmat
  #undef Rmatinv
}

__device__ __forceinline__ void 
compute_Qu_t_point(int ispec,int i,int j,int k,
                  realw duxdxl, realw duydxl, realw duzdxl, 
                  realw duxdyl, realw duydyl, realw duzdyl,
                  realw duxdzl, realw duydzl, realw duzdzl,
                  realw_const_p rtrans_inv,realw_const_p pml_d,
                  realw_p Qu_t)
{
  realw d1l,d2l,d3l,
        ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl;
  const int idx = ispec * NGLL3 + k * NGLL2 + j * NGLLX + i;
  
  #define Rmatinv(i,j) rtrans_inv[(idx * NDIM + j-1) * NDIM + i-1]
  #define Qmat(i,j) Qu_t[(idx * NDIM + j-1) * NDIM + i-1]
  ri1xl = Rmatinv(1,1);
  ri2xl = Rmatinv(2,1);
  ri3xl = Rmatinv(3,1);
  ri1yl = Rmatinv(1,2);
  ri2yl = Rmatinv(2,2);
  ri3yl = Rmatinv(3,2);
  ri1zl = Rmatinv(1,3);
  ri2zl = Rmatinv(2,3);
  ri3zl = Rmatinv(3,3);

  d1l = pml_d[idx*NDIM+1-1];
  d2l = pml_d[idx*NDIM+2-1];
  d3l = pml_d[idx*NDIM+3-1];

  Qmat(1,1) = -d1l*(ri1xl*duxdxl+ri1yl*duxdyl+ri1zl*duxdzl);
  Qmat(1,2) = -d1l*(ri1xl*duydxl+ri1yl*duydyl+ri1zl*duydzl);
  Qmat(1,3) = -d1l*(ri1xl*duzdxl+ri1yl*duzdyl+ri1zl*duzdzl);
  Qmat(2,1) = -d2l*(ri2xl*duxdxl+ri2yl*duxdyl+ri2zl*duxdzl);
  Qmat(2,2) = -d2l*(ri2xl*duydxl+ri2yl*duydyl+ri2zl*duydzl);
  Qmat(2,3) = -d2l*(ri2xl*duzdxl+ri2yl*duzdyl+ri2zl*duzdzl);
  Qmat(3,1) = -d3l*(ri3xl*duxdxl+ri3yl*duxdyl+ri3zl*duxdzl);
  Qmat(3,2) = -d3l*(ri3xl*duydxl+ri3yl*duydyl+ri3zl*duydzl);
  Qmat(3,3) = -d3l*(ri3xl*duzdxl+ri3yl*duzdyl+ri3zl*duzdzl);

  #undef Rmatinv
  #undef Qmat
}

__device__ __forceinline__ void 
compute_Qt_t(int ispec,int i,int j,int k,int tid,
             realw_const_p temp,realw fac1,
             realw fac2,realw fac3,
             realw_const_p r_trans, realw_const_p r_trans_inv, 
             realw_const_p pml_d,
             const int* ibool_CPML,
             realw_p Qt_t)
{
  realw d1l,d2l,d3l,
        r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,
        ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,
        pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl;
  const int idx = ispec * NGLL3 + k * NGLL2 + j * NGLLX + i;
  int iglob = ibool_CPML[idx] - 1;

  #define Rmat(i,j) r_trans[(idx * NDIM + j-1) * NDIM + i-1]
  #define Rmatinv(i,j) r_trans_inv[(idx * NDIM + j-1) * NDIM + i-1]
  r1xl = Rmat(1,1);
  r2xl = Rmat(2,1);
  r3xl = Rmat(3,1);
  r1yl = Rmat(1,2);
  r2yl = Rmat(2,2);
  r3yl = Rmat(3,2);
  r1zl = Rmat(1,3);
  r2zl = Rmat(2,3);
  r3zl = Rmat(3,3);

  ri1xl = Rmatinv(1,1);
  ri2xl = Rmatinv(2,1);
  ri3xl = Rmatinv(3,1);
  ri1yl = Rmatinv(1,2);
  ri2yl = Rmatinv(2,2);
  ri3yl = Rmatinv(3,2);
  ri1zl = Rmatinv(1,3);
  ri2zl = Rmatinv(2,3);
  ri3zl = Rmatinv(3,3);

  d1l = pml_d[1-1 + idx * NDIM];
  d2l = pml_d[2-1 + idx * NDIM];
  d3l = pml_d[3-1 + idx * NDIM];

  pmlxxl = r1xl*ri1xl;
  pmlxyl = r1xl*ri1yl;
  pmlxzl = r1xl*ri1zl;
  pmlyxl = r1yl*ri1xl;
  pmlyyl = r1yl*ri1yl;
  pmlyzl = r1yl*ri1zl;
  pmlzxl = r1zl*ri1xl;
  pmlzyl = r1zl*ri1yl;
  pmlzzl = r1zl*ri1zl;

  #define Qmat(a,b) Qt_t[(iglob * NDIM + b - 1) * NDIM + a-1 ]
  #define TMAT(a,b,c)  temp[((c-1) * NDIM + (b-1)) * 6 + a-1]

  realw val1,val2,val3;

  val1 = 
     d1l*(fac1*(TMAT(1,1,1)*pmlxxl+TMAT(1,2,1)*pmlxyl+TMAT(1,3,1)*pmlxzl+
                TMAT(2,1,1)*pmlyxl+TMAT(2,2,1)*pmlyyl+TMAT(2,3,1)*pmlyzl+
                TMAT(3,1,1)*pmlzxl+TMAT(3,2,1)*pmlzyl+TMAT(3,3,1)*pmlzzl) +
          fac2*(TMAT(1,1,2)*pmlxxl+TMAT(1,2,2)*pmlxyl+TMAT(1,3,2)*pmlxzl+
                TMAT(2,1,2)*pmlyxl+TMAT(2,2,2)*pmlyyl+TMAT(2,3,2)*pmlyzl+
                TMAT(3,1,2)*pmlzxl+TMAT(3,2,2)*pmlzyl+TMAT(3,3,2)*pmlzzl) +
          fac3*(TMAT(1,1,3)*pmlxxl+TMAT(1,2,3)*pmlxyl+TMAT(1,3,3)*pmlxzl+
                TMAT(2,1,3)*pmlyxl+TMAT(2,2,3)*pmlyyl+TMAT(2,3,3)*pmlyzl+
                TMAT(3,1,3)*pmlzxl+TMAT(3,2,3)*pmlzyl+TMAT(3,3,3)*pmlzzl));
  val2 = 
     d1l*(fac1*(TMAT(2,1,1)*pmlxxl+TMAT(2,2,1)*pmlxyl+TMAT(2,3,1)*pmlxzl+
                TMAT(4,1,1)*pmlyxl+TMAT(4,2,1)*pmlyyl+TMAT(4,3,1)*pmlyzl+
                TMAT(5,1,1)*pmlzxl+TMAT(5,2,1)*pmlzyl+TMAT(5,3,1)*pmlzzl) +
          fac2*(TMAT(2,1,2)*pmlxxl+TMAT(2,2,2)*pmlxyl+TMAT(2,3,2)*pmlxzl+
                TMAT(4,1,2)*pmlyxl+TMAT(4,2,2)*pmlyyl+TMAT(4,3,2)*pmlyzl+
                TMAT(5,1,2)*pmlzxl+TMAT(5,2,2)*pmlzyl+TMAT(5,3,2)*pmlzzl) +
          fac3*(TMAT(2,1,3)*pmlxxl+TMAT(2,2,3)*pmlxyl+TMAT(2,3,3)*pmlxzl+
                TMAT(4,1,3)*pmlyxl+TMAT(4,2,3)*pmlyyl+TMAT(4,3,3)*pmlyzl+
                TMAT(5,1,3)*pmlzxl+TMAT(5,2,3)*pmlzyl+TMAT(5,3,3)*pmlzzl));
  
  val3 = 
     d1l*(fac1*(TMAT(3,1,1)*pmlxxl+TMAT(3,2,1)*pmlxyl+TMAT(3,3,1)*pmlxzl+
                TMAT(5,1,1)*pmlyxl+TMAT(5,2,1)*pmlyyl+TMAT(5,3,1)*pmlyzl+
                TMAT(6,1,1)*pmlzxl+TMAT(6,2,1)*pmlzyl+TMAT(6,3,1)*pmlzzl) +
          fac2*(TMAT(3,1,2)*pmlxxl+TMAT(3,2,2)*pmlxyl+TMAT(3,3,2)*pmlxzl+
                TMAT(5,1,2)*pmlyxl+TMAT(5,2,2)*pmlyyl+TMAT(5,3,2)*pmlyzl+
                TMAT(6,1,2)*pmlzxl+TMAT(6,2,2)*pmlzyl+TMAT(6,3,2)*pmlzzl) +
          fac3*(TMAT(3,1,3)*pmlxxl+TMAT(3,2,3)*pmlxyl+TMAT(3,3,3)*pmlxzl+
                TMAT(5,1,3)*pmlyxl+TMAT(5,2,3)*pmlyyl+TMAT(5,3,3)*pmlyzl+
                TMAT(6,1,3)*pmlzxl+TMAT(6,2,3)*pmlzyl+TMAT(6,3,3)*pmlzzl));
  
  if(tid < NGLL3) {
    atomicAdd(&Qmat(1,1),val1);
    atomicAdd(&Qmat(1,2),val2);
    atomicAdd(&Qmat(1,3),val3);
  }


  pmlxxl = r2xl*ri2xl;
  pmlxyl = r2xl*ri2yl;
  pmlxzl = r2xl*ri2zl;
  pmlyxl = r2yl*ri2xl;
  pmlyyl = r2yl*ri2yl;
  pmlyzl = r2yl*ri2zl;
  pmlzxl = r2zl*ri2xl;
  pmlzyl = r2zl*ri2yl;
  pmlzzl = r2zl*ri2zl;

  val1 = 
     d2l*(fac1*(TMAT(1,1,1)*pmlxxl+TMAT(1,2,1)*pmlxyl+TMAT(1,3,1)*pmlxzl+
                TMAT(2,1,1)*pmlyxl+TMAT(2,2,1)*pmlyyl+TMAT(2,3,1)*pmlyzl+
                TMAT(3,1,1)*pmlzxl+TMAT(3,2,1)*pmlzyl+TMAT(3,3,1)*pmlzzl) +
          fac2*(TMAT(1,1,2)*pmlxxl+TMAT(1,2,2)*pmlxyl+TMAT(1,3,2)*pmlxzl+
                TMAT(2,1,2)*pmlyxl+TMAT(2,2,2)*pmlyyl+TMAT(2,3,2)*pmlyzl+
                TMAT(3,1,2)*pmlzxl+TMAT(3,2,2)*pmlzyl+TMAT(3,3,2)*pmlzzl) +
          fac3*(TMAT(1,1,3)*pmlxxl+TMAT(1,2,3)*pmlxyl+TMAT(1,3,3)*pmlxzl+
                TMAT(2,1,3)*pmlyxl+TMAT(2,2,3)*pmlyyl+TMAT(2,3,3)*pmlyzl+
                TMAT(3,1,3)*pmlzxl+TMAT(3,2,3)*pmlzyl+TMAT(3,3,3)*pmlzzl));
  val2 = 
     d2l*(fac1*(TMAT(2,1,1)*pmlxxl+TMAT(2,2,1)*pmlxyl+TMAT(2,3,1)*pmlxzl+
                TMAT(4,1,1)*pmlyxl+TMAT(4,2,1)*pmlyyl+TMAT(4,3,1)*pmlyzl+
                TMAT(5,1,1)*pmlzxl+TMAT(5,2,1)*pmlzyl+TMAT(5,3,1)*pmlzzl) +
          fac2*(TMAT(2,1,2)*pmlxxl+TMAT(2,2,2)*pmlxyl+TMAT(2,3,2)*pmlxzl+
                TMAT(4,1,2)*pmlyxl+TMAT(4,2,2)*pmlyyl+TMAT(4,3,2)*pmlyzl+
                TMAT(5,1,2)*pmlzxl+TMAT(5,2,2)*pmlzyl+TMAT(5,3,2)*pmlzzl) +
          fac3*(TMAT(2,1,3)*pmlxxl+TMAT(2,2,3)*pmlxyl+TMAT(2,3,3)*pmlxzl+
                TMAT(4,1,3)*pmlyxl+TMAT(4,2,3)*pmlyyl+TMAT(4,3,3)*pmlyzl+
                TMAT(5,1,3)*pmlzxl+TMAT(5,2,3)*pmlzyl+TMAT(5,3,3)*pmlzzl));
  val3 = 
     d2l*(fac1*(TMAT(3,1,1)*pmlxxl+TMAT(3,2,1)*pmlxyl+TMAT(3,3,1)*pmlxzl+
                TMAT(5,1,1)*pmlyxl+TMAT(5,2,1)*pmlyyl+TMAT(5,3,1)*pmlyzl+
                TMAT(6,1,1)*pmlzxl+TMAT(6,2,1)*pmlzyl+TMAT(6,3,1)*pmlzzl) +
          fac2*(TMAT(3,1,2)*pmlxxl+TMAT(3,2,2)*pmlxyl+TMAT(3,3,2)*pmlxzl+
                TMAT(5,1,2)*pmlyxl+TMAT(5,2,2)*pmlyyl+TMAT(5,3,2)*pmlyzl+
                TMAT(6,1,2)*pmlzxl+TMAT(6,2,2)*pmlzyl+TMAT(6,3,2)*pmlzzl) +
          fac3*(TMAT(3,1,3)*pmlxxl+TMAT(3,2,3)*pmlxyl+TMAT(3,3,3)*pmlxzl+
                TMAT(5,1,3)*pmlyxl+TMAT(5,2,3)*pmlyyl+TMAT(5,3,3)*pmlyzl+
                TMAT(6,1,3)*pmlzxl+TMAT(6,2,3)*pmlzyl+TMAT(6,3,3)*pmlzzl));
  if(tid < NGLL3) {
    atomicAdd(&Qmat(2,1),val1);
    atomicAdd(&Qmat(2,2),val2);
    atomicAdd(&Qmat(2,3),val3);
  }

  pmlxxl = r3xl*ri3xl;
  pmlxyl = r3xl*ri3yl;
  pmlxzl = r3xl*ri3zl;
  pmlyxl = r3yl*ri3xl;
  pmlyyl = r3yl*ri3yl;
  pmlyzl = r3yl*ri3zl;
  pmlzxl = r3zl*ri3xl;
  pmlzyl = r3zl*ri3yl;
  pmlzzl = r3zl*ri3zl;
  
  val1 =
     d3l*(fac1*(TMAT(1,1,1)*pmlxxl+TMAT(1,2,1)*pmlxyl+TMAT(1,3,1)*pmlxzl+
                TMAT(2,1,1)*pmlyxl+TMAT(2,2,1)*pmlyyl+TMAT(2,3,1)*pmlyzl+
                TMAT(3,1,1)*pmlzxl+TMAT(3,2,1)*pmlzyl+TMAT(3,3,1)*pmlzzl) +
          fac2*(TMAT(1,1,2)*pmlxxl+TMAT(1,2,2)*pmlxyl+TMAT(1,3,2)*pmlxzl+
                TMAT(2,1,2)*pmlyxl+TMAT(2,2,2)*pmlyyl+TMAT(2,3,2)*pmlyzl+
                TMAT(3,1,2)*pmlzxl+TMAT(3,2,2)*pmlzyl+TMAT(3,3,2)*pmlzzl) +
          fac3*(TMAT(1,1,3)*pmlxxl+TMAT(1,2,3)*pmlxyl+TMAT(1,3,3)*pmlxzl+
                TMAT(2,1,3)*pmlyxl+TMAT(2,2,3)*pmlyyl+TMAT(2,3,3)*pmlyzl+
                TMAT(3,1,3)*pmlzxl+TMAT(3,2,3)*pmlzyl+TMAT(3,3,3)*pmlzzl));
  val2 = 
     d3l*(fac1*(TMAT(2,1,1)*pmlxxl+TMAT(2,2,1)*pmlxyl+TMAT(2,3,1)*pmlxzl+
                TMAT(4,1,1)*pmlyxl+TMAT(4,2,1)*pmlyyl+TMAT(4,3,1)*pmlyzl+
                TMAT(5,1,1)*pmlzxl+TMAT(5,2,1)*pmlzyl+TMAT(5,3,1)*pmlzzl) +
          fac2*(TMAT(2,1,2)*pmlxxl+TMAT(2,2,2)*pmlxyl+TMAT(2,3,2)*pmlxzl+
                TMAT(4,1,2)*pmlyxl+TMAT(4,2,2)*pmlyyl+TMAT(4,3,2)*pmlyzl+
                TMAT(5,1,2)*pmlzxl+TMAT(5,2,2)*pmlzyl+TMAT(5,3,2)*pmlzzl) +
          fac3*(TMAT(2,1,3)*pmlxxl+TMAT(2,2,3)*pmlxyl+TMAT(2,3,3)*pmlxzl+
                TMAT(4,1,3)*pmlyxl+TMAT(4,2,3)*pmlyyl+TMAT(4,3,3)*pmlyzl+
                TMAT(5,1,3)*pmlzxl+TMAT(5,2,3)*pmlzyl+TMAT(5,3,3)*pmlzzl));
  val3 = 
     d3l*(fac1*(TMAT(3,1,1)*pmlxxl+TMAT(3,2,1)*pmlxyl+TMAT(3,3,1)*pmlxzl+
                TMAT(5,1,1)*pmlyxl+TMAT(5,2,1)*pmlyyl+TMAT(5,3,1)*pmlyzl+
                TMAT(6,1,1)*pmlzxl+TMAT(6,2,1)*pmlzyl+TMAT(6,3,1)*pmlzzl) +
          fac2*(TMAT(3,1,2)*pmlxxl+TMAT(3,2,2)*pmlxyl+TMAT(3,3,2)*pmlxzl+
                TMAT(5,1,2)*pmlyxl+TMAT(5,2,2)*pmlyyl+TMAT(5,3,2)*pmlyzl+
                TMAT(6,1,2)*pmlzxl+TMAT(6,2,2)*pmlzyl+TMAT(6,3,2)*pmlzzl) +
          fac3*(TMAT(3,1,3)*pmlxxl+TMAT(3,2,3)*pmlxyl+TMAT(3,3,3)*pmlxzl+
                TMAT(5,1,3)*pmlyxl+TMAT(5,2,3)*pmlyyl+TMAT(5,3,3)*pmlyzl+
                TMAT(6,1,3)*pmlzxl+TMAT(6,2,3)*pmlzyl+TMAT(6,3,3)*pmlzzl));
  if(tid < NGLL3) {
    atomicAdd(&Qmat(3,1),val1);
    atomicAdd(&Qmat(3,2),val2);
    atomicAdd(&Qmat(3,3),val3);
  }

  #undef Qmat 
  #undef Rmat 
  #undef Rmatinv
  #undef TMAT
}

__device__ __forceinline__ void  
mxm_3op(int i,int j,int k,realw_const_p hprimewgll,
        realw_const_p ux,realw_const_p uy,
        realw_const_p uz,realw_p tempx,realw_p tempy,
        realw_p tempz)
{
  realw sx{},sy{},sz{};
  #pragma unroll
  for(int l = 0; l < NGLLX; l ++) {
    sx += hprimewgll[i*NGLLX+l] * ux[k*NGLL2+j*NGLLX+l];
    sy += hprimewgll[j*NGLLX+l] * uy[k*NGLL2+l*NGLLX+i];
    sz += hprimewgll[k*NGLLX+l] * uz[l*NGLL2+j*NGLLX+i];
  }
  *tempx = sx;
  *tempy = sy;
  *tempz = sz;
}

__device__ __forceinline__ void 
compute_accel_adepml(int ispec,int i,int j,int k,realw_const_p temp,
                    realw fac1,realw fac2,realw fac3,realw_const_p r_trans, 
                    realw_const_p r_trans_inv, realw_const_p pml_kappa,
                    realw_const_p pml_d, realw_p ax,realw_p ay, realw_p az)
{
  realw k1l,k2l,k3l,
          r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,
          pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl;

  const int idx = ispec * NGLL3 + k * NGLL2 + j * NGLLX + i;
  #define Rmat(i,j) r_trans[(idx * NDIM + j-1) * NDIM + i-1]
  #define Rmatinv(i,j) r_trans_inv[(idx * NDIM + j-1) * NDIM + i-1]
  #define TMAT(a,b,c)  temp[((c-1) * NDIM + (b-1)) * 6 + a-1]

  r1xl = Rmat(1,1);
  r1yl = Rmat(1,2);
  r1zl = Rmat(1,3);
  r2xl = Rmat(2,1);
  r2yl = Rmat(2,2);
  r2zl = Rmat(2,3);
  r3xl = Rmat(3,1);
  r3yl = Rmat(3,2);
  r3zl = Rmat(3,3);

  ri1xl = Rmatinv(1,1);
  ri1yl = Rmatinv(1,2);
  ri1zl = Rmatinv(1,3);
  ri2xl = Rmatinv(2,1);
  ri2yl = Rmatinv(2,2);
  ri2zl = Rmatinv(2,3);
  ri3xl = Rmatinv(3,1);
  ri3yl = Rmatinv(3,2);
  ri3zl = Rmatinv(3,3);

  k1l = pml_kappa[1-1 + idx * NDIM];
  k2l = pml_kappa[2-1 + idx * NDIM];
  k3l = pml_kappa[3-1 + idx * NDIM];

  pmlxxl = r1xl*ri1xl/k1l + r2xl*ri2xl/k2l + r3xl*ri3xl/k3l;
  pmlxyl = r1xl*ri1yl/k1l + r2xl*ri2yl/k2l + r3xl*ri3yl/k3l;
  pmlxzl = r1xl*ri1zl/k1l + r2xl*ri2zl/k2l + r3xl*ri3zl/k3l;
  pmlyxl = r1yl*ri1xl/k1l + r2yl*ri2xl/k2l + r3yl*ri3xl/k3l;
  pmlyyl = r1yl*ri1yl/k1l + r2yl*ri2yl/k2l + r3yl*ri3yl/k3l;
  pmlyzl = r1yl*ri1zl/k1l + r2yl*ri2zl/k2l + r3yl*ri3zl/k3l;
  pmlzxl = r1zl*ri1xl/k1l + r2zl*ri2xl/k2l + r3zl*ri3xl/k3l;
  pmlzyl = r1zl*ri1yl/k1l + r2zl*ri2yl/k2l + r3zl*ri3yl/k3l;
  pmlzzl = r1zl*ri1zl/k1l + r2zl*ri2zl/k2l + r3zl*ri3zl/k3l;

  //! 1-xx 2-xy 3-xz 4-yy 5-yz 6-zz
  *ax = 
          fac1*(TMAT(1,1,1)*pmlxxl+TMAT(1,2,1)*pmlxyl+TMAT(1,3,1)*pmlxzl+
                TMAT(2,1,1)*pmlyxl+TMAT(2,2,1)*pmlyyl+TMAT(2,3,1)*pmlyzl+
                TMAT(3,1,1)*pmlzxl+TMAT(3,2,1)*pmlzyl+TMAT(3,3,1)*pmlzzl) +
          fac2*(TMAT(1,1,2)*pmlxxl+TMAT(1,2,2)*pmlxyl+TMAT(1,3,2)*pmlxzl+
                TMAT(2,1,2)*pmlyxl+TMAT(2,2,2)*pmlyyl+TMAT(2,3,2)*pmlyzl+
                TMAT(3,1,2)*pmlzxl+TMAT(3,2,2)*pmlzyl+TMAT(3,3,2)*pmlzzl) +
          fac3*(TMAT(1,1,3)*pmlxxl+TMAT(1,2,3)*pmlxyl+TMAT(1,3,3)*pmlxzl+
                TMAT(2,1,3)*pmlyxl+TMAT(2,2,3)*pmlyyl+TMAT(2,3,3)*pmlyzl+
                TMAT(3,1,3)*pmlzxl+TMAT(3,2,3)*pmlzyl+TMAT(3,3,3)*pmlzzl);
  *ay =
          fac1*(TMAT(2,1,1)*pmlxxl+TMAT(2,2,1)*pmlxyl+TMAT(2,3,1)*pmlxzl+
                TMAT(4,1,1)*pmlyxl+TMAT(4,2,1)*pmlyyl+TMAT(4,3,1)*pmlyzl+
                TMAT(5,1,1)*pmlzxl+TMAT(5,2,1)*pmlzyl+TMAT(5,3,1)*pmlzzl) +
          fac2*(TMAT(2,1,2)*pmlxxl+TMAT(2,2,2)*pmlxyl+TMAT(2,3,2)*pmlxzl+
                TMAT(4,1,2)*pmlyxl+TMAT(4,2,2)*pmlyyl+TMAT(4,3,2)*pmlyzl+
                TMAT(5,1,2)*pmlzxl+TMAT(5,2,2)*pmlzyl+TMAT(5,3,2)*pmlzzl) +
          fac3*(TMAT(2,1,3)*pmlxxl+TMAT(2,2,3)*pmlxyl+TMAT(2,3,3)*pmlxzl+
                TMAT(4,1,3)*pmlyxl+TMAT(4,2,3)*pmlyyl+TMAT(4,3,3)*pmlyzl+
                TMAT(5,1,3)*pmlzxl+TMAT(5,2,3)*pmlzyl+TMAT(5,3,3)*pmlzzl);
  *az = 
          fac1*(TMAT(3,1,1)*pmlxxl+TMAT(3,2,1)*pmlxyl+TMAT(3,3,1)*pmlxzl+
                TMAT(5,1,1)*pmlyxl+TMAT(5,2,1)*pmlyyl+TMAT(5,3,1)*pmlyzl+
                TMAT(6,1,1)*pmlzxl+TMAT(6,2,1)*pmlzyl+TMAT(6,3,1)*pmlzzl) +
          fac2*(TMAT(3,1,2)*pmlxxl+TMAT(3,2,2)*pmlxyl+TMAT(3,3,2)*pmlxzl+
                TMAT(5,1,2)*pmlyxl+TMAT(5,2,2)*pmlyyl+TMAT(5,3,2)*pmlyzl+
                TMAT(6,1,2)*pmlzxl+TMAT(6,2,2)*pmlzyl+TMAT(6,3,2)*pmlzzl) +
          fac3*(TMAT(3,1,3)*pmlxxl+TMAT(3,2,3)*pmlxyl+TMAT(3,3,3)*pmlxzl+
                TMAT(5,1,3)*pmlyxl+TMAT(5,2,3)*pmlyyl+TMAT(5,3,3)*pmlyzl+
                TMAT(6,1,3)*pmlzxl+TMAT(6,2,3)*pmlzyl+TMAT(6,3,3)*pmlzzl);
    
  #undef Rmat 
  #undef Rmatinv
  #undef TMAT
}

__device__ __forceinline__ void 
add_pml_physical_contribution(int ispec_pml, int igll, const int *phy_spec,
                              const int *phy_ijk,const int *ibool_CPML,
                              realw_const_p phy_norm,realw_const_p rtrans,
                              realw_const_p rtrans_inv,realw_const_p pml_d,
                              realw_const_p phy_jaco2Dw,realw_const_p sigma_ade,
                              realw_p Qt_t)
{

  realw sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yz,sigma_zz, 
  d1l,d2l,d3l,r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,
  ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,
  weight, tx, ty, tz;

  for(int iside = 0; iside < 2 * NDIM; iside ++) {
    // get current loc
    int iface = phy_spec[ispec_pml * 6 + iside] - 1;
    if(iface < 0) continue;

    const int idx1 = (iface*NGLL2+igll)*NDIM;
    int i = phy_ijk[idx1+0] - 1;
    int j = phy_ijk[idx1+1] - 1;
    int k = phy_ijk[idx1+2] - 1;

    const int idx = ispec_pml * NGLL3 + k * NGLL2 + j * NGLLX + i;
    
    int iglob = ibool_CPML[idx] - 1;
    realw nx = phy_norm[idx1+0];
    realw ny = phy_norm[idx1+1];
    realw nz = phy_norm[idx1+2];

    #define Rmat(i,j) rtrans[(idx * NDIM + j-1) * NDIM + i-1]
    #define Rmatinv(i,j) rtrans_inv[(idx * NDIM + j-1) * NDIM + i-1]
    #define Qmat(a,b) Qt_t[(iglob * NDIM + b - 1) * NDIM + a-1 ]

    r1xl = Rmat(1,1);
    r1yl = Rmat(1,2);
    r1zl = Rmat(1,3);
    r2xl = Rmat(2,1);
    r2yl = Rmat(2,2);
    r2zl = Rmat(2,3);
    r3xl = Rmat(3,1);
    r3yl = Rmat(3,2);
    r3zl = Rmat(3,3);

    ri1xl = Rmatinv(1,1);
    ri1yl = Rmatinv(1,2);
    ri1zl = Rmatinv(1,3);
    ri2xl = Rmatinv(2,1);
    ri2yl = Rmatinv(2,2);
    ri2zl = Rmatinv(2,3);
    ri3xl = Rmatinv(3,1);
    ri3yl = Rmatinv(3,2);
    ri3zl = Rmatinv(3,3);

    d1l = pml_d[1-1+idx*NDIM];
    d2l = pml_d[2-1+idx*NDIM];
    d3l = pml_d[3-1+idx*NDIM];

    sigma_xx = sigma_ade[1-1];
    sigma_xy = sigma_ade[2-1];
    sigma_xz = sigma_ade[3-1];
    sigma_yy = sigma_ade[4-1];
    sigma_yz = sigma_ade[5-1];
    sigma_zz = sigma_ade[6-1];

    weight = phy_jaco2Dw[iface * NGLL2 + igll];
    tx = d1l*(ri1xl*nx+ri1yl*ny+ri1zl*nz)*
                        (r1xl*sigma_xx+r1yl*sigma_xy+r1zl*sigma_xz);
    ty = d1l*(ri1xl*nx+ri1yl*ny+ri1zl*nz)*
                        (r1xl*sigma_xy+r1yl*sigma_yy+r1zl*sigma_yz);
    tz = d1l*(ri1xl*nx+ri1yl*ny+ri1zl*nz)*
                        (r1xl*sigma_xz+r1yl*sigma_yz+r1zl*sigma_zz);
    atomicAdd(&Qmat(1,1),-tx * weight);
    atomicAdd(&Qmat(1,2),-ty * weight);
    atomicAdd(&Qmat(1,3),-tz * weight);

    tx = d2l*(ri2xl*nx+ri2yl*ny+ri2zl*nz)*
                        (r2xl*sigma_xx+r2yl*sigma_xy+r2zl*sigma_xz);
    ty = d2l*(ri2xl*nx+ri2yl*ny+ri2zl*nz)*
                        (r2xl*sigma_xy+r2yl*sigma_yy+r2zl*sigma_yz);
    tz = d2l*(ri2xl*nx+ri2yl*ny+ri2zl*nz)*
                        (r2xl*sigma_xz+r2yl*sigma_yz+r2zl*sigma_zz);
    atomicAdd(&Qmat(2,1),-tx * weight);
    atomicAdd(&Qmat(2,2),-ty * weight);
    atomicAdd(&Qmat(2,3),-tz * weight);


    tx = d3l*(ri3xl*nx+ri3yl*ny+ri3zl*nz)*
                        (r3xl*sigma_xx+r3yl*sigma_xy+r3zl*sigma_xz);
    ty = d3l*(ri3xl*nx+ri3yl*ny+ri3zl*nz)*
                        (r3xl*sigma_xy+r3yl*sigma_yy+r3zl*sigma_yz);
    tz = d3l*(ri3xl*nx+ri3yl*ny+ri3zl*nz)*
                        (r3xl*sigma_xz+r3yl*sigma_yz+r3zl*sigma_zz);
    atomicAdd(&Qmat(3,1),-tx * weight);
    atomicAdd(&Qmat(3,2),-ty * weight);
    atomicAdd(&Qmat(3,3),-tz * weight);

    #undef Qmat 
    #undef Rmat
    #undef Rmatinv
  }
}

// nqdu added ADE-PML kernels
__global__ void 
kernel_forces_adepml(int nb_blocks_to_compute,
                      const int* d_ibool,
                     const int* d_phase_ispec_inner_elastic,
                     const int num_phase_ispec_elastic,
                     const int d_iphase,
                     const int* d_irregular_element_number,
                     const realw* d_displ,
                      realw_p d_accel,
                      const int ANISOTROPY,
                      const realw* d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                      const realw* d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                      const realw* d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                      const realw xix_regular,const realw jacobian_regular,
                      realw_const_p d_hprime_xx,
                      realw_const_p d_hprimewgll_xx,
                      realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,
                      realw_const_p d_kappav,realw_const_p d_muv,
                      const realw* d_c11store, const realw* d_c12store,const realw* d_c13store,
                      const realw* d_c14store, const realw* d_c15store,const realw* d_c16store,
                      const realw* d_c22store, const realw* d_c23store,const realw* d_c24store,
                      const realw* d_c25store, const realw* d_c26store,const realw* d_c33store,
                      const realw* d_c34store, const realw* d_c35store,const realw* d_c36store,
                      const realw* d_c44store, const realw* d_c45store,const realw* d_c46store,
                      const realw* d_c55store, const realw* d_c56store,const realw* d_c66store,
                      int is_wavediscon,realw_const_p displ_wd, 
                      const int* ispec_to_elem_wd, const int *ibool_wd,
                      const int *is_CPML,
                      const int *spec_to_CPML, const int *ibool_CPML,
                      realw_const_p rtrans,realw_const_p rtrans_inv,
                      realw_const_p pml_kappa,realw_const_p pml_d,
                      const int *phy_ijk, const int *phy_spec,
                      realw_const_p phy_norm,
                      realw_const_p phy_jaco2Dw,
                      const realw* Qu,realw_p Qu_t,realw_p Qt_t)
{
  // elastic compute kernel without attenuation for anisotropic elements, with PML

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3-1;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int iglob,offset;
  int working_element,ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;

  realw fac1 = 0,fac2 = 0,fac3 = 0;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  // realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;

  realw sum_terms1,sum_terms2,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  ispec_irreg = d_irregular_element_number[working_element] - 1;

  // local padded index
  offset = working_element*NGLL3_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;


  // loads hprime's into shared memory
  if (threadIdx.x < NGLL3) {
    // copy hprime from global memory to shared memory
    if(threadIdx.x < NGLL2) {
      sh_hprime_xx[tx] = d_hprime_xx[tx];
      sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];
    }
    realw ux = d_displ[iglob*NDIM];
    realw uy = d_displ[iglob*NDIM+1];
    realw uz = d_displ[iglob*NDIM+2];
    
    if(is_wavediscon)
      add_displ_discontinuity(working_element,tx,displ_wd,
                              ibool_wd,ispec_to_elem_wd,
                                &ux,&uy,&uz);

    sh_tempx[tx] = ux;
    sh_tempy[tx] = uy;
    sh_tempz[tx] = uz;
  }
  __syncthreads();

  // computes the spatial derivatives duxdxl ... depending on the regularity of the element
  get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,
    &gammaxl,&gammayl,&gammazl,&jacobianl,I,J,K,tx,
    &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,
    &tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx,
    &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
    d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,ispec_irreg,xix_regular,0);
  
  if(!is_CPML[working_element]) {
    // compute stress
    compute_stress(ANISOTROPY,d_c11store,d_c12store,d_c13store,d_c14store,d_c15store,
                  d_c16store,d_c22store,d_c23store, d_c24store, d_c25store,d_c26store,
                  d_c33store, d_c34store, d_c35store, d_c36store, d_c44store, 
                  d_c45store, d_c46store, d_c55store, d_c56store, d_c66store,
                  d_kappav,d_muv,duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl,
                  &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz,offset);
    

    // 1. cut-plane xi
    __syncthreads();
    get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                    xixl,xiyl,xizl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,1);
    sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

    // 2. cut-plane eta
    __syncthreads();
    get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                    etaxl,etayl,etazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,2);
    sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

    // 3. cut-plane gamma
    __syncthreads();
    get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_xy,sigma_xz,sigma_xz,sigma_yy,sigma_yz,sigma_yz,sigma_zz,
                    gammaxl,gammayl,gammazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,3);
    sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

    // gets double weights
    if(threadIdx.x < NGLL3) {
      fac1 = d_wgllwgll_yz[K*NGLLX+J];
      fac2 = d_wgllwgll_xz[K*NGLLX+I];
      fac3 = d_wgllwgll_xy[J*NGLLX+I];
    }

    sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
    sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
    sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

    // assembles acceleration array
    if (threadIdx.x < NGLL3) {
      
      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);
      
    } // threadIdx.x
  }
  else {
    realw dsigma[6];
    realw temp_adepml[3][3][6]{}, newtemp_adepml[3][3][6]{};
    int ispec_pml = spec_to_CPML[working_element] - 1;
    compute_sigma_adepml(ispec_pml,I,J,K,ANISOTROPY,d_c11store,d_c12store,d_c13store,
                        d_c14store,d_c15store,d_c16store,d_c22store,d_c23store,
                        d_c24store, d_c25store,d_c26store,d_c33store, d_c34store, 
                        d_c35store, d_c36store, d_c44store, d_c45store, d_c46store,
                        d_c55store, d_c56store, d_c66store,d_kappav,d_muv,offset,
                        duxdxl,duydxl,duzdxl,duxdyl,duydyl,duzdyl,duxdzl,duydzl,duzdzl,
                        rtrans,rtrans_inv,pml_kappa,dsigma,Qu
                      );
    if(threadIdx.x < NGLL3) {
      compute_Qu_t_point(ispec_pml,I,J,K,duxdxl, duydxl, duzdxl,
                          duxdyl, duydyl, duzdyl,
                          duxdzl, duydzl, duzdzl,rtrans_inv,
                        pml_d,Qu_t);
      #define SETTEMP(a,b) temp_adepml[b-1][a-1][l]
      for(int l = 0; l < 6; l ++) {
        SETTEMP(1,1) = dsigma[l] * xixl * jacobianl;
        SETTEMP(2,1) = dsigma[l] * xiyl * jacobianl;
        SETTEMP(3,1) = dsigma[l] * xizl * jacobianl;
        SETTEMP(1,2) = dsigma[l] * etaxl * jacobianl;
        SETTEMP(2,2) = dsigma[l] * etayl * jacobianl;
        SETTEMP(3,2) = dsigma[l] * etazl * jacobianl;
        SETTEMP(1,3) = dsigma[l] * gammaxl * jacobianl;
        SETTEMP(2,3) = dsigma[l] * gammayl * jacobianl;
        SETTEMP(3,3) = dsigma[l] * gammazl * jacobianl;
      }
      #undef SETTEMP
    }

    // compute new_temp
    for(int q = 0; q < 3; q ++) {
      for(int r = 0; r <  6; r ++) {
          if(threadIdx.x < NGLL3) {
            sh_tempx[tx] = temp_adepml[0][q][r];
            sh_tempy[tx] = temp_adepml[1][q][r];
            sh_tempz[tx] = temp_adepml[2][q][r];
          }
          __syncthreads();
          mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,
                  &newtemp_adepml[0][q][r],&newtemp_adepml[1][q][r],
                &newtemp_adepml[2][q][r]);
          __syncthreads();
      }
    }

    // update Qt_t
    compute_Qt_t(ispec_pml,I,J,K,threadIdx.x,&newtemp_adepml[0][0][0],
                  d_wgllwgll_yz[K*NGLLX+J],d_wgllwgll_xz[K*NGLLX+I],
                d_wgllwgll_xy[J*NGLLX+I],rtrans,rtrans_inv,
              pml_d,ibool_CPML,Qt_t);
    
    compute_accel_adepml(ispec_pml,I,J,K,&newtemp_adepml[0][0][0],
                d_wgllwgll_yz[K*NGLLX+J],d_wgllwgll_xz[K*NGLLX+I],
                d_wgllwgll_xy[J*NGLLX+I],rtrans,rtrans_inv,
              pml_kappa,pml_d,&sum_terms1,&sum_terms2,&sum_terms3);

    // assembles acceleration array
    if (threadIdx.x < NGLL3) {
      atomicAdd(&d_accel[iglob*3], -sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], -sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], -sum_terms3);
        
    } // threadIdx.x
    
    // PML contribution
    if(threadIdx.x < NGLL2) {
      add_pml_physical_contribution(ispec_pml,threadIdx.x,phy_spec,phy_ijk,
                              ibool_CPML,phy_norm,rtrans,rtrans_inv,
                            pml_d,phy_jaco2Dw,dsigma,Qt_t);
    }
  }

}


extern "C"
void compute_forces_viscoelastic_cuda_ade_(long* Mesh_pointer,
                                                int* iphase,
                                                realw* deltat,
                                                int* nspec_outer_elastic,
                                                int* nspec_inner_elastic,
                                                int* COMPUTE_AND_STORE_STRAIN,
                                                int* ATTENUATION,
                                                int* ANISOTROPY) {

  TRACE("\tcompute_forces_viscoelastic_cuda_ade");
  // EPIK_TRACER("compute_forces_viscoelastic_cuda");
  //printf("Running compute_forces\n");
  //double start_time = get_time();

  // compute_forces_viscoelastic_cuda_(Mesh_pointer,iphase,deltat,nspec_outer_elastic,nspec_inner_elastic,
  //                             COMPUTE_AND_STORE_STRAIN,ATTENUATION,ANISOTROPY);
  // return;

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int num_elements;

  if (*iphase == 1)
    num_elements = *nspec_outer_elastic;
  else
    num_elements = *nspec_inner_elastic;

  // checks if anything to do
  if (num_elements == 0) return;

  // gpu resources
  int blocksize = NGLL3_PADDED;
  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_elements,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);
  
  kernel_forces_adepml <<< grid,threads,0,mp->compute_stream>>> (
    num_elements,mp->d_ibool,mp->d_phase_ispec_inner_elastic,
    mp->num_phase_ispec_elastic,*iphase,mp->d_irregular_element_number,
    mp->d_displ,mp->d_accel,*ANISOTROPY,mp->d_xix,mp->d_xiy,mp->d_xiz,
    mp->d_etax,mp->d_etay,mp->d_etaz,mp->d_gammax,mp->d_gammay,mp->d_gammaz,
    mp->xix_regular,mp->jacobian_regular,mp->d_hprime_xx,mp->d_hprimewgll_xx,
    mp->d_wgllwgll_xy,mp->d_wgllwgll_xz,mp->d_wgllwgll_yz,mp->d_kappav,
    mp->d_muv,mp->d_c11store,mp->d_c12store,mp->d_c13store,mp->d_c14store,
    mp->d_c15store,mp->d_c16store,mp->d_c22store,mp->d_c23store,mp->d_c24store,
    mp->d_c25store,mp->d_c26store,mp->d_c33store,mp->d_c34store,mp->d_c35store,
    mp->d_c36store,mp->d_c44store,mp->d_c45store,mp->d_c46store,mp->d_c55store,
    mp->d_c56store,mp->d_c66store,mp->is_wavefield_discontinuity,mp->d_displ_wd,
    mp->d_ispec_to_elem_wd,mp->d_ibool_wd,mp->d_is_pml,mp->d_spec_to_CPML,
    mp->d_ibool_CPML,mp->d_r_trans,mp->d_r_trans_inv,mp->d_pml_kappa,
    mp->d_pml_d,mp->d_pml_physical_ijk,mp->d_pml_spec_physical,
    mp->d_pml_physical_normal,mp->d_pml_physical_jacobian2Dw,mp->d_Qu,
    mp->d_Qu_t,mp->d_Qt_t
  );

}

#endif