#include "mesh_constants_cuda.h"


__global__ static void 
kernel_prepare_boundary_matrix (realw* array_val,int ndim,int num_interfaces,
                         int max_nibool_interfaces, const int *nibool_interfaces,
                          const int *ibool_interfaces,
                          realw * buffer_send_matrix,int buffer2arr)
{
  int ipoin = threadIdx.x + blockIdx.x * blockDim.x;
  
  for( int iinterface=0; iinterface < num_interfaces; iinterface++) {
    if (ipoin < nibool_interfaces[iinterface]) {
      // entry in interface array
      int ientry = ipoin + max_nibool_interfaces*iinterface;
      // global index in wavefield
      int iglob = ibool_interfaces[ientry] - 1;
      if(buffer2arr == 0) {
        for(int i = 0; i < ndim; i ++)
          buffer_send_matrix[ndim*ientry + i] = array_val[ndim*iglob + i];
      }
      else {
        for(int i = 0; i < ndim; i ++)
          atomicAdd(&array_val[ndim*iglob + i],buffer_send_matrix[ndim*ientry + i]);
      }
      
    }
  }
}


extern "C"
void sync_accel_bdry_buffers_(long *Mesh_pointer,const int *iphase,realw_p buffer)
{
    Mesh *mp = (Mesh*)(*Mesh_pointer);
    // asynchronous transfer from device to host

    TRACE("\tsync_accel_bdry_buffers");
    if (mp->size_mpi_buffer == 0) return;

    int blocksize = BLOCKSIZE_TRANSFER;
    int buf2arr = (*iphase == 2);
    size_t size = mp->size_mpi_buffer* sizeof(realw);

    if(buf2arr) {
        cudaMemcpy(mp->d_send_accel_buffer,buffer,size,cudaMemcpyHostToDevice);
    }

    
    int size_padded = (mp->max_nibool_interfaces_ext_mesh + blocksize - 1) / blocksize;

    dim3 grid(size_padded,1,1);
    dim3 threads(blocksize,1,1);

    kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
        mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
        mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
        mp->d_send_accel_buffer,buf2arr
    );

    // waits until previous compute stream finishes
    cudaStreamSynchronize(mp->compute_stream);

    if(!buf2arr) {
        cudaMemcpy(buffer,mp->d_send_accel_buffer,size,cudaMemcpyDeviceToHost);
    }
}

extern "C"
void sync_ade_bdry_buffers_(long *Mesh_pointer,int *iphase,realw_p buffer)
{
    Mesh *mp = (Mesh*)(*Mesh_pointer);
    // asynchronous transfer from device to host

    TRACE("\tsync_accel_bdry_buffers");
    if (mp->size_mpi_buffer_pml == 0) return;

    int blocksize = BLOCKSIZE_TRANSFER;
    int buf2arr = (*iphase == 2);
    size_t size = mp->size_mpi_buffer_pml*sizeof(realw);

    if(buf2arr) {
        cudaMemcpy(mp->d_buffer_send_matrix_PML,buffer,size,cudaMemcpyHostToDevice);
    }

    
    int size_padded = (mp->max_nibool_interfaces_PML + blocksize - 1) / blocksize;

    dim3 grid(size_padded,1,1);
    dim3 threads(blocksize,1,1);

    kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
        mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
        mp->max_nibool_interfaces_PML,
        mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
        mp->d_buffer_send_matrix_PML,buf2arr
    );

    if(!buf2arr) {
        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);
        cudaMemcpy(buffer,mp->d_buffer_send_matrix_PML,size,cudaMemcpyDeviceToHost);
    }
}