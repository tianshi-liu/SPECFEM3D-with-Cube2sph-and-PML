#include "mesh_constants_cuda.h"


__global__ static void 
kernel_prepare_boundary_matrix (realw* array_val,int ndim,int num_interfaces,
                         int max_nibool_interfaces, const int *nibool_interfaces,
                          const int *ibool_interfaces,
                          realw * buffer_send_matrix,int recv_stage)
{
    int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

    // shape is like (num_interface,max_nibool,ndim)
    int iinterface = id / (max_nibool_interfaces * ndim);
    id = id % (max_nibool_interfaces * ndim);
    int ipoin = id / ndim;
    int i = id % ndim;

    if(iinterface >= num_interfaces || i >= ndim || ipoin >=nibool_interfaces[iinterface]) return;

    // entry in interface array
    int ientry = ipoin + max_nibool_interfaces*iinterface;
    // global index in wavefield
    int iglob = ibool_interfaces[ientry] - 1;

    if(recv_stage == 0) {
        buffer_send_matrix[ndim*ientry + i] = array_val[ndim*iglob + i];
    }
    else {
        atomicAdd(&array_val[ndim*iglob + i],buffer_send_matrix[ndim*ientry + i]);
    }
}


#ifndef USE_CUDA_AWARE_MPI

extern "C"
void sync_accel_bdry_buffers_(long *Mesh_pointer,const int *iphase,realw_p buffer)
{
    Mesh *mp = (Mesh*)(*Mesh_pointer);
    // asynchronous transfer from device to host

    TRACE("\tsync_accel_bdry_buffers");
    if (mp->size_mpi_buffer == 0) return;

    int blocksize = BLOCKSIZE_TRANSFER;
    int recv_stage = (*iphase == 2);
    size_t size = mp->size_mpi_buffer* sizeof(realw);

    int size_padded = mp->max_nibool_interfaces_ext_mesh*NDIM * mp->num_interfaces_ext_mesh;
    size_padded = (size_padded + blocksize - 1) / blocksize;
    int nx,ny;
    get_blocks_xy(size_padded,&nx,&ny);
    dim3 grid(nx,ny,1);
    dim3 threads(blocksize,1,1);

    if(recv_stage) {
        cudaMemcpyAsync(mp->d_send_accel_buffer,buffer,size,
                        cudaMemcpyHostToDevice,mp->compute_stream);
        kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
            mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
            mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
            mp->d_send_accel_buffer,1
        );
    }
    else {

        kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
            mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
            mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
            mp->d_send_accel_buffer,0
        );

        // waits until previous compute stream finishes
        cudaMemcpyAsync(buffer,mp->d_send_accel_buffer,size,
                        cudaMemcpyDeviceToHost,mp->compute_stream);
        cudaStreamSynchronize(mp->compute_stream);
    }

    cudaDeviceSynchronize();

    
}

extern "C"
void sync_ade_bdry_buffers_(long *Mesh_pointer,int *iphase,realw_p buffer)
{
    Mesh *mp = (Mesh*)(*Mesh_pointer);
    // asynchronous transfer from device to host

    TRACE("\tsync_ade_bdry_buffers");
    if (mp->size_mpi_buffer_pml == 0) return;

    int blocksize = BLOCKSIZE_TRANSFER;
    int recv_stage = (*iphase == 2);
    size_t size = mp->size_mpi_buffer_pml*sizeof(realw);

    int size_padded = mp->max_nibool_interfaces_PML*NDIM*NDIM* mp->num_interfaces_PML;
    size_padded = (size_padded + blocksize - 1) / blocksize;
    int nx,ny;
    get_blocks_xy(size_padded,&nx,&ny);
    dim3 grid(nx,ny,1);
    dim3 threads(blocksize,1,1);

    if(recv_stage) {
        cudaMemcpyAsync(mp->d_buffer_send_matrix_PML,buffer,size,
                            cudaMemcpyHostToDevice,mp->compute_stream);
        kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
            mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
            mp->max_nibool_interfaces_PML,
            mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
            mp->d_buffer_send_matrix_PML,recv_stage
        );
    }
    else {
        kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
            mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
            mp->max_nibool_interfaces_PML,
            mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
            mp->d_buffer_send_matrix_PML,recv_stage
        );

        // waits until previous compute stream finishes
        cudaMemcpyAsync(buffer,mp->d_buffer_send_matrix_PML,size,
                        cudaMemcpyDeviceToHost,mp->compute_stream);
        cudaStreamSynchronize(mp->compute_stream);
    }

    cudaDeviceSynchronize();
} 

#else

void assemble_asyn_send(int ndim,realw *d_buf_sd, realw *d_buf_rv,
                        int num_intfs,int max_nibool,const int* my_neighbors,
                        const int *nibool,MPI_Request *req_sd,
                        MPI_Request *req_rv,int tag)
{
    for(int it = 0; it < num_intfs; it ++) {
        int iloc = it * max_nibool * ndim;
        int npts = ndim * nibool[it];

        MPI_Isend(d_buf_sd + iloc,npts,MPI_FLOAT,
                  my_neighbors[it],tag,MPI_COMM_WORLD,
                  &req_sd[it]);
        MPI_Irecv(d_buf_rv + iloc,npts,MPI_FLOAT,
                  my_neighbors[it],tag,MPI_COMM_WORLD,
                  &req_rv[it]);
    }
}

extern "C"
void sync_accel_bdry_buffers_(long *Mesh_pointer,const int *iphase,
                              const int* my_neighbors,
                              const int *nibool)
{
    Mesh *mp = (Mesh*)(*Mesh_pointer);
    // asynchronous transfer from device to host

    TRACE("\tsync_accel_bdry_buffers");
    if (mp->size_mpi_buffer == 0) return;

    int blocksize = BLOCKSIZE_TRANSFER;
    int recv_stage = (*iphase == 2);

    int size_padded = mp->max_nibool_interfaces_ext_mesh*NDIM * mp->num_interfaces_ext_mesh;
    size_padded = (size_padded + blocksize - 1) / blocksize;
    int nx,ny;
    get_blocks_xy(size_padded,&nx,&ny);
    dim3 grid(nx,ny,1);
    dim3 threads(blocksize,1,1);

    if(recv_stage) {

        // wait for recv communication 
        MPI_Waitall(mp->num_interfaces_ext_mesh,mp->req_recv_ext,MPI_STATUS_IGNORE);

        kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
            mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
            mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
            mp->d_recv_accel_buffer,1
        );

        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);

        MPI_Waitall(mp->num_interfaces_ext_mesh,mp->req_send_ext,MPI_STATUS_IGNORE);
    }
    else {
        kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
            mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
            mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
            mp->d_send_accel_buffer,0
        );

        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);

        // send/recv cuda buffers
        assemble_asyn_send(NDIM,mp->d_send_accel_buffer,
                            mp->d_recv_accel_buffer,mp->num_interfaces_ext_mesh,
                            mp->max_nibool_interfaces_ext_mesh,my_neighbors,
                            nibool,mp->req_send_ext,mp->req_recv_ext,1242);
    }

    cudaDeviceSynchronize();
}

extern "C"
void sync_ade_bdry_buffers_(long *Mesh_pointer,const int *iphase,
                              const int* my_neighbors,
                              const int *nibool)
{
    Mesh *mp = (Mesh*)(*Mesh_pointer);
    // asynchronous transfer from device to host

    TRACE("\tsync_ade_bdry_buffers");
    if (mp->size_mpi_buffer_pml == 0) return;

    int blocksize = BLOCKSIZE_TRANSFER;
    int recv_stage = (*iphase == 2);

    int size_padded = mp->max_nibool_interfaces_PML*NDIM*NDIM* mp->num_interfaces_PML;
    size_padded = (size_padded + blocksize - 1) / blocksize;
    int nx,ny;
    get_blocks_xy(size_padded,&nx,&ny);
    dim3 grid(nx,ny,1);
    dim3 threads(blocksize,1,1);

    if(recv_stage) {

        // wait for recv communication 
        MPI_Waitall(mp->num_interfaces_PML,mp->req_recv_PML,MPI_STATUS_IGNORE);

        kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
            mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
            mp->max_nibool_interfaces_PML,
            mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
            mp->d_buffer_recv_matrix_PML,recv_stage
        );

        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);

        MPI_Waitall(mp->num_interfaces_PML,mp->req_send_PML,MPI_STATUS_IGNORE);
    }
    else {
        kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
            mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
            mp->max_nibool_interfaces_PML,
            mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
            mp->d_buffer_send_matrix_PML,0
        );

        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);

        // send/recv cuda buffers
        assemble_asyn_send(NDIM*NDIM,mp->d_buffer_send_matrix_PML,
                            mp->d_buffer_recv_matrix_PML,mp->num_interfaces_PML,
                            mp->max_nibool_interfaces_PML,my_neighbors,
                            nibool,mp->req_send_PML,mp->req_recv_PML,1567);
    }

    cudaDeviceSynchronize();
}

#endif