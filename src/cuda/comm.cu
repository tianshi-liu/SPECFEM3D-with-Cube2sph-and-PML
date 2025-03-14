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


// #ifndef USE_CUDA_AWARE_MPI

// extern "C"
// void sync_accel_bdry_buffers_(long *Mesh_pointer,const int *iphase,realw_p buffer)
// {
//     Mesh *mp = (Mesh*)(*Mesh_pointer);
//     // asynchronous transfer from device to host

//     TRACE("\tsync_accel_bdry_buffers");
//     if (mp->size_mpi_buffer == 0) return;

//     int blocksize = BLOCKSIZE_TRANSFER;
//     int recv_stage = (*iphase == 2);
//     size_t size = mp->size_mpi_buffer* sizeof(realw);

//     int size_padded = mp->max_nibool_interfaces_ext_mesh*NDIM * mp->num_interfaces_ext_mesh;
//     size_padded = (size_padded + blocksize - 1) / blocksize;
//     int nx,ny;
//     get_blocks_xy(size_padded,&nx,&ny);
//     dim3 grid(nx,ny,1);
//     dim3 threads(blocksize,1,1);

//     if(recv_stage) {
//         cudaMemcpyAsync(mp->d_send_accel_buffer,buffer,size,
//                         cudaMemcpyHostToDevice,mp->compute_stream);
//         kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
//             mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
//             mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
//             mp->d_send_accel_buffer,1
//         );
//     }
//     else {

//         kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
//             mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
//             mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
//             mp->d_send_accel_buffer,0
//         );

//         // waits until previous compute stream finishes
//         cudaMemcpyAsync(buffer,mp->d_send_accel_buffer,size,
//                         cudaMemcpyDeviceToHost,mp->compute_stream);
//         cudaStreamSynchronize(mp->compute_stream);
//     }

//     cudaDeviceSynchronize();

    
// }

// extern "C"
// void sync_ade_bdry_buffers_(long *Mesh_pointer,int *iphase,realw_p buffer)
// {
//     Mesh *mp = (Mesh*)(*Mesh_pointer);
//     // asynchronous transfer from device to host

//     TRACE("\tsync_ade_bdry_buffers");
//     if (mp->size_mpi_buffer_pml == 0) return;

//     int blocksize = BLOCKSIZE_TRANSFER;
//     int recv_stage = (*iphase == 2);
//     size_t size = mp->size_mpi_buffer_pml*sizeof(realw);

//     int size_padded = mp->max_nibool_interfaces_PML*NDIM*NDIM* mp->num_interfaces_PML;
//     size_padded = (size_padded + blocksize - 1) / blocksize;
//     int nx,ny;
//     get_blocks_xy(size_padded,&nx,&ny);
//     dim3 grid(nx,ny,1);
//     dim3 threads(blocksize,1,1);

//     if(recv_stage) {
//         cudaMemcpyAsync(mp->d_buffer_send_matrix_PML,buffer,size,
//                             cudaMemcpyHostToDevice,mp->compute_stream);
//         kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
//             mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
//             mp->max_nibool_interfaces_PML,
//             mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
//             mp->d_buffer_send_matrix_PML,recv_stage
//         );
//     }
//     else {
//         kernel_prepare_boundary_matrix <<<grid,threads,0,mp->compute_stream>>>(
//             mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
//             mp->max_nibool_interfaces_PML,
//             mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
//             mp->d_buffer_send_matrix_PML,recv_stage
//         );

//         // waits until previous compute stream finishes
//         cudaMemcpyAsync(buffer,mp->d_buffer_send_matrix_PML,size,
//                         cudaMemcpyDeviceToHost,mp->compute_stream);
//         cudaStreamSynchronize(mp->compute_stream);
//     }

//     cudaDeviceSynchronize();
// } 

//#else

/**
 * @brief assemble mpi buffers across different procs
 * @param buf_sd/rv mpi buffers, shape(num_intfs,max_nibool,ndim), it could be host/device memory
 * @param my_neighbors shape(num_intfs) all neigbors need to exchange information
 * @param num_intfs no. of neighbors
 * @param nibool shape(num_intfs) no. of points needto exchange for each neighbor
 * @param req_sd/rv MPI asynchronize requests, shape(num_intfs)
 */
void assemble_asyn_send(int ndim,realw *buf_sd, realw *buf_rv,
                        int num_intfs,int max_nibool,const int* my_neighbors,
                        const int *nibool,MPI_Request *req_sd,
                        MPI_Request *req_rv,int tag)
{
    for(int it = 0; it < num_intfs; it ++) {
        int iloc = it * max_nibool * ndim;
        int npts = ndim * nibool[it];

        MPI_Isend(buf_sd + iloc,npts,MPI_FLOAT,
                  my_neighbors[it],tag,MPI_COMM_WORLD,
                  &req_sd[it]);
        MPI_Irecv(buf_rv + iloc,npts,MPI_FLOAT,
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

    // buffer pointer
    realw *buf_sd = mp->d_send_accel_buffer, *buf_rv = mp->d_recv_accel_buffer;

    if(recv_stage) {

        // wait for recv communication 
        MPI_Waitall(mp->num_interfaces_ext_mesh,mp->req_recv_ext,MPI_STATUS_IGNORE);

#ifndef USE_CUDA_AWARE_MPI
        // copy from pinned host to device
        size_t size_cp = mp->size_mpi_buffer * sizeof(realw);
        cudaMemcpyAsync(mp->d_recv_accel_buffer,mp->h_recv_accel_buffer,
                        size_cp,cudaMemcpyHostToDevice,
                        mp->compute_stream);
#endif
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

#ifndef USE_CUDA_AWARE_MPI
        size_t size_cp = mp->size_mpi_buffer * sizeof(realw);
        cudaMemcpyAsync(mp->h_send_accel_buffer,mp->d_send_accel_buffer,
                        size_cp,cudaMemcpyDeviceToHost,
                        mp->compute_stream);
        buf_sd = mp->h_send_accel_buffer;
        buf_rv = mp->h_recv_accel_buffer;
#endif

        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);

        // send/recv cuda buffers
        assemble_asyn_send(NDIM,buf_sd,buf_rv,mp->num_interfaces_ext_mesh,
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

    // buffer pointer
    realw *buf_sd = mp->d_buffer_send_matrix_PML, 
          *buf_rv = mp->d_buffer_recv_matrix_PML;

    if(recv_stage) {

        // wait for recv communication 
        MPI_Waitall(mp->num_interfaces_PML,mp->req_recv_PML,MPI_STATUS_IGNORE);

#ifndef USE_CUDA_AWARE_MPI
        // copy from pinned host to device
        size_t size_cp = mp->size_mpi_buffer_pml * sizeof(realw);
        cudaMemcpyAsync(mp->d_buffer_recv_matrix_PML,mp->h_buffer_recv_matrix_PML,
                        size_cp,cudaMemcpyHostToDevice,mp->compute_stream);
#endif

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

#ifndef USE_CUDA_AWARE_MPI
        size_t size_cp = mp->size_mpi_buffer_pml * sizeof(realw);
        cudaMemcpyAsync(mp->h_buffer_send_matrix_PML,mp->d_buffer_send_matrix_PML,
                        size_cp,cudaMemcpyDeviceToHost,
                        mp->compute_stream);
        buf_sd = mp->h_buffer_send_matrix_PML;
        buf_rv = mp->h_buffer_recv_matrix_PML;
#endif

        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);

        // send/recv cuda buffers
        assemble_asyn_send(NDIM*NDIM,buf_sd,buf_rv,mp->num_interfaces_PML,
                            mp->max_nibool_interfaces_PML,my_neighbors,
                            nibool,mp->req_send_PML,mp->req_recv_PML,1567);
    }

    cudaDeviceSynchronize();
}



extern "C"
void sync_ade_bdry(long *Mesh_pointer,const int *iphase,
                    const int* my_neighbors,
                    const int* my_neighbors_PML,
                    const int *nibool)
{
    Mesh *mp = (Mesh*)(*Mesh_pointer);
    // asynchronous transfer from device to host

    TRACE("\tsync_ade_bdry_buffers");
    if (mp->size_mpi_buffer == 0) return;

    int blocksize = BLOCKSIZE_TRANSFER;
    int recv_stage = (*iphase == 2);

    // PML 
    int size_padded = mp->max_nibool_interfaces_PML*NDIM*NDIM* mp->num_interfaces_PML;
    size_padded = (size_padded + blocksize - 1) / blocksize;
    int nx,ny;
    get_blocks_xy(size_padded,&nx,&ny);
    dim3 grid_pml(nx,ny,1);
    dim3 threads(blocksize,1,1);

    // accel 
    size_padded = mp->max_nibool_interfaces_ext_mesh*NDIM* mp->num_interfaces_ext_mesh;
    size_padded = (size_padded + blocksize - 1) / blocksize;
    get_blocks_xy(size_padded,&nx,&ny);
    dim3 grid_ext(nx,ny,1);

    // buffer pointer
    realw *buf_sd_PML = mp->d_buffer_send_matrix_PML, 
          *buf_rv_PML = mp->d_buffer_recv_matrix_PML,
          *buf_sd = mp->d_send_accel_buffer,
          *buf_rv = mp->d_recv_accel_buffer;

    if(recv_stage) {

        // wait for accel recv communication 
        MPI_Waitall(mp->num_interfaces_ext_mesh,mp->req_recv_ext,MPI_STATUS_IGNORE);
        MPI_Waitall(mp->num_interfaces_PML,mp->req_recv_PML,MPI_STATUS_IGNORE);

        MPI_Waitall(mp->num_interfaces_ext_mesh,mp->req_send_ext,MPI_STATUS_IGNORE);
        MPI_Waitall(mp->num_interfaces_PML,mp->req_send_PML,MPI_STATUS_IGNORE);

#ifndef USE_CUDA_AWARE_MPI
        size_t size_cp = mp->size_mpi_buffer * sizeof(realw);
        cudaMemcpyAsync(mp->d_recv_accel_buffer,mp->h_recv_accel_buffer,
                        size_cp,cudaMemcpyHostToDevice,
                        mp->compute_stream);
        
        // copy from pinned host to device
        size_cp = mp->size_mpi_buffer_pml * sizeof(realw);
        cudaMemcpyAsync(mp->d_buffer_recv_matrix_PML,mp->h_buffer_recv_matrix_PML,
                        size_cp,cudaMemcpyHostToDevice,mp->compute_stream);
#endif

        kernel_prepare_boundary_matrix <<<grid_ext,threads,0,mp->compute_stream>>>(
            mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
            mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
            mp->d_recv_accel_buffer,1
        );

        kernel_prepare_boundary_matrix <<<grid_pml,threads,0,mp->compute_stream>>>(
            mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
            mp->max_nibool_interfaces_PML,
            mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
            mp->d_buffer_recv_matrix_PML,recv_stage
        );

        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);
    }
    else {

        kernel_prepare_boundary_matrix <<<grid_ext,threads,0,mp->compute_stream>>>(
            mp->d_accel,NDIM,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
            mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh,
            mp->d_send_accel_buffer,0
        );

        kernel_prepare_boundary_matrix <<<grid_pml,threads,0,mp->compute_stream>>>(
            mp->d_Qt_t,NDIM*NDIM,mp->num_interfaces_PML,
            mp->max_nibool_interfaces_PML,
            mp->d_nibool_interfaces_PML,mp->d_ibool_interfaces_PML,
            mp->d_buffer_send_matrix_PML,0
        );

#ifndef USE_CUDA_AWARE_MPI
        size_t size_cp = mp->size_mpi_buffer_pml * sizeof(realw);
        cudaMemcpyAsync(mp->h_buffer_send_matrix_PML,mp->d_buffer_send_matrix_PML,
                        size_cp,cudaMemcpyDeviceToHost,
                        mp->compute_stream);
        buf_sd_PML = mp->h_buffer_send_matrix_PML;
        buf_rv_PML = mp->h_buffer_recv_matrix_PML;

        size_cp = mp->size_mpi_buffer * sizeof(realw);
        cudaMemcpyAsync(mp->h_send_accel_buffer,mp->d_send_accel_buffer,
            size_cp,cudaMemcpyDeviceToHost,
            mp->compute_stream);
        buf_sd = mp->h_send_accel_buffer;
        buf_rv = mp->h_recv_accel_buffer;
#endif

        // waits until previous compute stream finishes
        cudaStreamSynchronize(mp->compute_stream);

        // send/recv accel buffers
        assemble_asyn_send(NDIM,buf_sd,buf_rv,mp->num_interfaces_ext_mesh,
            mp->max_nibool_interfaces_ext_mesh,my_neighbors,
            nibool,mp->req_send_ext,mp->req_recv_ext,1242);

        // send/recv qt_t buffers
        assemble_asyn_send(NDIM*NDIM,buf_sd_PML,buf_rv_PML,mp->num_interfaces_PML,
                            mp->max_nibool_interfaces_PML,my_neighbors,
                            nibool,mp->req_send_PML,mp->req_recv_PML,1567);
    }

    cudaDeviceSynchronize();
}

// #endif