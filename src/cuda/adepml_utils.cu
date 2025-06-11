#include "mesh_constants_cuda.h"

__global__ static void 
kernel_update_adepml_accel(const int *CPML_to_glob,int nglob_CPML,
                          const realw *Qt,const realw* rvolume,
                          realw *__restrict__ accel)
{
    #define MXM(i,j) Qt[((iglob_CPML * NDIM) + j-1) * NDIM + i-1]
    int iglob_CPML = threadIdx.x + blockDim.x * blockIdx.x;
    if(iglob_CPML < nglob_CPML) {
        int iglob = CPML_to_glob[iglob_CPML] - 1;
        realw r = rvolume[iglob_CPML];
        realw ax = (MXM(1,1)+MXM(2,1)+MXM(3,1))  / r;
        realw ay = (MXM(1,2)+MXM(2,2)+MXM(3,2))  / r;
        realw az = (MXM(1,3)+MXM(2,3)+MXM(3,3))  / r;
        // atomicAdd(&accel[iglob*NDIM+0],ax);
        // atomicAdd(&accel[iglob*NDIM+1],ay);
        // atomicAdd(&accel[iglob*NDIM+2],az);
        accel[iglob*NDIM+0] +=  ax;
        accel[iglob*NDIM+1] +=  ay;
        accel[iglob*NDIM+2] +=  az;
    }
    #undef MXM
}

extern "C"
void include_adepml_accel_aux_device_(long *Mesh_pointer)
{
    TRACE("\tinclude_adepml_accel_aux_device");
    Mesh* mp = (Mesh*)(*Mesh_pointer);
    int nglob_CPML = mp->nglob_CPML;
    int nblocks = (nglob_CPML + BLOCKSIZE_KERNEL1 - 1) / BLOCKSIZE_KERNEL1;

    kernel_update_adepml_accel <<< nblocks,BLOCKSIZE_KERNEL1,0,mp->compute_stream >>> (
        mp->d_CPML_to_glob,nglob_CPML,mp->d_Qt,mp->d_rvolume,mp->d_accel
    );

}

__global__ static void 
kernel_update_Qu_conv(const realw *coeff_exp1, const realw* coeff_exp2,
                        int nspec_CPML, const realw* Qu_t,realw* __restrict__ Qu)
{
    int idx = blockIdx.x + gridDim.x*blockIdx.y; // = 
    int igll3 = threadIdx.x;
    if(idx < nspec_CPML * NDIM && igll3 < NGLL3){
        int ispec_CPML = idx / NDIM, idim = idx % NDIM;
        // Qu(3,3,NGLL3,nspec_CPML) coeff(3,NGLL3,nspec_CPML)
        for(int i = 0; i < NDIM; i ++) {
            int idx_q = ((ispec_CPML * NGLL3 + igll3) * NDIM + i) * NDIM + idim;
            int idx_c = (ispec_CPML * NGLL3 + igll3) * NDIM + idim;
            Qu[idx_q] = Qu[idx_q] * coeff_exp1[idx_c] + 
                        Qu_t[idx_q] * coeff_exp2[idx_c];
        }
    }
}

extern "C"
void update_qu_conv_device_(long *Mesh_pointer)
{
    TRACE("\tupdate_Qu_conv_device");
    Mesh* mp = (Mesh*)(*Mesh_pointer);
    int nspec_CPML = mp->nspec_pml;
    int nbx,nby;
    get_blocks_xy(nspec_CPML*NDIM,&nbx,&nby);
    dim3 grid(nbx,nby,1);
    dim3 block(NGLL3_PADDED,1,1);

    kernel_update_Qu_conv <<< grid,block,0,mp->compute_stream >>> (
        mp->d_coeff_exp1,mp->d_coeff_exp2,nspec_CPML,
        mp->d_Qu_t,mp->d_Qu
    );

    cudaMemsetAsync(mp->d_Qu_t,0,
                    sizeof(realw)*nspec_CPML*NDIM*NDIM*NGLL3,
                mp->compute_stream);

}

/**
 * do i=1,nglob_CMPL
 * Qt_t(:,:,i) = Qt_t(:,:,i)*rvolume(i)
 */
__global__ static  void 
kernel_update_Qt_conv1(int nglob_CPML, const realw* rvolume,
                      realw* __restrict__ Qt_t)
{
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    const int NDIM2 = NDIM * NDIM;
    if(idx < nglob_CPML * NDIM2) {
        int iglob_CPML = idx / NDIM2;
        int i = idx % NDIM2;
        int idx_q = (iglob_CPML * NDIM2 + i);
        Qt_t[idx_q] *= rvolume[iglob_CPML];
    }
}

/**
 *   do i = 1, 3 Qt(:,i,:)=Qt(:,i,:)*coeff_glob_exp1(:,:) + coeff_glob_exp2(:,:)*Qt_t(:,i,:)
 */
__global__ static  void 
kernel_update_Qt_conv2(int nglob_CPML, const realw *coeff_glob_exp1,
                        const realw* coeff_glob_exp2,realw* __restrict__ Qt,
                        const realw* Qt_t)
{
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if(idx < nglob_CPML * NDIM) {
        int iglob_CPML = idx / NDIM, idim = idx % NDIM;
        int idx_c = iglob_CPML * NDIM + idim;
        for(int i = 0; i < NDIM; i ++) {
            int idx_q = (iglob_CPML*NDIM+i)*NDIM + idim;
            Qt[idx_q] = Qt[idx_q] * coeff_glob_exp1[idx_c] + 
                        Qt_t[idx_q] * coeff_glob_exp2[idx_c];
        }
    }
}

extern "C" void
update_qt_conv_device_(long *Mesh_pointer)
{
    TRACE("\tupdate_Qt_conv_device");
    Mesh* mp = (Mesh*)(*Mesh_pointer);
    int nglob_CPML = mp->nglob_CPML;

    // GPU resources
    int nb1 = (nglob_CPML*NDIM + BLOCKSIZE_KERNEL1 - 1) / BLOCKSIZE_KERNEL1;
    int nb2 = (nglob_CPML*NDIM*NDIM + BLOCKSIZE_KERNEL1 - 1) / BLOCKSIZE_KERNEL1;
    
    // kernel_update_Qt_conv <<< nb,BLOCKSIZE_KERNEL1,0,mp->compute_stream >>> (
    //     nglob_CPML,mp->d_coeff_glob_exp1,mp->d_coeff_glob_exp2,mp->d_rvolume,
    //     mp->d_Qt,mp->d_Qt_t
    // );
    kernel_update_Qt_conv1 <<< nb2,BLOCKSIZE_KERNEL1,0,mp->compute_stream >>> (
        nglob_CPML,mp->d_rvolume,mp->d_Qt_t
    );
    kernel_update_Qt_conv2 <<< nb1,BLOCKSIZE_KERNEL1,0,mp->compute_stream >>> (
        nglob_CPML,mp->d_coeff_glob_exp1,mp->d_coeff_glob_exp2,
        mp->d_Qt,mp->d_Qt_t
    );

    cudaMemsetAsync(mp->d_Qt_t,0,sizeof(realw)*nglob_CPML*NDIM*NDIM,mp->compute_stream);
}