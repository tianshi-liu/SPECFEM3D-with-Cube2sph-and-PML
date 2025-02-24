#include <cuda_runtime.h>

const int  NGLL =  5;
const int  NGLLX = 5;
const int  NGLL2 = 25;
const int  NGLL3 = 125;
const int NGLL3PAD = 128;

typedef float realw;
typedef const realw* __restrict__ realw_const_p;
typedef realw* __restrict__ realw_p;

__device__ __forceinline__ void 
mxm_opt(const realw *u, const realw *hprime,int I,int J, int K,
        int comp,realw_p temp1)
{
    *temp1 = 0.;
    switch (comp)
    {
    case 1:
        for(int l = 0; l < NGLL; l ++) {
            *temp1 += u[K*NGLL2+J*NGLL+l] * hprime[I*NGLL+l];
        }
        break;
    case 2:
        for(int l = 0; l < NGLL; l ++) {
            *temp1 += u[K*NGLL2+l*NGLL+I] * hprime[J*NGLL+l];
        }
        break;
    default:
        for(int l = 0; l < NGLL; l ++) {
            *temp1 += u[l*NGLL2+J*NGLL+I] * hprime[K*NGLL+l];
        }
        break;
    }
}


__global__ void 
kernel_smooth_sph_pde(int num_elmts,int iphase,int num_phase,
                    const int* phase_ispec_inner_elastic,
                    const realw* xix,const realw* xiy,const realw* xiz,
                    const realw* etax,const realw* etay,const realw* etaz,
                    const realw* gamx,const realw* gamy,const realw* gamz,
                    const realw *jaco, const realw* rotate, 
                    const realw* wgllwgll_xy,const realw* wgllwgll_xz,
                    const realw* wgllwgll_yz,const realw* hprimeT,
                    const realw* hprime_wgll,const int *ibool,const realw cv,const realw ch,
                    const realw* dat,realw_p ddat_glob)
{
    int bx = blockIdx.x + gridDim.x*blockIdx.y;
    if(bx >= num_elmts) return;
    int ispec = phase_ispec_inner_elastic[bx + (iphase - 1) * num_phase] - 1;
    int tx = threadIdx.x;
    if(tx >= NGLL3) tx = NGLL3 - 1; 
    int offset = ispec * NGLL3PAD + tx;
    int iglob = ibool[offset] - 1;

    // ijk
    int K,J,I;
    K = (tx/NGLL2);
    J = ((tx-K*NGLL2)/NGLLX);
    I = (tx-K*NGLL2-J*NGLLX);

    // shared memory
    __shared__ realw sh_u[NGLL3];
    __shared__ realw sh_hprimeT[NGLL2], sh_hprimewgll[NGLL2];

    // copy shared memory
    if(threadIdx.x < NGLL3) {
        sh_u[tx] = dat[ispec * NGLL3 + tx];
        if(threadIdx.x < NGLL2) {
            sh_hprimeT[tx] = hprimeT[tx];
            sh_hprimewgll[tx] = hprime_wgll[tx];
        }
    }
    __syncthreads();

    // compute gradient
    realw dudx{},dudy{},dudz{};
    realw temp1{},temp2{},temp3{};
    realw rl,rxl,ryl,rzl,xixl,xiyl,xizl,etaxl,etayl,
          etazl,gamxl,gamyl,gamzl,jacobianl;
    rxl = rotate[(1-1) * offset];
    ryl = rotate[(2-1) * offset];
    rzl = rotate[(3-1) * offset];
    xixl = xix[offset];
    xiyl = xiy[offset];
    xizl = xiz[offset];
    etaxl = etax[offset];
    etayl = etay[offset];
    etazl = etaz[offset];
    gamxl = gamx[offset];
    gamyl = gamy[offset];
    gamzl = gamz[offset];
    jacobianl = jaco[offset];
    mxm_opt(sh_u,sh_hprimeT,I,J,K,1,&temp1);
    mxm_opt(sh_u,sh_hprimeT,I,J,K,2,&temp2);
    mxm_opt(sh_u,sh_hprimeT,I,J,K,3,&temp3);
    dudx = temp1 * xixl + temp2 * etaxl + temp3 * gamxl;
    dudy = temp1 * xiyl + temp2 * etayl + temp3 * gamyl;
    dudz = temp1 * xizl + temp2 * etazl + temp3 * gamzl;

    // compute new terms
    if(threadIdx.x < NGLL3) {
        sh_u[tx] = ((cv-ch) * (rxl*xixl+ryl*xiyl+rzl*xizl) * 
                    (rxl*dudx+ryl*dudy+rzl*dudz) +
                    ch * (xixl*dudx+xiyl*dudy + 
                    xizl*dudz)) * jacobianl;
    }
    __syncthreads();
    mxm_opt(sh_u,sh_hprimewgll,I,J,K,1,&temp1);

    __syncthreads();
    if(threadIdx.x < NGLL3) {
        sh_u[tx] = ((cv-ch) * (rxl*etaxl+ryl*etayl+rzl*etazl) * 
                    (rxl*dudx+ryl*dudy+rzl*dudz) +
                    ch * (etaxl*dudx+etayl*dudy +
                    etazl*dudz)) * jacobianl;
    }
    __syncthreads();
    mxm_opt(sh_u,sh_hprimewgll,I,J,K,2,&temp2);

    __syncthreads();
    if(threadIdx.x < NGLL3) {
        sh_u[tx] = ((cv-ch) * (rxl*gamxl+ryl*gamyl+rzl*gamzl) * 
                    (rxl*dudx+ryl*dudy+rzl*dudz) +
                    ch * (gamxl*dudx+gamyl*dudy + 
                    gamzl*dudz)) * jacobianl;
    }
    __syncthreads();
    mxm_opt(sh_u,sh_hprimewgll,I,J,K,3,&temp3);

    // update ddat_glob
    realw fac1 = wgllwgll_yz[J*NGLL + K];
    realw fac2 = wgllwgll_xz[I*NGLL + K];
    realw fac3 = wgllwgll_xy[I*NGLL + J];
    realw sum_term = fac1*temp1+ fac2*temp2 + fac3*temp3;
    atomicAdd(&ddat_glob[iglob],-sum_term);
}


