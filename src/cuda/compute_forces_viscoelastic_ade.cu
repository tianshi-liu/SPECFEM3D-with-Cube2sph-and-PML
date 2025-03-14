

#include "utils.cuh"


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

/**
 * sig[1-6]: xx,xy,xz,yy,yz,zz
 */
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
                    realw_const_p pml_kappa,
                    realw_p sig1,realw_p sig2,realw_p sig3,
                    realw_p sig4,realw_p sig5,realw_p sig6,
                    const realw* Qu)
{
  const int idx = ispec * NGLL3 + k * NGLL2 + j * NGLLX + i;

  #define Rmat(i,j)  rtrans[(idx*NDIM + j-1) * NDIM + i-1]
  #define Rmatinv(i,j) rtrans_inv[(idx*NDIM + j-1) * NDIM + i-1]
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
    *sig1 = lambdalplus2mul*e11+lambdal*(e22+e33); // xx
    *sig2 = mul*e12; // xy
    *sig3 = mul*e13;; // xz
    *sig4 = lambdalplus2mul*e22+lambdal*(e11+e33); // yy
    *sig5 = mul*e23;; // yz
    *sig6 = lambdalplus2mul*e33+lambdal*(e11+e22); //zz
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

    *sig1 = c11 * e11 + c16 * e12 + c12 * e22 +
        c15 * e13 + c14 * e23 + c13 * e33;
    *sig4 = c12 * e11 + c26 * e12 + c22 * e22 +
        c25 * e13 + c24 * e23 + c23 * e33;
    *sig6 = c13 * e11 + c36 * e12 + c23 * e22 +
        c35 * e13 + c34 * e23 + c33 * e33;
    *sig2 = c16 * e11 + c66 * e12 + c26 * e22 +
        c56 * e13 + c46 * e23 + c36 * e33;
    *sig3 = c15 * e11 + c56 * e12 + c25 * e22 +
        c55 * e13 + c45 * e23 + c35 * e33;
    *sig5 = c14 * e11 + c46 * e12 + c24 * e22 +
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
  realw d1l,d2l,d3l;
  const int idx = ispec * NGLL3 + k * NGLL2 + j * NGLLX + i;
  
  #define Rmatinv(i,j) rtrans_inv[(idx * NDIM + j-1) * NDIM + i-1]
  #define Qmat(i,j) Qu_t[(idx * NDIM + j-1) * NDIM + i-1]
  const realw ri1xl = Rmatinv(1,1);
  const realw ri2xl = Rmatinv(2,1);
  const realw ri3xl = Rmatinv(3,1);
  const realw ri1yl = Rmatinv(1,2);
  const realw ri2yl = Rmatinv(2,2);
  const realw ri3yl = Rmatinv(3,2);
  const realw ri1zl = Rmatinv(1,3);
  const realw ri2zl = Rmatinv(2,3);
  const realw ri3zl = Rmatinv(3,3);

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

  realw rt[9],rtinv[9];
  for(int i = 0; i < 9; i ++) {
    rt[i] = r_trans[idx*NDIM*NDIM + i];
    rtinv[i] = r_trans_inv[idx*NDIM*NDIM + i];
  }

  // #define Rmat(i,j) r_trans[(idx*NDIM + j-1) * NDIM + i-1]
  // #define Rmatinv(i,j) r_trans_inv[(idx*NDIM + j-1) * NDIM + i-1]
  #define Rmat(i,j) rt[(j-1) * NDIM + i-1]
  #define Rmatinv(i,j) rtinv[(j-1) * NDIM + i-1]

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

  d1l = pml_d[1-1 + idx*NDIM];
  d2l = pml_d[2-1 + idx*NDIM];
  d3l = pml_d[3-1 + idx*NDIM];

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
compute_newcomp(realw_p sh_tempx,realw_p sh_tempy,realw_p sh_tempz,
               realw_p sh_hprimewgll_xx,realw sigma_xx,realw sigma_xy,
               realw sigma_xz,realw sigma_yy,realw sigma_yz,
               realw sigma_zz,realw xixj,realw xiyj,realw xizj,
               realw etaxj,realw etayj, realw etazj, 
               realw gamxj,realw gamyj,realw gamzj,int tx,
               int I,int J,int K, realw_p newtemp_adepml)
{
  const int N = 18;
  realw ux,uy,uz;


  //   if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xx*xixj;
  //     sh_tempy[tx] = sigma_xx*xiyj;
  //     sh_tempz[tx] = sigma_xx*xizj;
  // }
  // __syncthreads();
  // sum_hprimewgll_xi(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[0 + 0 * 6] = ux;
  // newtemp_adepml[0 + 1 * 6] = uy;
  // newtemp_adepml[0 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xy*xixj;
  //     sh_tempy[tx] = sigma_xy*xiyj;
  //     sh_tempz[tx] = sigma_xy*xizj;
  // }
  // __syncthreads();
  // sum_hprimewgll_xi(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[1 + 0 * 6] = ux;
  // newtemp_adepml[1 + 1 * 6] = uy;
  // newtemp_adepml[1 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xz*xixj;
  //     sh_tempy[tx] = sigma_xz*xiyj;
  //     sh_tempz[tx] = sigma_xz*xizj;
  // }
  // __syncthreads();
  // sum_hprimewgll_xi(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[2 + 0 * 6] = ux;
  // newtemp_adepml[2 + 1 * 6] = uy;
  // newtemp_adepml[2 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_yy*xixj;
  //     sh_tempy[tx] = sigma_yy*xiyj;
  //     sh_tempz[tx] = sigma_yy*xizj;
  // }
  // __syncthreads();
  // sum_hprimewgll_xi(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[3 + 0 * 6] = ux;
  // newtemp_adepml[3 + 1 * 6] = uy;
  // newtemp_adepml[3 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_yz*xixj;
  //     sh_tempy[tx] = sigma_yz*xiyj;
  //     sh_tempz[tx] = sigma_yz*xizj;
  // }
  // __syncthreads();
  // sum_hprimewgll_xi(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[4 + 0 * 6] = ux;
  // newtemp_adepml[4 + 1 * 6] = uy;
  // newtemp_adepml[4 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_zz*xixj;
  //     sh_tempy[tx] = sigma_zz*xiyj;
  //     sh_tempz[tx] = sigma_zz*xizj;
  // }
  // __syncthreads();
  // sum_hprimewgll_xi(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[5 + 0 * 6] = ux;
  // newtemp_adepml[5 + 1 * 6] = uy;
  // newtemp_adepml[5 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xx*etaxj;
  //     sh_tempy[tx] = sigma_xx*etayj;
  //     sh_tempz[tx] = sigma_xx*etazj;
  // }
  // __syncthreads();
  // sum_hprimewgll_eta(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[18 + 0 * 6] = ux;
  // newtemp_adepml[18 + 1 * 6] = uy;
  // newtemp_adepml[18 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xy*etaxj;
  //     sh_tempy[tx] = sigma_xy*etayj;
  //     sh_tempz[tx] = sigma_xy*etazj;
  // }
  // __syncthreads();
  // sum_hprimewgll_eta(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[19 + 0 * 6] = ux;
  // newtemp_adepml[19 + 1 * 6] = uy;
  // newtemp_adepml[19 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xz*etaxj;
  //     sh_tempy[tx] = sigma_xz*etayj;
  //     sh_tempz[tx] = sigma_xz*etazj;
  // }
  // __syncthreads();
  // sum_hprimewgll_eta(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[20 + 0 * 6] = ux;
  // newtemp_adepml[20 + 1 * 6] = uy;
  // newtemp_adepml[20 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_yy*etaxj;
  //     sh_tempy[tx] = sigma_yy*etayj;
  //     sh_tempz[tx] = sigma_yy*etazj;
  // }
  // __syncthreads();
  // sum_hprimewgll_eta(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[21 + 0 * 6] = ux;
  // newtemp_adepml[21 + 1 * 6] = uy;
  // newtemp_adepml[21 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_yz*etaxj;
  //     sh_tempy[tx] = sigma_yz*etayj;
  //     sh_tempz[tx] = sigma_yz*etazj;
  // }
  // __syncthreads();
  // sum_hprimewgll_eta(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[22 + 0 * 6] = ux;
  // newtemp_adepml[22 + 1 * 6] = uy;
  // newtemp_adepml[22 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_zz*etaxj;
  //     sh_tempy[tx] = sigma_zz*etayj;
  //     sh_tempz[tx] = sigma_zz*etazj;
  // }
  // __syncthreads();
  // sum_hprimewgll_eta(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[23 + 0 * 6] = ux;
  // newtemp_adepml[23 + 1 * 6] = uy;
  // newtemp_adepml[23 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xx*gamxj;
  //     sh_tempy[tx] = sigma_xx*gamyj;
  //     sh_tempz[tx] = sigma_xx*gamzj;
  // }
  // __syncthreads();
  // sum_hprimewgll_gamma(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[36 + 0 * 6] = ux;
  // newtemp_adepml[36 + 1 * 6] = uy;
  // newtemp_adepml[36 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xy*gamxj;
  //     sh_tempy[tx] = sigma_xy*gamyj;
  //     sh_tempz[tx] = sigma_xy*gamzj;
  // }
  // __syncthreads();
  // sum_hprimewgll_gamma(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[37 + 0 * 6] = ux;
  // newtemp_adepml[37 + 1 * 6] = uy;
  // newtemp_adepml[37 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_xz*gamxj;
  //     sh_tempy[tx] = sigma_xz*gamyj;
  //     sh_tempz[tx] = sigma_xz*gamzj;
  // }
  // __syncthreads();
  // sum_hprimewgll_gamma(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[38 + 0 * 6] = ux;
  // newtemp_adepml[38 + 1 * 6] = uy;
  // newtemp_adepml[38 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_yy*gamxj;
  //     sh_tempy[tx] = sigma_yy*gamyj;
  //     sh_tempz[tx] = sigma_yy*gamzj;
  // }
  // __syncthreads();
  // sum_hprimewgll_gamma(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[39 + 0 * 6] = ux;
  // newtemp_adepml[39 + 1 * 6] = uy;
  // newtemp_adepml[39 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_yz*gamxj;
  //     sh_tempy[tx] = sigma_yz*gamyj;
  //     sh_tempz[tx] = sigma_yz*gamzj;
  // }
  // __syncthreads();
  // sum_hprimewgll_gamma(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[40 + 0 * 6] = ux;
  // newtemp_adepml[40 + 1 * 6] = uy;
  // newtemp_adepml[40 + 2 * 6] = uz;

  // if(tx < NGLL3) {
  //     sh_tempx[tx] = sigma_zz*gamxj;
  //     sh_tempy[tx] = sigma_zz*gamyj;
  //     sh_tempz[tx] = sigma_zz*gamzj;
  // }
  // __syncthreads();
  // sum_hprimewgll_gamma(I,J,K,&ux,&uy,&uz,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);
  // __syncthreads();
  // newtemp_adepml[41 + 0 * 6] = ux;
  // newtemp_adepml[41 + 1 * 6] = uy;
  // newtemp_adepml[41 + 2 * 6] = uz;

  if(tx < NGLL3) {
    sh_tempx[tx] = sigma_xx*xixj;
    sh_tempy[tx] = sigma_xx*etaxj;
    sh_tempz[tx] = sigma_xx*gamxj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+0] = ux;
  newtemp_adepml[1*N+0] = uy;
  newtemp_adepml[2*N+0] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_xy*xixj;
      sh_tempy[tx] = sigma_xy*etaxj;
      sh_tempz[tx] = sigma_xy*gamxj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+1] = ux;
  newtemp_adepml[1*N+1] = uy;
  newtemp_adepml[2*N+1] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_xz*xixj;
      sh_tempy[tx] = sigma_xz*etaxj;
      sh_tempz[tx] = sigma_xz*gamxj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+2] = ux;
  newtemp_adepml[1*N+2] = uy;
  newtemp_adepml[2*N+2] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_yy*xixj;
      sh_tempy[tx] = sigma_yy*etaxj;
      sh_tempz[tx] = sigma_yy*gamxj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+3] = ux;
  newtemp_adepml[1*N+3] = uy;
  newtemp_adepml[2*N+3] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_yz*xixj;
      sh_tempy[tx] = sigma_yz*etaxj;
      sh_tempz[tx] = sigma_yz*gamxj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+4] = ux;
  newtemp_adepml[1*N+4] = uy;
  newtemp_adepml[2*N+4] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_zz*xixj;
      sh_tempy[tx] = sigma_zz*etaxj;
      sh_tempz[tx] = sigma_zz*gamxj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+5] = ux;
  newtemp_adepml[1*N+5] = uy;
  newtemp_adepml[2*N+5] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_xx*xiyj;
      sh_tempy[tx] = sigma_xx*etayj;
      sh_tempz[tx] = sigma_xx*gamyj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+6] = ux;
  newtemp_adepml[1*N+6] = uy;
  newtemp_adepml[2*N+6] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_xy*xiyj;
      sh_tempy[tx] = sigma_xy*etayj;
      sh_tempz[tx] = sigma_xy*gamyj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+7] = ux;
  newtemp_adepml[1*N+7] = uy;
  newtemp_adepml[2*N+7] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_xz*xiyj;
      sh_tempy[tx] = sigma_xz*etayj;
      sh_tempz[tx] = sigma_xz*gamyj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+8] = ux;
  newtemp_adepml[1*N+8] = uy;
  newtemp_adepml[2*N+8] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_yy*xiyj;
      sh_tempy[tx] = sigma_yy*etayj;
      sh_tempz[tx] = sigma_yy*gamyj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+9] = ux;
  newtemp_adepml[1*N+9] = uy;
  newtemp_adepml[2*N+9] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_yz*xiyj;
      sh_tempy[tx] = sigma_yz*etayj;
      sh_tempz[tx] = sigma_yz*gamyj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+10] = ux;
  newtemp_adepml[1*N+10] = uy;
  newtemp_adepml[2*N+10] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_zz*xiyj;
      sh_tempy[tx] = sigma_zz*etayj;
      sh_tempz[tx] = sigma_zz*gamyj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+11] = ux;
  newtemp_adepml[1*N+11] = uy;
  newtemp_adepml[2*N+11] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_xx*xizj;
      sh_tempy[tx] = sigma_xx*etazj;
      sh_tempz[tx] = sigma_xx*gamzj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+12] = ux;
  newtemp_adepml[1*N+12] = uy;
  newtemp_adepml[2*N+12] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_xy*xizj;
      sh_tempy[tx] = sigma_xy*etazj;
      sh_tempz[tx] = sigma_xy*gamzj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+13] = ux;
  newtemp_adepml[1*N+13] = uy;
  newtemp_adepml[2*N+13] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_xz*xizj;
      sh_tempy[tx] = sigma_xz*etazj;
      sh_tempz[tx] = sigma_xz*gamzj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+14] = ux;
  newtemp_adepml[1*N+14] = uy;
  newtemp_adepml[2*N+14] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_yy*xizj;
      sh_tempy[tx] = sigma_yy*etazj;
      sh_tempz[tx] = sigma_yy*gamzj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+15] = ux;
  newtemp_adepml[1*N+15] = uy;
  newtemp_adepml[2*N+15] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_yz*xizj;
      sh_tempy[tx] = sigma_yz*etazj;
      sh_tempz[tx] = sigma_yz*gamzj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+16] = ux;
  newtemp_adepml[1*N+16] = uy;
  newtemp_adepml[2*N+16] = uz;

  if(tx < NGLL3) {
      sh_tempx[tx] = sigma_zz*xizj;
      sh_tempy[tx] = sigma_zz*etazj;
      sh_tempz[tx] = sigma_zz*gamzj;
  }
  __syncthreads();
    mxm_3op(I,J,K,sh_hprimewgll_xx,sh_tempx,sh_tempy,sh_tempz,&ux,&uy,&uz);
  __syncthreads();
  newtemp_adepml[0*N+17] = ux;
  newtemp_adepml[1*N+17] = uy;
  newtemp_adepml[2*N+17] = uz;
}

__device__ __forceinline__ void 
compute_accel_adepml(int ispec,int i,int j,int k,realw_const_p temp,
                    realw fac1,realw fac2,realw fac3,realw_const_p r_trans, 
                    realw_const_p r_trans_inv, realw_const_p pml_kappa,
                    realw_const_p pml_d, realw_p ax,realw_p ay, realw_p az)
{
  realw  pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl;

  const int idx = ispec * NGLL3 + k * NGLL2 + j * NGLLX + i;
  #define Rmat(i,j) r_trans[(idx * NDIM + j-1) * NDIM + i-1]
  #define Rmatinv(i,j) r_trans_inv[(idx * NDIM + j-1) * NDIM + i-1]
  #define TMAT(a,b,c)  temp[((c-1) * NDIM + (b-1)) * 6 + a-1]

  #define r1xl Rmat(1,1)
  #define r1yl Rmat(1,2)
  #define r1zl Rmat(1,3)
  #define r2xl Rmat(2,1)
  #define r2yl Rmat(2,2)
  #define r2zl Rmat(2,3)
  #define r3xl Rmat(3,1)
  #define r3yl Rmat(3,2)
  #define r3zl Rmat(3,3)


  #define ri1xl Rmatinv(1,1)
  #define ri1yl Rmatinv(1,2)
  #define ri1zl Rmatinv(1,3)
  #define ri2xl Rmatinv(2,1)
  #define ri2yl Rmatinv(2,2)
  #define ri2zl Rmatinv(2,3)
  #define ri3xl Rmatinv(3,1)
  #define ri3yl Rmatinv(3,2)
  #define ri3zl Rmatinv(3,3)

  realw k1l = 1./ pml_kappa[1-1 + idx * NDIM];
  realw k2l = 1./ pml_kappa[2-1 + idx * NDIM];
  realw k3l = 1./ pml_kappa[3-1 + idx * NDIM];

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
  #undef r1xl 
  #undef r1yl 
  #undef r1zl 
  #undef r2xl 
  #undef r2yl 
  #undef r2zl 
  #undef r3xl 
  #undef r3yl 
  #undef r3zl 
  #undef ri1xl
  #undef ri1yl
  #undef ri1zl
  #undef ri2xl
  #undef ri2yl
  #undef ri2zl
  #undef ri3xl
  #undef ri3yl
  #undef ri3zl
}

__device__ __forceinline__ void 
add_pml_physical_contribution(int ispec_pml,const int *phy_spec,
                              const int *phy_ijk,const int *ibool_CPML,
                              realw_const_p phy_norm,realw_const_p rtrans,
                              realw_const_p rtrans_inv,realw_const_p pml_d,
                              realw_const_p phy_jaco2Dw,realw sigma_xx,
                              realw sigma_xy,realw sigma_xz,realw sigma_yy,
                              realw sigma_yz,realw sigma_zz,
                              realw_p Qt_t)
{
  realw  d1l,d2l,d3l,r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,
  ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,
  weight, tx, ty, tz;

  for(int id = threadIdx.x; id < 6*NGLL2; id += blockDim.x) {
    int iside = id % 6;
    int igll = id / 6;
    if(iside >=6 || igll >=NGLL2) continue;
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
    sh_tempx[tx] = d_displ[iglob*NDIM];
    sh_tempy[tx] = d_displ[iglob*NDIM+1];
    sh_tempz[tx] = d_displ[iglob*NDIM+2];
    
    if(is_wavediscon)
      add_displ_discontinuity(working_element,tx,displ_wd,
                              ibool_wd,ispec_to_elem_wd,
                                &sh_tempx[tx],&sh_tempy[tx],
                                &sh_tempz[tx]);
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
    realw newtemp_adepml[3][3][6];
    int ispec_pml = spec_to_CPML[working_element] - 1;

    compute_sigma_adepml(ispec_pml,I,J,K,ANISOTROPY,d_c11store,d_c12store,d_c13store,
                        d_c14store,d_c15store,d_c16store,d_c22store,d_c23store,
                        d_c24store, d_c25store,d_c26store,d_c33store, d_c34store, 
                        d_c35store, d_c36store, d_c44store, d_c45store, d_c46store,
                        d_c55store, d_c56store, d_c66store,d_kappav,d_muv,offset,
                        duxdxl,duydxl,duzdxl,duxdyl,duydyl,duzdyl,duxdzl,duydzl,duzdzl,
                        rtrans,rtrans_inv,pml_kappa,&sigma_xx,&sigma_xy,&sigma_xz,
                        &sigma_yy,&sigma_yz,&sigma_zz,Qu
                      );

    if(threadIdx.x < NGLL3) {
      compute_Qu_t_point(ispec_pml,I,J,K,duxdxl, duydxl, duzdxl,
                          duxdyl, duydyl, duzdyl,
                          duxdzl, duydzl, duzdzl,rtrans_inv,
                        pml_d,Qu_t);
    }

    xixl *= jacobianl;
    xiyl *= jacobianl;
    xizl *= jacobianl;
    etaxl *= jacobianl;
    etayl *= jacobianl;
    etazl *= jacobianl;
    gammaxl *= jacobianl;
    gammayl *= jacobianl;
    gammazl *= jacobianl;

    // compute new_temp
    compute_newcomp(sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx,
                    sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yz,
                    sigma_zz,xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,
                  gammayl,gammazl,threadIdx.x,I,J,K,&newtemp_adepml[0][0][0]);

    // update Qt_t
    fac1 = d_wgllwgll_yz[K*NGLLX+J];
    fac2 = d_wgllwgll_xz[K*NGLLX+I];
    fac3 = d_wgllwgll_xy[J*NGLLX+I];
    compute_Qt_t(ispec_pml,I,J,K,threadIdx.x,&newtemp_adepml[0][0][0],
                  fac1,fac2,fac3,rtrans,rtrans_inv,
              pml_d,ibool_CPML,Qt_t);
    
    compute_accel_adepml(ispec_pml,I,J,K,&newtemp_adepml[0][0][0],
                          fac1,fac2,fac3,rtrans,rtrans_inv,
                         pml_kappa,pml_d,&sum_terms1,&sum_terms2,
                         &sum_terms3);

    // assembles acceleration array
    if (threadIdx.x < NGLL3) {
      atomicAdd(&d_accel[iglob*3], -sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], -sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], -sum_terms3);
        
    } // threadIdx.x
    
    // PML contribution
    add_pml_physical_contribution(ispec_pml,phy_spec,phy_ijk,
                            ibool_CPML,phy_norm,rtrans,rtrans_inv,
                          pml_d,phy_jaco2Dw,sigma_xx,sigma_xy,sigma_xz,
                          sigma_yy,sigma_yz,sigma_zz,Qt_t);
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