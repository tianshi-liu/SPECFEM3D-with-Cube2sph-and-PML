#include "mesh_constants_cuda.h"
/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2

/* ----------------------------------------------------------------------------------------------- */


// updates stress

__device__ __forceinline__ void compute_element_att_stress(int tx,int working_element,const int NSPEC,
                                           realw* R_xx,realw* R_yy,realw* R_xy,
                                           realw* R_xz,realw* R_yz,realw* Rxx_loc,realw* Ryy_loc,realw* Rxy_loc,
                                           realw* Rxz_loc,realw* Ryz_loc,
                                           realw* sigma_xx,realw* sigma_yy,realw* sigma_zz,
                                           realw* sigma_xy,realw* sigma_xz,realw* sigma_yz) {

  int offset_sls;
  realw rxx_sum,ryy_sum,rxy_sum,rxz_sum,ryz_sum;

  rxx_sum = 0.f;
  ryy_sum = 0.f;
  rxy_sum = 0.f;
  rxz_sum = 0.f;
  ryz_sum = 0.f;

  for(int i_sls = 0; i_sls < N_SLS; i_sls++){
    // index
    offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);
    //loads R_** values in local registers, which will be used later in the code
    Rxx_loc[i_sls] = get_global_cr( &R_xx[offset_sls] );
    Ryy_loc[i_sls] = get_global_cr( &R_yy[offset_sls] );
    Rxy_loc[i_sls] = get_global_cr( &R_xy[offset_sls] );
    Rxz_loc[i_sls] = get_global_cr( &R_xz[offset_sls] );
    Ryz_loc[i_sls] = get_global_cr( &R_yz[offset_sls] );
    rxx_sum += Rxx_loc[i_sls];
    ryy_sum += Ryy_loc[i_sls];
    rxy_sum += Rxy_loc[i_sls];
    rxz_sum += Rxz_loc[i_sls];
    ryz_sum += Ryz_loc[i_sls];
}
    *sigma_xx = *sigma_xx - rxx_sum;//Rxx_loc[i_sls];
    *sigma_yy = *sigma_yy - ryy_sum;//Ryy_loc[i_sls];
    *sigma_zz = *sigma_zz + rxx_sum + ryy_sum;//Rxx_loc[i_sls] + Ryy_loc[i_sls];
    *sigma_xy = *sigma_xy - rxy_sum;//Rxy_loc[i_sls];
    *sigma_xz = *sigma_xz - rxz_sum;//Rxz_loc[i_sls];
    *sigma_yz = *sigma_yz - ryz_sum;//Ryz_loc[i_sls];

  return;
}

/* ----------------------------------------------------------------------------------------------- */

// updates R_memory

__device__  __forceinline__ void compute_element_att_memory(int tx,int working_element,const int NSPEC,
                                          realw mul,
                                          realw_const_p factor_common,
                                          realw_const_p alphaval,realw_const_p betaval,realw_const_p gammaval,
                                          realw_p R_xx,realw_p R_yy,realw_p R_xy,realw_p R_xz,realw_p R_yz,
                                          realw* Rxx_loc,realw* Ryy_loc,realw* Rxy_loc,realw* Rxz_loc,realw* Ryz_loc,
                                          realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,
                                          realw_p epsilondev_xz,realw_p epsilondev_yz,
                                          realw epsilondev_xx_loc,realw epsilondev_yy_loc,realw epsilondev_xy_loc,
                                          realw epsilondev_xz_loc,realw epsilondev_yz_loc
                                          ){

  int ijk_ispec;
  int offset_sls,offset_common;
  realw alphaval_loc,betaval_loc,gammaval_loc;
  realw factor_loc;
  realw Sn_xx,Sn_yy,Sn_xy,Sn_xz,Sn_yz;

  // indices
  ijk_ispec = tx + NGLL3 * working_element;

//  mul = get_global_cr( &d_muv[offset_align] );
   Sn_xx   = get_global_cr( &epsilondev_xx[ijk_ispec] ); //(i,j,k,ispec)
   Sn_yy   = get_global_cr( &epsilondev_yy[ijk_ispec] );
   Sn_xy   = get_global_cr( &epsilondev_xy[ijk_ispec] );
   Sn_xz   = get_global_cr( &epsilondev_xz[ijk_ispec] );
   Sn_yz   = get_global_cr( &epsilondev_yz[ijk_ispec] );

  // use Runge-Kutta scheme to march in time
  for(int i_sls = 0; i_sls < N_SLS; i_sls++){

    // indices
    offset_common = i_sls + N_SLS*(tx + NGLL3*working_element); // (i_sls,i,j,k,ispec)
    offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);   // (i,j,k,ispec,i_sls)

    factor_loc = mul*get_global_cr( &factor_common[offset_common] ); //mustore(i,j,k,ispec) * factor_common(i_sls,i,j,k,ispec)

    alphaval_loc = alphaval[i_sls]; // (i_sls)
    betaval_loc = factor_loc*betaval[i_sls];
    gammaval_loc = factor_loc*gammaval[i_sls];

    R_xx[offset_sls] = alphaval_loc * Rxx_loc[i_sls] + betaval_loc * Sn_xx + gammaval_loc *  epsilondev_xx_loc;
    R_yy[offset_sls] = alphaval_loc * Ryy_loc[i_sls] + betaval_loc * Sn_yy + gammaval_loc *  epsilondev_yy_loc;
    R_xy[offset_sls] = alphaval_loc * Rxy_loc[i_sls] + betaval_loc * Sn_xy + gammaval_loc *  epsilondev_xy_loc;
    R_xz[offset_sls] = alphaval_loc * Rxz_loc[i_sls] + betaval_loc * Sn_xz + gammaval_loc *  epsilondev_xz_loc;
    R_yz[offset_sls] = alphaval_loc * Ryz_loc[i_sls] + betaval_loc * Sn_yz + gammaval_loc *  epsilondev_yz_loc;

  }
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// pre-computes gravity term

__device__ __forceinline__ void compute_element_gravity(int tx,int working_element,
                                        const int* iglob,
                                        realw_const_p d_minus_g,
                                        realw_const_p d_minus_deriv_gravity,
                                        realw_const_p d_rhostore,
                                        realw_const_p wgll_cube,
                                        realw jacobianl,
                                        realw* sh_displx,
                                        realw* sh_disply,
                                        realw* sh_displz,
                                        realw* sigma_xx,
                                        realw* sigma_yy,
                                        realw* sigma_xz,
                                        realw* sigma_yz,
                                        realw* rho_s_H1,
                                        realw* rho_s_H2,
                                        realw* rho_s_H3){

  realw minus_g,minus_dg;
  realw rhol;
  realw gzl; // gxl,gyl,
  realw sx_l,sy_l,sz_l;
  realw Hxxl,Hyyl,Hzzl; //,Hxyl,Hxzl,Hyzl;
  realw factor;

  // compute non-symmetric terms for gravity

  // get g, rho and dg/dr=dg
  minus_g = d_minus_g[*iglob];
  minus_dg = d_minus_deriv_gravity[*iglob];

  // Cartesian components of the gravitational acceleration
  //gxl = 0.f;
  //gyl = 0.f;
  gzl = minus_g;

  // Cartesian components of gradient of gravitational acceleration
  // H = grad g
  // assumes g only acts in negative z-direction
  Hxxl = 0.f;
  Hyyl = 0.f;
  Hzzl = minus_dg;
  //Hxyl = 0.f;
  //Hxzl = 0.f;
  //Hyzl = 0.f;

  rhol = get_global_cr( &d_rhostore[working_element*NGLL3_PADDED + tx] );

  // get displacement and multiply by density to compute G tensor
  // G = rho [ sg - (s * g) I  ]
  sx_l = rhol * sh_displx[tx]; // d_displ[iglob*3];
  sy_l = rhol * sh_disply[tx]; // d_displ[iglob*3 + 1];
  sz_l = rhol * sh_displz[tx]; // d_displ[iglob*3 + 2];

  // compute G tensor from s . g and add to sigma (not symmetric)
  //sigma_xx += sy_l*gyl + sz_l*gzl;
  *sigma_xx += sz_l*gzl;
  //sigma_yy += sx_l*gxl + sz_l*gzl;
  *sigma_yy += sz_l*gzl;
  //sigma_zz += sx_l*gxl + sy_l*gyl;

  //sigma_xy -= sx_l*gyl;
  //sigma_yx -= sy_l*gxl;

  *sigma_xz -= sx_l*gzl;
  //sigma_zx -= sz_l*gxl;

  *sigma_yz -= sy_l*gzl;
  //sigma_zy -= sz_l*gyl;

  // precompute vector
  factor = jacobianl * wgll_cube[tx];

  //rho_s_H1 = fac1 * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl);
  //rho_s_H2 = fac1 * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl);
  //rho_s_H3 = fac1 * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl);

  // only non-zero z-direction
  *rho_s_H1 = factor * sx_l * Hxxl ; // 0.f;
  *rho_s_H2 = factor * sy_l * Hyyl ; // 0.f;
  *rho_s_H3 = factor * sz_l * Hzzl ;

  // debug
  //*rho_s_H1 = 0.f;
  //*rho_s_H2 = 0.f;
  //*rho_s_H3 = 0.f ;

}

/* ----------------------------------------------------------------------------------------------- */

// loads displacement into shared memory for element

template<int FORWARD_OR_ADJOINT>
__device__  __forceinline__ void load_shared_memory_displ(const int* tx, const int* iglob,
                                                          realw_p d_displ,
                                                          realw* sh_displx,
                                                          realw* sh_disply,
                                                          realw* sh_displz){

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
#ifdef USE_TEXTURES_FIELDS
  sh_displx[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3);
  sh_disply[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 1);
  sh_displz[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 2);
#else
  // changing iglob indexing to match fortran row changes fast style
  sh_displx[(*tx)] = d_displ[(*iglob)*3];
  sh_disply[(*tx)] = d_displ[(*iglob)*3 + 1];
  sh_displz[(*tx)] = d_displ[(*iglob)*3 + 2];
#endif

}

/* ----------------------------------------------------------------------------------------------- */

// loads displacement + viscosity * velocity into shared memory for element

template<int FORWARD_OR_ADJOINT>
__device__  __forceinline__ void load_shared_memory_displ_visco(const int* tx, const int* iglob,
                                                          realw_p d_displ,
                                                          realw_const_p d_veloc,
                                                          realw visco,
                                                          realw* sh_displx,
                                                          realw* sh_disply,
                                                          realw* sh_displz){

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
#ifdef USE_TEXTURES_FIELDS
  sh_displx[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3) + visco*texfetch_veloc<FORWARD_OR_ADJOINT>((*iglob)*3);
  sh_disply[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 1) + visco*texfetch_veloc<FORWARD_OR_ADJOINT>((*iglob)*3 + 1);
  sh_displz[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*3 + 2) + visco*texfetch_veloc<FORWARD_OR_ADJOINT>((*iglob)*3 + 2);
#else
  // changing iglob indexing to match fortran row changes fast style
  sh_displx[(*tx)] = d_displ[(*iglob)*3] + visco * d_veloc[(*iglob)*3];
  sh_disply[(*tx)] = d_displ[(*iglob)*3 + 1] + visco * d_veloc[(*iglob)*3 + 1];
  sh_displz[(*tx)] = d_displ[(*iglob)*3 + 2] + visco * d_veloc[(*iglob)*3 + 2];
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// loads hprime into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprime(const int* tx,
                                                           realw_const_p d_hprime_xx,
                                                           realw* sh_hprime_xx){

  // each thread reads its corresponding value
  // (might be faster sometimes...)
#ifdef USE_TEXTURES_CONSTANTS
  // hprime
  sh_hprime_xx[(*tx)] = tex1Dfetch(d_hprime_xx_tex,tx + d_hprime_xx_tex_offset);
#else
  // hprime
  sh_hprime_xx[(*tx)] = d_hprime_xx[(*tx)];
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// loads hprimewgll into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprimewgll(const int* tx,
                                                               realw_const_p d_hprimewgll_xx,
                                                               realw* sh_hprimewgll_xx ){

  // each thread reads its corresponding value
  // weighted hprime
  sh_hprimewgll_xx[(*tx)] = d_hprimewgll_xx[(*tx)];
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprime_xi(int I, int J, int K,
                                              realw* tempxl,realw* tempyl,realw* tempzl,
                                              realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprime ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 1. cut-plane along xi-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+I];

    sumx += sh_tempx[K*NGLL2+J*NGLLX+l] * fac;
    sumy += sh_tempy[K*NGLL2+J*NGLLX+l] * fac;
    sumz += sh_tempz[K*NGLL2+J*NGLLX+l] * fac;
  }

// counts:
// + NGLLX * ( 2 + 3*6) FLOP = 100 FLOP
//
// + 0 BYTE

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprime_eta(int I, int J, int K,
                                               realw* tempxl,realw* tempyl,realw* tempzl,
                                               realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprime ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 2. cut-plane along eta-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+J];

    sumx += sh_tempx[K*NGLL2+l*NGLLX+I] * fac;
    sumy += sh_tempy[K*NGLL2+l*NGLLX+I] * fac;
    sumz += sh_tempz[K*NGLL2+l*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprime_gamma(int I, int J, int K,
                                                 realw* tempxl,realw* tempyl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprime ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 3. cut-plane along gamma-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+K];

    sumx += sh_tempx[l*NGLL2+J*NGLLX+I] * fac;
    sumy += sh_tempy[l*NGLL2+J*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLL2+J*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprimewgll_xi(int I, int J, int K,
                                                   realw* tempxl,realw* tempyl,realw* tempzl,
                                                   realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprimewgll ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 1. cut-plane along xi-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[I*NGLLX+l]; //  d_hprimewgll_xx[I*NGLLX+l];

    sumx += sh_tempx[K*NGLL2+J*NGLLX+l] * fac;
    sumy += sh_tempy[K*NGLL2+J*NGLLX+l] * fac;
    sumz += sh_tempz[K*NGLL2+J*NGLLX+l] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprimewgll_eta(int I, int J, int K,
                                               realw* tempxl,realw* tempyl,realw* tempzl,
                                               realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprimewgll ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 2. cut-plane along eta-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[J*NGLLX+l]; // d_hprimewgll_xx[J*NGLLX+l];

    sumx += sh_tempx[K*NGLL2+l*NGLLX+I] * fac;
    sumy += sh_tempy[K*NGLL2+l*NGLLX+I] * fac;
    sumz += sh_tempz[K*NGLL2+l*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprimewgll_gamma(int I, int J, int K,
                                                 realw* tempxl,realw* tempyl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempy,realw* sh_tempz, realw* sh_hprimewgll ){

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumy = 0.f;
  realw sumz = 0.f;

  // 3. cut-plane along gamma-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[K*NGLLX+l]; // d_hprimewgll_xx[K*NGLLX+l];

    sumx += sh_tempx[l*NGLL2+J*NGLLX+I] * fac;
    sumy += sh_tempy[l*NGLL2+J*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLL2+J*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempyl = sumy;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */

// computes the spatial derivatives


__device__  __forceinline__ void
  get_spatial_derivatives(realw* xixl,realw* xiyl,realw* xizl,realw* etaxl,realw* etayl,realw* etazl,
                          realw* gammaxl,realw* gammayl,realw* gammazl,realw* jacobianl,int I,int J,int K,int tx,
                          realw* tempx1l,realw* tempy1l,realw* tempz1l,realw* tempx2l,realw* tempy2l,realw* tempz2l,
                          realw* tempx3l,realw* tempy3l,realw* tempz3l,realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,realw* sh_hprime_xx,
                          realw* duxdxl,realw* duxdyl,realw* duxdzl,realw* duydxl,realw* duydyl,realw* duydzl,realw* duzdxl,realw* duzdyl,realw* duzdzl,
                          realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                          int ispec_irreg, realw xix_regular,int ipass){

  // computes first matrix products
  if (ispec_irreg >= 0 && ipass==0 ){ //irregular_element
    // local padded index
    int offset = ispec_irreg*NGLL3_PADDED + tx;

    *xixl = get_global_cr( &d_xix[offset]);
    *xiyl = get_global_cr(&d_xiy[offset]);
    *xizl = get_global_cr(&d_xiz[offset]);
    *etaxl = get_global_cr(&d_etax[offset]);
    *etayl = get_global_cr(&d_etay[offset]);
    *etazl = get_global_cr(&d_etaz[offset]);
    *gammaxl = get_global_cr(&d_gammax[offset]);
    *gammayl = get_global_cr(&d_gammay[offset]);
    *gammazl = get_global_cr(&d_gammaz[offset]);

    *jacobianl = 1.f / ((*xixl)*((*etayl)*(*gammazl)-(*etazl)*(*gammayl))
                      -(*xiyl)*((*etaxl)*(*gammazl)-(*etazl)*(*gammaxl))
                      +(*xizl)*((*etaxl)*(*gammayl)-(*etayl)*(*gammaxl)));
  }

  // 1. cut-plane
  sum_hprime_xi(I,J,K,tempx1l,tempy1l,tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 2. cut-plane
  sum_hprime_eta(I,J,K,tempx2l,tempy2l,tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,K,tempx3l,tempy3l,tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);


    // synchronize all the threads (one thread for each of the NGLL grid points of the
    // current spectral element) because we need the whole element to be ready in order
    // to be able to compute the matrix products along cut planes of the 3D element below
    __syncthreads();

  if (ispec_irreg >= 0 ){ //irregular_element

    // compute derivatives of ux, uy and uz with respect to x, y and z
    (*duxdxl) = (*xixl)*(*tempx1l) + (*etaxl)*(*tempx2l) + (*gammaxl)*(*tempx3l);
    (*duxdyl) = (*xiyl)*(*tempx1l) + (*etayl)*(*tempx2l) + (*gammayl)*(*tempx3l);
    (*duxdzl) = (*xizl)*(*tempx1l) + (*etazl)*(*tempx2l) + (*gammazl)*(*tempx3l);

    (*duydxl) = (*xixl)*(*tempy1l) + (*etaxl)*(*tempy2l) + (*gammaxl)*(*tempy3l);
    (*duydyl) = (*xiyl)*(*tempy1l) + (*etayl)*(*tempy2l) + (*gammayl)*(*tempy3l);
    (*duydzl) = (*xizl)*(*tempy1l) + (*etazl)*(*tempy2l) + (*gammazl)*(*tempy3l);

    (*duzdxl) = (*xixl)*(*tempz1l) + (*etaxl)*(*tempz2l) + (*gammaxl)*(*tempz3l);
    (*duzdyl) = (*xiyl)*(*tempz1l) + (*etayl)*(*tempz2l) + (*gammayl)*(*tempz3l);
    (*duzdzl) = (*xizl)*(*tempz1l) + (*etazl)*(*tempz2l) + (*gammazl)*(*tempz3l);
  }
  else{
    // compute derivatives of ux, uy and uz with respect to x, y and z
    (*duxdxl) = xix_regular*(*tempx1l);
    (*duxdyl) = xix_regular*(*tempx2l);
    (*duxdzl) = xix_regular*(*tempx3l);

    (*duydxl) = xix_regular*(*tempy1l);
    (*duydyl) = xix_regular*(*tempy2l);
    (*duydzl) = xix_regular*(*tempy3l);

    (*duzdxl) = xix_regular*(*tempz1l);
    (*duzdyl) = xix_regular*(*tempz2l);
    (*duzdzl) = xix_regular*(*tempz3l);
  }
  // counts:
  // + 9 * 5 FLOP = 45 FLOP
  //
  // + 0 BYTE


}

/* ----------------------------------------------------------------------------------------------- */
// computes dot product between tensor and derivatives

__device__  __forceinline__ void
  get_dot_product(realw jacobianl,realw sigma_xx,realw sigma_xy,realw sigma_yx,realw sigma_xz,realw sigma_zx,realw sigma_yy,realw sigma_yz,realw sigma_zy,realw sigma_zz,
                  realw Dxl,realw Dyl,realw Dzl,realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,int tx,
                  int ispec_irreg,realw xix_regular,realw jacobian_regular,int component){

  // fills shared memory arrays
  if (threadIdx.x < NGLL3) {

    if (ispec_irreg>=0){ //irregular element
      sh_tempx[tx] = jacobianl * (sigma_xx*Dxl + sigma_yx*Dyl + sigma_zx*Dzl); // sh_tempx1
      sh_tempy[tx] = jacobianl * (sigma_xy*Dxl + sigma_yy*Dyl + sigma_zy*Dzl); // sh_tempy1
      sh_tempz[tx] = jacobianl * (sigma_xz*Dxl + sigma_yz*Dyl + sigma_zz*Dzl); // sh_tempz1
    }
    else if (component==1){
      sh_tempx[tx] = jacobian_regular * (sigma_xx*xix_regular); // sh_tempx1
      sh_tempy[tx] = jacobian_regular * (sigma_xy*xix_regular); // sh_tempy1
      sh_tempz[tx] = jacobian_regular * (sigma_xz*xix_regular); // sh_tempz1
    }
    else if (component==2){
      sh_tempx[tx] = jacobian_regular * (sigma_yx*xix_regular); // sh_tempx1
      sh_tempy[tx] = jacobian_regular * (sigma_yy*xix_regular); // sh_tempy1
      sh_tempz[tx] = jacobian_regular * (sigma_yz*xix_regular); // sh_tempz1

    }else{
      sh_tempx[tx] = jacobian_regular * (sigma_zx*xix_regular); // sh_tempx1
      sh_tempy[tx] = jacobian_regular * (sigma_zy*xix_regular); // sh_tempy1
      sh_tempz[tx] = jacobian_regular * (sigma_zz*xix_regular); // sh_tempz1
    }
  }
  __syncthreads();
}

