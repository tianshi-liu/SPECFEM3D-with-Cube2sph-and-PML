subroutine compute_sigma_adepml(ispec,i,j,k,lambdalplus2mul,lambdal,mul,&
                                dsigma,&
                                duxdxl, duydxl, duzdxl, &
                                duxdyl, duydyl, duzdyl, &
                                duxdzl, duydzl, duzdzl)
  use constants, only: CUSTOM_REAL
  use pml_par, only: Qu,r_trans, r_trans_inv, pml_kappa
  implicit none
  real(kind=CUSTOM_REAL) :: k1l,k2l,k3l,&
          r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,&
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,&
          pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl,&
          Qu1xl,Qu1yl,Qu1zl,Qu2xl,Qu2yl,Qu2zl,Qu3xl,Qu3yl,Qu3zl,&
          e11,e12,e13,e22,e23,e33
  real(kind=CUSTOM_REAL), intent(in) :: lambdalplus2mul,lambdal,mul, &
                                duxdxl, duydxl, duzdxl, &
                                duxdyl, duydyl, duzdyl, &
                                duxdzl, duydzl, duzdzl
  real(kind=CUSTOM_REAL), intent(out) :: dsigma(6)
  integer, intent(in) :: ispec,i,j,k
  r1xl = r_trans(1,1,i,j,k,ispec)
  r2xl = r_trans(2,1,i,j,k,ispec)
  r3xl = r_trans(3,1,i,j,k,ispec)
  r1yl = r_trans(1,2,i,j,k,ispec)
  r2yl = r_trans(2,2,i,j,k,ispec)
  r3yl = r_trans(3,2,i,j,k,ispec)
  r1zl = r_trans(1,3,i,j,k,ispec)
  r2zl = r_trans(2,3,i,j,k,ispec)
  r3zl = r_trans(3,3,i,j,k,ispec)

  ri1xl = r_trans_inv(1,1,i,j,k,ispec)
  ri2xl = r_trans_inv(2,1,i,j,k,ispec)
  ri3xl = r_trans_inv(3,1,i,j,k,ispec)
  ri1yl = r_trans_inv(1,2,i,j,k,ispec)
  ri2yl = r_trans_inv(2,2,i,j,k,ispec)
  ri3yl = r_trans_inv(3,2,i,j,k,ispec)
  ri1zl = r_trans_inv(1,3,i,j,k,ispec)
  ri2zl = r_trans_inv(2,3,i,j,k,ispec)
  ri3zl = r_trans_inv(3,3,i,j,k,ispec)

  k1l = pml_kappa(1,i,j,k,ispec)
  k2l = pml_kappa(2,i,j,k,ispec)
  k3l = pml_kappa(3,i,j,k,ispec)
  !b1l = pml_beta(1,i,j,k,ispec)
  !b2l = pml_beta(2,i,j,k,ispec)
  !b3l = pml_beta(3,i,j,k,ispec)
  !d1l = pml_d(1,i,j,k,ispec)
  !d2l = pml_d(2,i,j,k,ispec)
  !d3l = pml_d(3,i,j,k,ispec)

  pmlxxl = r1xl*ri1xl/k1l + r2xl*ri2xl/k2l + r3xl*ri3xl/k3l
  pmlxyl = r1xl*ri1yl/k1l + r2xl*ri2yl/k2l + r3xl*ri3yl/k3l
  pmlxzl = r1xl*ri1zl/k1l + r2xl*ri2zl/k2l + r3xl*ri3zl/k3l
  pmlyxl = r1yl*ri1xl/k1l + r2yl*ri2xl/k2l + r3yl*ri3xl/k3l
  pmlyyl = r1yl*ri1yl/k1l + r2yl*ri2yl/k2l + r3yl*ri3yl/k3l
  pmlyzl = r1yl*ri1zl/k1l + r2yl*ri2zl/k2l + r3yl*ri3zl/k3l
  pmlzxl = r1zl*ri1xl/k1l + r2zl*ri2xl/k2l + r3zl*ri3xl/k3l
  pmlzyl = r1zl*ri1yl/k1l + r2zl*ri2yl/k2l + r3zl*ri3yl/k3l
  pmlzzl = r1zl*ri1zl/k1l + r2zl*ri2zl/k2l + r3zl*ri3zl/k3l

  Qu1xl = Qu(1,1,i,j,k,ispec)
  Qu2xl = Qu(2,1,i,j,k,ispec)
  Qu3xl = Qu(3,1,i,j,k,ispec)
  Qu1yl = Qu(1,2,i,j,k,ispec)
  Qu2yl = Qu(2,2,i,j,k,ispec)
  Qu3yl = Qu(3,2,i,j,k,ispec)
  Qu1zl = Qu(1,3,i,j,k,ispec)
  Qu2zl = Qu(2,3,i,j,k,ispec)
  Qu3zl = Qu(3,3,i,j,k,ispec)

  e11 = pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
               r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl
  e22 = pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
               r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl
  e33 = pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
               r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl
  e12 = pmlxxl*duydxl+pmlxyl*duydyl+pmlxzl*duydzl+ &
               r1xl*Qu1yl+r2xl*Qu2yl+r3xl*Qu3yl+ &
        pmlyxl*duxdxl+pmlyyl*duxdyl+pmlyzl*duxdzl+ &
               r1yl*Qu1xl+r2yl*Qu2xl+r3yl*Qu3xl
  e13 = pmlxxl*duzdxl+pmlxyl*duzdyl+pmlxzl*duzdzl+ &
               r1xl*Qu1zl+r2xl*Qu2zl+r3xl*Qu3zl+ &
        pmlzxl*duxdxl+pmlzyl*duxdyl+pmlzzl*duxdzl+ &
               r1zl*Qu1xl+r2zl*Qu2xl+r3zl*Qu3xl
  e23 = pmlyxl*duzdxl+pmlyyl*duzdyl+pmlyzl*duzdzl+ &
               r1yl*Qu1zl+r2yl*Qu2zl+r3yl*Qu3zl+ &
        pmlzxl*duydxl+pmlzyl*duydyl+pmlzzl*duydzl+ &
               r1zl*Qu1yl+r2zl*Qu2yl+r3zl*Qu3yl
  dsigma(1) = lambdalplus2mul*e11+lambdal*(e22+e33)
  dsigma(2) = mul*e12
  dsigma(3) = mul*e13
  dsigma(4) = lambdalplus2mul*e22+lambdal*(e11+e33)
  dsigma(5) = mul*e23
  dsigma(6) = lambdalplus2mul*e33+lambdal*(e11+e22)
  !dsigma(1) = lambdalplus2mul*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
  !         r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
  !           lambdal*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
  !         r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
  !           lambdal*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
  !         r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)
  !dsigma(2) = mul*(pmlxxl*duydxl+pmlxyl*duydyl+pmlxzl*duydzl+ &
  !         r1xl*Qu1yl+r2xl*Qu2yl+r3xl*Qu3yl) + &
  !           mul*(pmlyxl*duxdxl+pmlyyl*duxdyl+pmlyzl*duxdzl+ &
  !         r1yl*Qu1xl+r2yl*Qu2xl+r3yl*Qu3xl)
  !dsigma(3) = mul*(pmlxxl*duzdxl+pmlxyl*duzdyl+pmlxzl*duzdzl+ &
  !         r1xl*Qu1zl+r2xl*Qu2zl+r3xl*Qu3zl) + &
  !           mul*(pmlzxl*duxdxl+pmlzyl*duxdyl+pmlzzl*duxdzl+ &
  !         r1zl*Qu1xl+r2zl*Qu2xl+r3zl*Qu3xl)
  !dsigma(4) = lambdal*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
  !         r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
  !           lambdalplus2mul*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
  !         r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
  !           lambdal*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
  !         r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)
  !dsigma(5) = mul*(pmlyxl*duzdxl+pmlyyl*duzdyl+pmlyzl*duzdzl+ &
  !         r1yl*Qu1zl+r2yl*Qu2zl+r3yl*Qu3zl) + &
  !           mul*(pmlzxl*duydxl+pmlzyl*duydyl+pmlzzl*duydzl+ &
  !         r1zl*Qu1yl+r2zl*Qu2yl+r3zl*Qu3yl)
  !dsigma(6) = lambdal*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
  !         r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
  !           lambdal*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
  !         r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
  !           lambdalplus2mul*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
  !         r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)

end subroutine compute_sigma_adepml

subroutine compute_Qu_t(veloc)
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD,FOUR_THIRDS
  use specfem_par, only: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
               NGLOB_AB,hprime_xxT,hprime_yyT,hprime_zzT,ibool, &
               irregular_element_number,xix_regular
  use pml_par, only: CPML_to_spec,nspec_CPML
  implicit none
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: veloc
  ! local parameters
  integer :: ispec, iglob, ispec_irreg
  integer :: i,j,k
  integer :: ispec_CPML
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  !real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  !real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  !real(kind=CUSTOM_REAL) :: kappal
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            dummyx_loc,dummyy_loc,dummyz_loc, &
            tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: zero_array
  ! generate an array equal to zero
  zero_array(:,:,:) = 0._CUSTOM_REAL
  do ispec_CPML = 1, nspec_CPML
    ispec = CPML_to_spec(ispec_CPML)
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dummyx_loc(i,j,k) = veloc(1,iglob)
          dummyy_loc(i,j,k) = veloc(2,iglob)
          dummyz_loc(i,j,k) = veloc(3,iglob)
        enddo
      enddo
    enddo

    call compute_strain_in_element( &
                 tempx1,tempx2,tempx3,zero_array,zero_array,zero_array, &
                 tempy1,tempy2,tempy3,zero_array,zero_array,zero_array, &
                 tempz1,tempz2,tempz3,zero_array,zero_array,zero_array, &
                 dummyx_loc,dummyy_loc,dummyz_loc, &
                 hprime_xxT,hprime_yyT,hprime_zzT)
    ispec_irreg = irregular_element_number(ispec)
    !if (ispec_irreg == 0) jacobianl = jacobian_regular
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          if (ispec_irreg /= 0) then !irregular element

            xixl = xix(i,j,k,ispec_irreg)
            xiyl = xiy(i,j,k,ispec_irreg)
            xizl = xiz(i,j,k,ispec_irreg)
            etaxl = etax(i,j,k,ispec_irreg)
            etayl = etay(i,j,k,ispec_irreg)
            etazl = etaz(i,j,k,ispec_irreg)
            gammaxl = gammax(i,j,k,ispec_irreg)
            gammayl = gammay(i,j,k,ispec_irreg)
            gammazl = gammaz(i,j,k,ispec_irreg)
            !jacobianl = jacobian(i,j,k,ispec_irreg)

            duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
            duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
            duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

            duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
            duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
            duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

            duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
            duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
            duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

          else !regular element   
            duxdxl = xix_regular*tempx1(i,j,k)
            duxdyl = xix_regular*tempx2(i,j,k)
            duxdzl = xix_regular*tempx3(i,j,k)

            duydxl = xix_regular*tempy1(i,j,k)
            duydyl = xix_regular*tempy2(i,j,k)
            duydzl = xix_regular*tempy3(i,j,k)

            duzdxl = xix_regular*tempz1(i,j,k)
            duzdyl = xix_regular*tempz2(i,j,k)
            duzdzl = xix_regular*tempz3(i,j,k)

          endif

          call compute_Qu_t_point(ispec_CPML,i,j,k,&
                                  duxdxl, duydxl, duzdxl, &
                                  duxdyl, duydyl, duzdyl, &
                                  duxdzl, duydzl, duzdzl)
        enddo
      enddo
    enddo  
  enddo

  contains
   subroutine compute_strain_in_element(tempx1_att,tempx2_att,tempx3_att,tempx1,tempx2,tempx3, &
                                            tempy1_att,tempy2_att,tempy3_att,tempy1,tempy2,tempy3, &
                                            tempz1_att,tempz2_att,tempz3_att,tempz1,tempz2,tempz3, &
                                            dummyx_loc,dummyy_loc,dummyz_loc,hprime_xxT,hprime_yyT,hprime_zzT)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_att,tempx2_att,tempx3_att, &
                                                          tempy1_att,tempy2_att,tempy3_att, &
                                                          tempz1_att,tempz2_att,tempz3_att

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3,&
                                                          tempy1,tempy2,tempy3,&
                                                          tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yyT
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zzT

  ! local variables
  integer :: i,j,k,l
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  tempx1_att(:,:,:) = tempx1(:,:,:)
  tempx2_att(:,:,:) = tempx2(:,:,:)
  tempx3_att(:,:,:) = tempx3(:,:,:)

  tempy1_att(:,:,:) = tempy1(:,:,:)
  tempy2_att(:,:,:) = tempy2(:,:,:)
  tempy3_att(:,:,:) = tempy3(:,:,:)

  tempz1_att(:,:,:) = tempz1(:,:,:)
  tempz2_att(:,:,:) = tempz2(:,:,:)
  tempz3_att(:,:,:) = tempz3(:,:,:)

  ! use first order Taylor expansion of displacement for local storage of
  ! stresses
  ! at this current time step, to fix attenuation in a consistent way
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

        ! we can merge these loops because NGLLX = NGLLY = NGLLZ
        do l=1,NGLLX
          hp1 = hprime_xxT(l,i)
          tempx1_att(i,j,k) = tempx1_att(i,j,k) + dummyx_loc(l,j,k) * hp1
          tempy1_att(i,j,k) = tempy1_att(i,j,k) + dummyy_loc(l,j,k) * hp1
          tempz1_att(i,j,k) = tempz1_att(i,j,k) + dummyz_loc(l,j,k) * hp1

          hp2 = hprime_yyT(l,j)
          tempx2_att(i,j,k) = tempx2_att(i,j,k) + dummyx_loc(i,l,k) * hp2
          tempy2_att(i,j,k) = tempy2_att(i,j,k) + dummyy_loc(i,l,k) * hp2
          tempz2_att(i,j,k) = tempz2_att(i,j,k) + dummyz_loc(i,l,k) * hp2

          hp3 = hprime_zzT(l,k)
          tempx3_att(i,j,k) = tempx3_att(i,j,k) + dummyx_loc(i,j,l) * hp3
          tempy3_att(i,j,k) = tempy3_att(i,j,k) + dummyy_loc(i,j,l) * hp3
          tempz3_att(i,j,k) = tempz3_att(i,j,k) + dummyz_loc(i,j,l) * hp3
        enddo

      enddo
    enddo
  enddo

  end subroutine compute_strain_in_element
end subroutine compute_Qu_t

subroutine compute_Qu_t_point(ispec,i,j,k,&
                              duxdxl, duydxl, duzdxl, &
                              duxdyl, duydyl, duzdyl, &
                              duxdzl, duydzl, duzdzl)
  use constants, only: CUSTOM_REAL
  use pml_par, only: Qu_t, r_trans_inv, pml_d
  implicit none
  real(kind=CUSTOM_REAL) :: d1l,d2l,d3l,&
          !r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,&
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl
  real(kind=CUSTOM_REAL), intent(in) :: &
                                duxdxl, duydxl, duzdxl, &
                                duxdyl, duydyl, duzdyl, &
                                duxdzl, duydzl, duzdzl
  integer, intent(in) :: ispec,i,j,k
  !r1xl = r_trans(1,1,i,j,k,ispec)
  !r2xl = r_trans(2,1,i,j,k,ispec)
  !r3xl = r_trans(3,1,i,j,k,ispec)
  !r1yl = r_trans(1,2,i,j,k,ispec)
  !r2yl = r_trans(2,2,i,j,k,ispec)
  !r3yl = r_trans(3,2,i,j,k,ispec)
  !r1zl = r_trans(1,3,i,j,k,ispec)
  !r2zl = r_trans(2,3,i,j,k,ispec)
  !r3zl = r_trans(3,3,i,j,k,ispec)

  ri1xl = r_trans_inv(1,1,i,j,k,ispec)
  ri2xl = r_trans_inv(2,1,i,j,k,ispec)
  ri3xl = r_trans_inv(3,1,i,j,k,ispec)
  ri1yl = r_trans_inv(1,2,i,j,k,ispec)
  ri2yl = r_trans_inv(2,2,i,j,k,ispec)
  ri3yl = r_trans_inv(3,2,i,j,k,ispec)
  ri1zl = r_trans_inv(1,3,i,j,k,ispec)
  ri2zl = r_trans_inv(2,3,i,j,k,ispec)
  ri3zl = r_trans_inv(3,3,i,j,k,ispec)

  !k1l = pml_kappa(1,i,j,k,ispec)
  !k2l = pml_kappa(2,i,j,k,ispec)
  !k3l = pml_kappa(3,i,j,k,ispec)
  !b1l = pml_beta(1,i,j,k,ispec)
  !b2l = pml_beta(2,i,j,k,ispec)
  !b3l = pml_beta(3,i,j,k,ispec)
  d1l = pml_d(1,i,j,k,ispec)
  d2l = pml_d(2,i,j,k,ispec)
  d3l = pml_d(3,i,j,k,ispec)

  Qu_t(1,1,i,j,k,ispec)=-d1l*(ri1xl*duxdxl+ri1yl*duxdyl+ri1zl*duxdzl)
  Qu_t(1,2,i,j,k,ispec)=-d1l*(ri1xl*duydxl+ri1yl*duydyl+ri1zl*duydzl)
  Qu_t(1,3,i,j,k,ispec)=-d1l*(ri1xl*duzdxl+ri1yl*duzdyl+ri1zl*duzdzl)
  Qu_t(2,1,i,j,k,ispec)=-d2l*(ri2xl*duxdxl+ri2yl*duxdyl+ri2zl*duxdzl)
  Qu_t(2,2,i,j,k,ispec)=-d2l*(ri2xl*duydxl+ri2yl*duydyl+ri2zl*duydzl)
  Qu_t(2,3,i,j,k,ispec)=-d2l*(ri2xl*duzdxl+ri2yl*duzdyl+ri2zl*duzdzl)
  Qu_t(3,1,i,j,k,ispec)=-d3l*(ri3xl*duxdxl+ri3yl*duxdyl+ri3zl*duxdzl)
  Qu_t(3,2,i,j,k,ispec)=-d3l*(ri3xl*duydxl+ri3yl*duydyl+ri3zl*duydzl)
  Qu_t(3,3,i,j,k,ispec)=-d3l*(ri3xl*duzdxl+ri3yl*duzdyl+ri3zl*duzdzl)
end subroutine compute_Qu_t_point


subroutine compute_Qt_t(ispec,i,j,k,temp,fac1,fac2,fac3)
  use constants, only: CUSTOM_REAL
  use pml_par, only: Qt_t, r_trans, r_trans_inv, &
          pml_d, ibool_CPML
  implicit none
  !real(kind=CUSTOM_REAL), dimension(3,3), intent(out) :: Qt_tl
  real(kind=CUSTOM_REAL) :: d1l,d2l,d3l,&
          r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,&
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,&
          pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl
  integer :: ispec,i,j,k
  real(kind=CUSTOM_REAL), intent(in) :: fac1,fac2,fac3
  real(kind=CUSTOM_REAL), dimension(6,3,3), intent(in) :: temp
  integer :: iglob
  iglob = ibool_CPML(i,j,k,ispec)
  r1xl = r_trans(1,1,i,j,k,ispec)
  r2xl = r_trans(2,1,i,j,k,ispec)
  r3xl = r_trans(3,1,i,j,k,ispec)
  r1yl = r_trans(1,2,i,j,k,ispec)
  r2yl = r_trans(2,2,i,j,k,ispec)
  r3yl = r_trans(3,2,i,j,k,ispec)
  r1zl = r_trans(1,3,i,j,k,ispec)
  r2zl = r_trans(2,3,i,j,k,ispec)
  r3zl = r_trans(3,3,i,j,k,ispec)

  ri1xl = r_trans_inv(1,1,i,j,k,ispec)
  ri2xl = r_trans_inv(2,1,i,j,k,ispec)
  ri3xl = r_trans_inv(3,1,i,j,k,ispec)
  ri1yl = r_trans_inv(1,2,i,j,k,ispec)
  ri2yl = r_trans_inv(2,2,i,j,k,ispec)
  ri3yl = r_trans_inv(3,2,i,j,k,ispec)
  ri1zl = r_trans_inv(1,3,i,j,k,ispec)
  ri2zl = r_trans_inv(2,3,i,j,k,ispec)
  ri3zl = r_trans_inv(3,3,i,j,k,ispec)

  !k1l = pml_kappa(1,i,j,k,ispec)
  !k2l = pml_kappa(2,i,j,k,ispec)
  !k3l = pml_kappa(3,i,j,k,ispec)
  !b1l = pml_beta(1,i,j,k,ispec)
  !b2l = pml_beta(2,i,j,k,ispec)
  !b3l = pml_beta(3,i,j,k,ispec)
  d1l = pml_d(1,i,j,k,ispec)
  d2l = pml_d(2,i,j,k,ispec)
  d3l = pml_d(3,i,j,k,ispec)

  pmlxxl = r1xl*ri1xl
  pmlxyl = r1xl*ri1yl
  pmlxzl = r1xl*ri1zl
  pmlyxl = r1yl*ri1xl
  pmlyyl = r1yl*ri1yl
  pmlyzl = r1yl*ri1zl
  pmlzxl = r1zl*ri1xl
  pmlzyl = r1zl*ri1yl
  pmlzzl = r1zl*ri1zl

  Qt_t(1,1,iglob) = Qt_t(1,1,iglob)+&
     d1l*(fac1*(temp(1,1,1)*pmlxxl+temp(1,2,1)*pmlxyl+temp(1,3,1)*pmlxzl+&
                temp(2,1,1)*pmlyxl+temp(2,2,1)*pmlyyl+temp(2,3,1)*pmlyzl+&
                temp(3,1,1)*pmlzxl+temp(3,2,1)*pmlzyl+temp(3,3,1)*pmlzzl) +&
          fac2*(temp(1,1,2)*pmlxxl+temp(1,2,2)*pmlxyl+temp(1,3,2)*pmlxzl+&
                temp(2,1,2)*pmlyxl+temp(2,2,2)*pmlyyl+temp(2,3,2)*pmlyzl+&
                temp(3,1,2)*pmlzxl+temp(3,2,2)*pmlzyl+temp(3,3,2)*pmlzzl) +&
          fac3*(temp(1,1,3)*pmlxxl+temp(1,2,3)*pmlxyl+temp(1,3,3)*pmlxzl+&
                temp(2,1,3)*pmlyxl+temp(2,2,3)*pmlyyl+temp(2,3,3)*pmlyzl+&
                temp(3,1,3)*pmlzxl+temp(3,2,3)*pmlzyl+temp(3,3,3)*pmlzzl))
  Qt_t(1,2,iglob) = Qt_t(1,2,iglob)+&
     d1l*(fac1*(temp(2,1,1)*pmlxxl+temp(2,2,1)*pmlxyl+temp(2,3,1)*pmlxzl+&
                temp(4,1,1)*pmlyxl+temp(4,2,1)*pmlyyl+temp(4,3,1)*pmlyzl+&
                temp(5,1,1)*pmlzxl+temp(5,2,1)*pmlzyl+temp(5,3,1)*pmlzzl) +&
          fac2*(temp(2,1,2)*pmlxxl+temp(2,2,2)*pmlxyl+temp(2,3,2)*pmlxzl+&
                temp(4,1,2)*pmlyxl+temp(4,2,2)*pmlyyl+temp(4,3,2)*pmlyzl+&
                temp(5,1,2)*pmlzxl+temp(5,2,2)*pmlzyl+temp(5,3,2)*pmlzzl) +&
          fac3*(temp(2,1,3)*pmlxxl+temp(2,2,3)*pmlxyl+temp(2,3,3)*pmlxzl+&
                temp(4,1,3)*pmlyxl+temp(4,2,3)*pmlyyl+temp(4,3,3)*pmlyzl+&
                temp(5,1,3)*pmlzxl+temp(5,2,3)*pmlzyl+temp(5,3,3)*pmlzzl))
  Qt_t(1,3,iglob) = Qt_t(1,3,iglob)+&
     d1l*(fac1*(temp(3,1,1)*pmlxxl+temp(3,2,1)*pmlxyl+temp(3,3,1)*pmlxzl+&
                temp(5,1,1)*pmlyxl+temp(5,2,1)*pmlyyl+temp(5,3,1)*pmlyzl+&
                temp(6,1,1)*pmlzxl+temp(6,2,1)*pmlzyl+temp(6,3,1)*pmlzzl) +&
          fac2*(temp(3,1,2)*pmlxxl+temp(3,2,2)*pmlxyl+temp(3,3,2)*pmlxzl+&
                temp(5,1,2)*pmlyxl+temp(5,2,2)*pmlyyl+temp(5,3,2)*pmlyzl+&
                temp(6,1,2)*pmlzxl+temp(6,2,2)*pmlzyl+temp(6,3,2)*pmlzzl) +&
          fac3*(temp(3,1,3)*pmlxxl+temp(3,2,3)*pmlxyl+temp(3,3,3)*pmlxzl+&
                temp(5,1,3)*pmlyxl+temp(5,2,3)*pmlyyl+temp(5,3,3)*pmlyzl+&
                temp(6,1,3)*pmlzxl+temp(6,2,3)*pmlzyl+temp(6,3,3)*pmlzzl))

  pmlxxl = r2xl*ri2xl
  pmlxyl = r2xl*ri2yl
  pmlxzl = r2xl*ri2zl
  pmlyxl = r2yl*ri2xl
  pmlyyl = r2yl*ri2yl
  pmlyzl = r2yl*ri2zl
  pmlzxl = r2zl*ri2xl
  pmlzyl = r2zl*ri2yl
  pmlzzl = r2zl*ri2zl
  Qt_t(2,1,iglob) = Qt_t(2,1,iglob)+&
     d2l*(fac1*(temp(1,1,1)*pmlxxl+temp(1,2,1)*pmlxyl+temp(1,3,1)*pmlxzl+&
                temp(2,1,1)*pmlyxl+temp(2,2,1)*pmlyyl+temp(2,3,1)*pmlyzl+&
                temp(3,1,1)*pmlzxl+temp(3,2,1)*pmlzyl+temp(3,3,1)*pmlzzl) +&
          fac2*(temp(1,1,2)*pmlxxl+temp(1,2,2)*pmlxyl+temp(1,3,2)*pmlxzl+&
                temp(2,1,2)*pmlyxl+temp(2,2,2)*pmlyyl+temp(2,3,2)*pmlyzl+&
                temp(3,1,2)*pmlzxl+temp(3,2,2)*pmlzyl+temp(3,3,2)*pmlzzl) +&
          fac3*(temp(1,1,3)*pmlxxl+temp(1,2,3)*pmlxyl+temp(1,3,3)*pmlxzl+&
                temp(2,1,3)*pmlyxl+temp(2,2,3)*pmlyyl+temp(2,3,3)*pmlyzl+&
                temp(3,1,3)*pmlzxl+temp(3,2,3)*pmlzyl+temp(3,3,3)*pmlzzl))
  Qt_t(2,2,iglob) = Qt_t(2,2,iglob)+&
     d2l*(fac1*(temp(2,1,1)*pmlxxl+temp(2,2,1)*pmlxyl+temp(2,3,1)*pmlxzl+&
                temp(4,1,1)*pmlyxl+temp(4,2,1)*pmlyyl+temp(4,3,1)*pmlyzl+&
                temp(5,1,1)*pmlzxl+temp(5,2,1)*pmlzyl+temp(5,3,1)*pmlzzl) +&
          fac2*(temp(2,1,2)*pmlxxl+temp(2,2,2)*pmlxyl+temp(2,3,2)*pmlxzl+&
                temp(4,1,2)*pmlyxl+temp(4,2,2)*pmlyyl+temp(4,3,2)*pmlyzl+&
                temp(5,1,2)*pmlzxl+temp(5,2,2)*pmlzyl+temp(5,3,2)*pmlzzl) +&
          fac3*(temp(2,1,3)*pmlxxl+temp(2,2,3)*pmlxyl+temp(2,3,3)*pmlxzl+&
                temp(4,1,3)*pmlyxl+temp(4,2,3)*pmlyyl+temp(4,3,3)*pmlyzl+&
                temp(5,1,3)*pmlzxl+temp(5,2,3)*pmlzyl+temp(5,3,3)*pmlzzl))
  Qt_t(2,3,iglob) = Qt_t(2,3,iglob)+&
     d2l*(fac1*(temp(3,1,1)*pmlxxl+temp(3,2,1)*pmlxyl+temp(3,3,1)*pmlxzl+&
                temp(5,1,1)*pmlyxl+temp(5,2,1)*pmlyyl+temp(5,3,1)*pmlyzl+&
                temp(6,1,1)*pmlzxl+temp(6,2,1)*pmlzyl+temp(6,3,1)*pmlzzl) +&
          fac2*(temp(3,1,2)*pmlxxl+temp(3,2,2)*pmlxyl+temp(3,3,2)*pmlxzl+&
                temp(5,1,2)*pmlyxl+temp(5,2,2)*pmlyyl+temp(5,3,2)*pmlyzl+&
                temp(6,1,2)*pmlzxl+temp(6,2,2)*pmlzyl+temp(6,3,2)*pmlzzl) +&
          fac3*(temp(3,1,3)*pmlxxl+temp(3,2,3)*pmlxyl+temp(3,3,3)*pmlxzl+&
                temp(5,1,3)*pmlyxl+temp(5,2,3)*pmlyyl+temp(5,3,3)*pmlyzl+&
                temp(6,1,3)*pmlzxl+temp(6,2,3)*pmlzyl+temp(6,3,3)*pmlzzl))

  pmlxxl = r3xl*ri3xl
  pmlxyl = r3xl*ri3yl
  pmlxzl = r3xl*ri3zl
  pmlyxl = r3yl*ri3xl
  pmlyyl = r3yl*ri3yl
  pmlyzl = r3yl*ri3zl
  pmlzxl = r3zl*ri3xl
  pmlzyl = r3zl*ri3yl
  pmlzzl = r3zl*ri3zl
  Qt_t(3,1,iglob) = Qt_t(3,1,iglob)+&
     d3l*(fac1*(temp(1,1,1)*pmlxxl+temp(1,2,1)*pmlxyl+temp(1,3,1)*pmlxzl+&
                temp(2,1,1)*pmlyxl+temp(2,2,1)*pmlyyl+temp(2,3,1)*pmlyzl+&
                temp(3,1,1)*pmlzxl+temp(3,2,1)*pmlzyl+temp(3,3,1)*pmlzzl) +&
          fac2*(temp(1,1,2)*pmlxxl+temp(1,2,2)*pmlxyl+temp(1,3,2)*pmlxzl+&
                temp(2,1,2)*pmlyxl+temp(2,2,2)*pmlyyl+temp(2,3,2)*pmlyzl+&
                temp(3,1,2)*pmlzxl+temp(3,2,2)*pmlzyl+temp(3,3,2)*pmlzzl) +&
          fac3*(temp(1,1,3)*pmlxxl+temp(1,2,3)*pmlxyl+temp(1,3,3)*pmlxzl+&
                temp(2,1,3)*pmlyxl+temp(2,2,3)*pmlyyl+temp(2,3,3)*pmlyzl+&
                temp(3,1,3)*pmlzxl+temp(3,2,3)*pmlzyl+temp(3,3,3)*pmlzzl))
  Qt_t(3,2,iglob) = Qt_t(3,2,iglob)+&
     d3l*(fac1*(temp(2,1,1)*pmlxxl+temp(2,2,1)*pmlxyl+temp(2,3,1)*pmlxzl+&
                temp(4,1,1)*pmlyxl+temp(4,2,1)*pmlyyl+temp(4,3,1)*pmlyzl+&
                temp(5,1,1)*pmlzxl+temp(5,2,1)*pmlzyl+temp(5,3,1)*pmlzzl) +&
          fac2*(temp(2,1,2)*pmlxxl+temp(2,2,2)*pmlxyl+temp(2,3,2)*pmlxzl+&
                temp(4,1,2)*pmlyxl+temp(4,2,2)*pmlyyl+temp(4,3,2)*pmlyzl+&
                temp(5,1,2)*pmlzxl+temp(5,2,2)*pmlzyl+temp(5,3,2)*pmlzzl) +&
          fac3*(temp(2,1,3)*pmlxxl+temp(2,2,3)*pmlxyl+temp(2,3,3)*pmlxzl+&
                temp(4,1,3)*pmlyxl+temp(4,2,3)*pmlyyl+temp(4,3,3)*pmlyzl+&
                temp(5,1,3)*pmlzxl+temp(5,2,3)*pmlzyl+temp(5,3,3)*pmlzzl))
  Qt_t(3,3,iglob) = Qt_t(3,3,iglob)+&
     d3l*(fac1*(temp(3,1,1)*pmlxxl+temp(3,2,1)*pmlxyl+temp(3,3,1)*pmlxzl+&
                temp(5,1,1)*pmlyxl+temp(5,2,1)*pmlyyl+temp(5,3,1)*pmlyzl+&
                temp(6,1,1)*pmlzxl+temp(6,2,1)*pmlzyl+temp(6,3,1)*pmlzzl) +&
          fac2*(temp(3,1,2)*pmlxxl+temp(3,2,2)*pmlxyl+temp(3,3,2)*pmlxzl+&
                temp(5,1,2)*pmlyxl+temp(5,2,2)*pmlyyl+temp(5,3,2)*pmlyzl+&
                temp(6,1,2)*pmlzxl+temp(6,2,2)*pmlzyl+temp(6,3,2)*pmlzzl) +&
          fac3*(temp(3,1,3)*pmlxxl+temp(3,2,3)*pmlxyl+temp(3,3,3)*pmlxzl+&
                temp(5,1,3)*pmlyxl+temp(5,2,3)*pmlyyl+temp(5,3,3)*pmlyzl+&
                temp(6,1,3)*pmlzxl+temp(6,2,3)*pmlzyl+temp(6,3,3)*pmlzzl))
end subroutine compute_Qt_t

subroutine compute_accel_adepml(ispec,i,j,k,temp,fac1,fac2,fac3)
  use constants, only: CUSTOM_REAL
  use pml_par, only: r_trans, r_trans_inv, pml_kappa, pml_d
  use pml_par, only: accel_elastic_CPML
  implicit none
  !real(kind=CUSTOM_REAL), dimension(3,3), intent(out) :: Qt_tl
  real(kind=CUSTOM_REAL) :: k1l,k2l,k3l,d1l,d2l,d3l,&
          r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,&
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,&
          pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl
  integer :: ispec,i,j,k
  real(kind=CUSTOM_REAL), intent(in) :: fac1,fac2,fac3
  real(kind=CUSTOM_REAL), dimension(6,3,3), intent(in) :: temp
  !integer :: iglob
  !iglob = ibool_CPML(i,j,k,ispec)
  r1xl = r_trans(1,1,i,j,k,ispec)
  r1yl = r_trans(1,2,i,j,k,ispec)
  r1zl = r_trans(1,3,i,j,k,ispec)
  r2xl = r_trans(2,1,i,j,k,ispec)
  r2yl = r_trans(2,2,i,j,k,ispec)
  r2zl = r_trans(2,3,i,j,k,ispec)
  r3xl = r_trans(3,1,i,j,k,ispec)
  r3yl = r_trans(3,2,i,j,k,ispec)
  r3zl = r_trans(3,3,i,j,k,ispec)

  ri1xl = r_trans_inv(1,1,i,j,k,ispec)
  ri1yl = r_trans_inv(1,2,i,j,k,ispec)
  ri1zl = r_trans_inv(1,3,i,j,k,ispec)
  ri2xl = r_trans_inv(2,1,i,j,k,ispec)
  ri2yl = r_trans_inv(2,2,i,j,k,ispec)
  ri2zl = r_trans_inv(2,3,i,j,k,ispec)
  ri3xl = r_trans_inv(3,1,i,j,k,ispec)
  ri3yl = r_trans_inv(3,2,i,j,k,ispec)
  ri3zl = r_trans_inv(3,3,i,j,k,ispec)

  k1l = pml_kappa(1,i,j,k,ispec)
  k2l = pml_kappa(2,i,j,k,ispec)
  k3l = pml_kappa(3,i,j,k,ispec)
  !b1l = pml_beta(1,i,j,k,ispec)
  !b2l = pml_beta(2,i,j,k,ispec)
  !b3l = pml_beta(3,i,j,k,ispec)
  d1l = pml_d(1,i,j,k,ispec)
  d2l = pml_d(2,i,j,k,ispec)
  d3l = pml_d(3,i,j,k,ispec)

  pmlxxl = r1xl*ri1xl/k1l + r2xl*ri2xl/k2l + r3xl*ri3xl/k3l
  pmlxyl = r1xl*ri1yl/k1l + r2xl*ri2yl/k2l + r3xl*ri3yl/k3l
  pmlxzl = r1xl*ri1zl/k1l + r2xl*ri2zl/k2l + r3xl*ri3zl/k3l
  pmlyxl = r1yl*ri1xl/k1l + r2yl*ri2xl/k2l + r3yl*ri3xl/k3l
  pmlyyl = r1yl*ri1yl/k1l + r2yl*ri2yl/k2l + r3yl*ri3yl/k3l
  pmlyzl = r1yl*ri1zl/k1l + r2yl*ri2zl/k2l + r3yl*ri3zl/k3l
  pmlzxl = r1zl*ri1xl/k1l + r2zl*ri2xl/k2l + r3zl*ri3xl/k3l
  pmlzyl = r1zl*ri1yl/k1l + r2zl*ri2yl/k2l + r3zl*ri3yl/k3l
  pmlzzl = r1zl*ri1zl/k1l + r2zl*ri2zl/k2l + r3zl*ri3zl/k3l

  ! 1-xx 2-xy 3-xz 4-yy 5-yz 6-zz
  accel_elastic_CPML(1,i,j,k) = &
          fac1*(temp(1,1,1)*pmlxxl+temp(1,2,1)*pmlxyl+temp(1,3,1)*pmlxzl+&
                temp(2,1,1)*pmlyxl+temp(2,2,1)*pmlyyl+temp(2,3,1)*pmlyzl+&
                temp(3,1,1)*pmlzxl+temp(3,2,1)*pmlzyl+temp(3,3,1)*pmlzzl) +&
          fac2*(temp(1,1,2)*pmlxxl+temp(1,2,2)*pmlxyl+temp(1,3,2)*pmlxzl+&
                temp(2,1,2)*pmlyxl+temp(2,2,2)*pmlyyl+temp(2,3,2)*pmlyzl+&
                temp(3,1,2)*pmlzxl+temp(3,2,2)*pmlzyl+temp(3,3,2)*pmlzzl) +&
          fac3*(temp(1,1,3)*pmlxxl+temp(1,2,3)*pmlxyl+temp(1,3,3)*pmlxzl+&
                temp(2,1,3)*pmlyxl+temp(2,2,3)*pmlyyl+temp(2,3,3)*pmlyzl+&
                temp(3,1,3)*pmlzxl+temp(3,2,3)*pmlzyl+temp(3,3,3)*pmlzzl)
  accel_elastic_CPML(2,i,j,k) = &
          fac1*(temp(2,1,1)*pmlxxl+temp(2,2,1)*pmlxyl+temp(2,3,1)*pmlxzl+&
                temp(4,1,1)*pmlyxl+temp(4,2,1)*pmlyyl+temp(4,3,1)*pmlyzl+&
                temp(5,1,1)*pmlzxl+temp(5,2,1)*pmlzyl+temp(5,3,1)*pmlzzl) +&
          fac2*(temp(2,1,2)*pmlxxl+temp(2,2,2)*pmlxyl+temp(2,3,2)*pmlxzl+&
                temp(4,1,2)*pmlyxl+temp(4,2,2)*pmlyyl+temp(4,3,2)*pmlyzl+&
                temp(5,1,2)*pmlzxl+temp(5,2,2)*pmlzyl+temp(5,3,2)*pmlzzl) +&
          fac3*(temp(2,1,3)*pmlxxl+temp(2,2,3)*pmlxyl+temp(2,3,3)*pmlxzl+&
                temp(4,1,3)*pmlyxl+temp(4,2,3)*pmlyyl+temp(4,3,3)*pmlyzl+&
                temp(5,1,3)*pmlzxl+temp(5,2,3)*pmlzyl+temp(5,3,3)*pmlzzl)
  accel_elastic_CPML(3,i,j,k) = &
          fac1*(temp(3,1,1)*pmlxxl+temp(3,2,1)*pmlxyl+temp(3,3,1)*pmlxzl+&
                temp(5,1,1)*pmlyxl+temp(5,2,1)*pmlyyl+temp(5,3,1)*pmlyzl+&
                temp(6,1,1)*pmlzxl+temp(6,2,1)*pmlzyl+temp(6,3,1)*pmlzzl) +&
          fac2*(temp(3,1,2)*pmlxxl+temp(3,2,2)*pmlxyl+temp(3,3,2)*pmlxzl+&
                temp(5,1,2)*pmlyxl+temp(5,2,2)*pmlyyl+temp(5,3,2)*pmlyzl+&
                temp(6,1,2)*pmlzxl+temp(6,2,2)*pmlzyl+temp(6,3,2)*pmlzzl) +&
          fac3*(temp(3,1,3)*pmlxxl+temp(3,2,3)*pmlxyl+temp(3,3,3)*pmlxzl+&
                temp(5,1,3)*pmlyxl+temp(5,2,3)*pmlyyl+temp(5,3,3)*pmlyzl+&
                temp(6,1,3)*pmlzxl+temp(6,2,3)*pmlzyl+temp(6,3,3)*pmlzzl)
end subroutine compute_accel_adepml

subroutine compute_add_pml_physical_element(ispec,sigma_elem)
  use constants, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NDIM, NGLLSQUARE
  use pml_par, only: pml_spec_physical
  use pml_par, only: r_trans, r_trans_inv, pml_d, &
           Qt_t, pml_physical_ispec, ibool_CPML, &
           pml_physical_ijk, pml_physical_normal, pml_physical_jacobian2Dw
  !use specfem_par, only: ibool
  !!local variables
  implicit none
  integer :: iface, ispec, i, j, k, igll, iglob, iside
  real(kind=CUSTOM_REAL) :: nx, ny, nz, &
          sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yz,sigma_zz, &
          d1l,d2l,d3l,r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,&
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,&
          weight, tx, ty, tz
  real(kind=CUSTOM_REAL), dimension(6,NGLLX,NGLLY,NGLLZ) :: sigma_elem
  do iside = 1, 2*NDIM
    iface = pml_spec_physical(iside,ispec)
    if (iface .eq. 0) cycle
    if (.not. (ispec .eq. pml_physical_ispec(iface))) &
      print *, 'wrong interface element index'
    do igll = 1, NGLLSQUARE
      i = pml_physical_ijk(1,igll,iface)
      j = pml_physical_ijk(2,igll,iface)
      k = pml_physical_ijk(3,igll,iface)
      iglob = ibool_CPML(i,j,k,ispec)
      nx = pml_physical_normal(1,igll,iface)
      ny = pml_physical_normal(2,igll,iface)
      nz = pml_physical_normal(3,igll,iface)
      
      r1xl = r_trans(1,1,i,j,k,ispec)
      r1yl = r_trans(1,2,i,j,k,ispec)
      r1zl = r_trans(1,3,i,j,k,ispec)
      r2xl = r_trans(2,1,i,j,k,ispec)
      r2yl = r_trans(2,2,i,j,k,ispec)
      r2zl = r_trans(2,3,i,j,k,ispec)
      r3xl = r_trans(3,1,i,j,k,ispec)
      r3yl = r_trans(3,2,i,j,k,ispec)
      r3zl = r_trans(3,3,i,j,k,ispec)

      ri1xl = r_trans_inv(1,1,i,j,k,ispec)
      ri1yl = r_trans_inv(1,2,i,j,k,ispec)
      ri1zl = r_trans_inv(1,3,i,j,k,ispec)
      ri2xl = r_trans_inv(2,1,i,j,k,ispec)
      ri2yl = r_trans_inv(2,2,i,j,k,ispec)
      ri2zl = r_trans_inv(2,3,i,j,k,ispec)
      ri3xl = r_trans_inv(3,1,i,j,k,ispec)
      ri3yl = r_trans_inv(3,2,i,j,k,ispec)
      ri3zl = r_trans_inv(3,3,i,j,k,ispec)

      d1l = pml_d(1,i,j,k,ispec)
      d2l = pml_d(2,i,j,k,ispec)
      d3l = pml_d(3,i,j,k,ispec)

      sigma_xx = sigma_elem(1,i,j,k)
      sigma_xy = sigma_elem(2,i,j,k)
      sigma_xz = sigma_elem(3,i,j,k)
      sigma_yy = sigma_elem(4,i,j,k)
      sigma_yz = sigma_elem(5,i,j,k)
      sigma_zz = sigma_elem(6,i,j,k)

      weight = pml_physical_jacobian2Dw(igll,iface)
      tx = d1l*(ri1xl*nx+ri1yl*ny+ri1zl*nz)*&
                          (r1xl*sigma_xx+r1yl*sigma_xy+r1zl*sigma_xz)
      ty = d1l*(ri1xl*nx+ri1yl*ny+ri1zl*nz)*&
                          (r1xl*sigma_xy+r1yl*sigma_yy+r1zl*sigma_yz)
      tz = d1l*(ri1xl*nx+ri1yl*ny+ri1zl*nz)*&
                          (r1xl*sigma_xz+r1yl*sigma_yz+r1zl*sigma_zz)
      Qt_t(1,1,iglob) = Qt_t(1,1,iglob) - tx * weight
      Qt_t(1,2,iglob) = Qt_t(1,2,iglob) - ty * weight
      Qt_t(1,3,iglob) = Qt_t(1,3,iglob) - tz * weight
      tx = d2l*(ri2xl*nx+ri2yl*ny+ri2zl*nz)*&
                          (r2xl*sigma_xx+r2yl*sigma_xy+r2zl*sigma_xz)
      ty = d2l*(ri2xl*nx+ri2yl*ny+ri2zl*nz)*&
                          (r2xl*sigma_xy+r2yl*sigma_yy+r2zl*sigma_yz)
      tz = d2l*(ri2xl*nx+ri2yl*ny+ri2zl*nz)*&
                          (r2xl*sigma_xz+r2yl*sigma_yz+r2zl*sigma_zz)
      Qt_t(2,1,iglob) = Qt_t(2,1,iglob) - tx * weight
      Qt_t(2,2,iglob) = Qt_t(2,2,iglob) - ty * weight
      Qt_t(2,3,iglob) = Qt_t(2,3,iglob) - tz * weight
      tx = d3l*(ri3xl*nx+ri3yl*ny+ri3zl*nz)*&
                          (r3xl*sigma_xx+r3yl*sigma_xy+r3zl*sigma_xz)
      ty = d3l*(ri3xl*nx+ri3yl*ny+ri3zl*nz)*&
                          (r3xl*sigma_xy+r3yl*sigma_yy+r3zl*sigma_yz)
      tz = d3l*(ri3xl*nx+ri3yl*ny+ri3zl*nz)*&
                          (r3xl*sigma_xz+r3yl*sigma_yz+r3zl*sigma_zz)
      Qt_t(3,1,iglob) = Qt_t(3,1,iglob) - tx * weight
      Qt_t(3,2,iglob) = Qt_t(3,2,iglob) - ty * weight
      Qt_t(3,3,iglob) = Qt_t(3,3,iglob) - tz * weight
    enddo  ! igll = 1, NGLLSQUARE
  enddo ! iside = 1, NDIM
end subroutine compute_add_pml_physical_element

subroutine update_Qt_conv()
  !use specfem_par, only: deltat, deltatover2
  use pml_par, only: Qt_t,Qt,nglob_CPML,rvolume,&
                     coeff_glob_exp1,coeff_glob_exp2, &
                     nglob_pml_in, pml_in_iglob
  implicit none
  integer :: i
  do i = 1,nglob_CPML
    Qt_t(:,:,i) = Qt_t(:,:,i)*rvolume(i)
  enddo
  do i = 1, 3
    Qt(:,i,:)=Qt(:,i,:)*coeff_glob_exp1(:,:)+coeff_glob_exp2(:,:)*Qt_t(:,i,:)
  enddo
  !do i = 1, nglob_pml_in
  !  Qt(:,:,pml_in_iglob(i)) = 0.0
  !enddo
  Qt_t(:,:,:) = 0.0
  !Qt_t_ext(:,:,:) = 0.0
end subroutine update_Qt_conv

subroutine update_Qu_conv()
  !use specfem_par, only: deltat, deltatover2
  use pml_par, only: Qu_t,Qu,coeff_exp1,coeff_exp2
  implicit none
  integer :: i
  do i = 1, 3
    Qu(:,i,:,:,:,:)=Qu(:,i,:,:,:,:)*coeff_exp1(:,:,:,:,:)+&
                  coeff_exp2(:,:,:,:,:)*Qu_t(:,i,:,:,:,:)
  enddo
  Qu_t(:,:,:,:,:,:) = 0.0
end subroutine update_Qu_conv


subroutine include_adepml_accel_aux()
  use pml_par, only: nglob_CPML,CPML_to_glob,Qt,rvolume
  use specfem_par_elastic, only: rmassx, rmassy, rmassz, accel
  integer :: iglob_CPML, iglob
  do iglob_CPML = 1, nglob_CPML
    iglob = CPML_to_glob(iglob_CPML)
    accel(1,iglob)=accel(1,iglob)+(Qt(1,1,iglob_CPML)+Qt(2,1,iglob_CPML)+ &
                     Qt(3,1,iglob_CPML))*rmassx(iglob)/rvolume(iglob_CPML)
    accel(2,iglob)=accel(2,iglob)+(Qt(1,2,iglob_CPML)+Qt(2,2,iglob_CPML)+ &
                     Qt(3,2,iglob_CPML))*rmassy(iglob)/rvolume(iglob_CPML)
    accel(3,iglob)=accel(3,iglob)+(Qt(1,3,iglob_CPML)+Qt(2,3,iglob_CPML)+ &
                     Qt(3,3,iglob_CPML))*rmassz(iglob)/rvolume(iglob_CPML)
    
  enddo
end subroutine include_adepml_accel_aux

subroutine prepare_adepml_update_coeff()
  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
  use specfem_par, only: deltat
  use pml_par, only: coeff_exp1,coeff_exp2,coeff_glob_exp1,&
    coeff_glob_exp2,pml_beta,nspec_CPML,ibool_CPML
  implicit none
  integer :: i,j,k,iglob,ispec
  if (nspec_CPML == 0) return
  coeff_exp1(:,:,:,:,:)=exp(-pml_beta(:,:,:,:,:)*deltat)
  !coeff_exp3(:,:,:,:,:)=sqrt(coeff_exp1(:,:,:,:,:))
  coeff_exp2(:,:,:,:,:)=(1._CUSTOM_REAL-coeff_exp1(:,:,:,:,:))/pml_beta(:,:,:,:,:)
  do ispec = 1, nspec_CPML
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_CPML(i,j,k,ispec)
          coeff_glob_exp1(:,iglob)=coeff_exp1(:,i,j,k,ispec)
          !coeff_glob_exp2(:,iglob)=coeff_exp2(:,i,j,k,ispec)
          coeff_glob_exp2(:,iglob)=coeff_exp2(:,i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo
  if (maxval(coeff_exp1) > 1.0) print *, 'bad value coeff_exp1'
  !if (maxval(coeff_exp2) > 1.0) print *, 'bad value coeff_exp3'
  if (maxval(coeff_glob_exp1) > 1.0) print *, 'bad value coeff_glob_exp1'
  !if (maxval(coeff_glob_exp3) > 1.0) print *, 'bad value coeff_glob_exp3'
end subroutine prepare_adepml_update_coeff
