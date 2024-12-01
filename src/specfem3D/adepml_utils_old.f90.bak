subroutine compute_sigma_adepml(ispec,i,j,k,lambdalplus2mul,lambdal,mul, &
                                sigma_xx,sigma_xy,sigma_xz, &
                                sigma_yy,sigma_yz,sigma_zz, &
                                duxdxl, duydxl, duzdxl, &
                                duxdyl, duydyl, duzdyl, &
                                duxdzl, duydzl, duzdzl)
  use constants, only: CUSTOM_REAL
  use pml_par, only: Qu,Qu_t, r_trans, r_trans_inv, pml_kappa, pml_d
  implicit none
  real(kind=CUSTOM_REAL) :: k1l,k2l,k3l,d1l,d2l,d3l,&
          r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,&
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,&
          pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl,&
          Qu1xl,Qu1yl,Qu1zl,Qu2xl,Qu2yl,Qu2zl,Qu3xl,Qu3yl,Qu3zl
  real(kind=CUSTOM_REAL), intent(in) :: lambdalplus2mul,lambdal,mul, &
                                duxdxl, duydxl, duzdxl, &
                                duxdyl, duydyl, duzdyl, &
                                duxdzl, duydzl, duzdzl
  real(kind=CUSTOM_REAL), intent(out) :: sigma_xx,sigma_xy,sigma_xz, &
                                         sigma_yy,sigma_yz,sigma_zz
  integer, intent(in) :: ispec,i,j,k
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

  Qu1xl = Qu(1,1,i,j,k,ispec)
  Qu1yl = Qu(1,2,i,j,k,ispec)
  Qu1zl = Qu(1,3,i,j,k,ispec)
  Qu2xl = Qu(2,1,i,j,k,ispec)
  Qu2yl = Qu(2,2,i,j,k,ispec)
  Qu2zl = Qu(2,3,i,j,k,ispec)
  Qu3xl = Qu(3,1,i,j,k,ispec)
  Qu3yl = Qu(3,2,i,j,k,ispec)
  Qu3zl = Qu(3,3,i,j,k,ispec)
  
  sigma_xx = lambdalplus2mul*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
           r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
             lambdal*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
           r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
             lambdal*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
           r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)
  sigma_xy = mul*(pmlxxl*duydxl+pmlxyl*duydyl+pmlxzl*duydzl+ &
           r1xl*Qu1yl+r2xl*Qu2yl+r3xl*Qu3yl) + &
             mul*(pmlyxl*duxdxl+pmlyyl*duxdyl+pmlyzl*duxdzl+ &
           r1yl*Qu1xl+r2yl*Qu2xl+r3yl*Qu3xl)
  sigma_xz = mul*(pmlxxl*duzdxl+pmlxyl*duzdyl+pmlxzl*duzdzl+ &
           r1xl*Qu1zl+r2xl*Qu2zl+r3xl*Qu3zl) + &
             mul*(pmlzxl*duxdxl+pmlzyl*duxdyl+pmlzzl*duxdzl+ &
           r1zl*Qu1xl+r2zl*Qu2xl+r3zl*Qu3xl)
  sigma_yy = lambdal*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
           r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
             lambdalplus2mul*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
           r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
             lambdal*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
           r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)
  sigma_yz = mul*(pmlyxl*duzdxl+pmlyyl*duzdyl+pmlyzl*duzdzl+ &
           r1yl*Qu1zl+r2yl*Qu2zl+r3yl*Qu3zl) + &
             mul*(pmlzxl*duydxl+pmlzyl*duydyl+pmlzzl*duydzl+ &
           r1zl*Qu1yl+r2zl*Qu2yl+r3zl*Qu3yl)
  sigma_zz = lambdal*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
           r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
             lambdal*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
           r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
             lambdalplus2mul*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
           r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)
  
  Qu_t(1,1,i,j,k,ispec)=-d1l*(ri1xl*duxdxl+ri1yl*duxdyl+ri1zl*duxdzl)
  Qu_t(1,2,i,j,k,ispec)=-d1l*(ri1xl*duydxl+ri1yl*duydyl+ri1zl*duydzl)
  Qu_t(1,3,i,j,k,ispec)=-d1l*(ri1xl*duzdxl+ri1yl*duzdyl+ri1zl*duzdzl)
  Qu_t(2,1,i,j,k,ispec)=-d2l*(ri2xl*duxdxl+ri2yl*duxdyl+ri2zl*duxdzl)
  Qu_t(2,2,i,j,k,ispec)=-d2l*(ri2xl*duydxl+ri2yl*duydyl+ri2zl*duydzl)
  Qu_t(2,3,i,j,k,ispec)=-d2l*(ri2xl*duzdxl+ri2yl*duzdyl+ri2zl*duzdzl)
  Qu_t(3,1,i,j,k,ispec)=-d3l*(ri3xl*duxdxl+ri3yl*duxdyl+ri3zl*duxdzl)
  Qu_t(3,2,i,j,k,ispec)=-d3l*(ri3xl*duydxl+ri3yl*duydyl+ri3zl*duydzl)
  Qu_t(3,3,i,j,k,ispec)=-d3l*(ri3xl*duzdxl+ri3yl*duzdyl+ri3zl*duzdzl)
end subroutine compute_sigma_adepml


subroutine compute_accel_adepml(ispec,i,j,k,temp,jacobianl,&
                fac1,fac2,fac3,wgll3_xyzl)
  use constants, only: CUSTOM_REAL
  use pml_par, only: Qt, Qt_t, r_trans, r_trans_inv, pml_kappa, &
          pml_d, ibool_CPML
  use pml_par, only: accel_elastic_CPML
  implicit none
  !real(kind=CUSTOM_REAL), dimension(3,3), intent(out) :: Qt_tl
  real(kind=CUSTOM_REAL) :: k1l,k2l,k3l,d1l,d2l,d3l,&
          r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,&
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,&
          pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl,&
          Qt1xl,Qt1yl,Qt1zl,Qt2xl,Qt2yl,Qt2zl,Qt3xl,Qt3yl,Qt3zl
  integer :: ispec,i,j,k
  real(kind=CUSTOM_REAL), intent(in) :: jacobianl,fac1,fac2,fac3,wgll3_xyzl
  real(kind=CUSTOM_REAL), dimension(6,3,3), intent(in) :: temp
  integer :: iglob
  iglob = ibool_CPML(i,j,k,ispec)
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

  Qt1xl = Qt(1,1,iglob)
  Qt1yl = Qt(1,2,iglob)
  Qt1zl = Qt(1,3,iglob)
  Qt2xl = Qt(2,1,iglob)
  Qt2yl = Qt(2,2,iglob)
  Qt2zl = Qt(2,3,iglob)
  Qt3xl = Qt(3,1,iglob)
  Qt3yl = Qt(3,2,iglob)
  Qt3zl = Qt(3,3,iglob)
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
                temp(3,1,3)*pmlzxl+temp(3,2,3)*pmlzyl+temp(3,3,3)*pmlzzl) -&
          wgll3_xyzl*jacobianl*(Qt1xl+Qt2xl+Qt3xl)
  accel_elastic_CPML(2,i,j,k) = &
          fac1*(temp(2,1,1)*pmlxxl+temp(2,2,1)*pmlxyl+temp(2,3,1)*pmlxzl+&
                temp(4,1,1)*pmlyxl+temp(4,2,1)*pmlyyl+temp(4,3,1)*pmlyzl+&
                temp(5,1,1)*pmlzxl+temp(5,2,1)*pmlzyl+temp(5,3,1)*pmlzzl) +&
          fac2*(temp(2,1,2)*pmlxxl+temp(2,2,2)*pmlxyl+temp(2,3,2)*pmlxzl+&
                temp(4,1,2)*pmlyxl+temp(4,2,2)*pmlyyl+temp(4,3,2)*pmlyzl+&
                temp(5,1,2)*pmlzxl+temp(5,2,2)*pmlzyl+temp(5,3,2)*pmlzzl) +&
          fac3*(temp(2,1,3)*pmlxxl+temp(2,2,3)*pmlxyl+temp(2,3,3)*pmlxzl+&
                temp(4,1,3)*pmlyxl+temp(4,2,3)*pmlyyl+temp(4,3,3)*pmlyzl+&
                temp(5,1,3)*pmlzxl+temp(5,2,3)*pmlzyl+temp(5,3,3)*pmlzzl) -&
          wgll3_xyzl*jacobianl*(Qt1yl+Qt2yl+Qt3yl)
  accel_elastic_CPML(3,i,j,k) = &
          fac1*(temp(3,1,1)*pmlxxl+temp(3,2,1)*pmlxyl+temp(3,3,1)*pmlxzl+&
                temp(5,1,1)*pmlyxl+temp(5,2,1)*pmlyyl+temp(5,3,1)*pmlyzl+&
                temp(6,1,1)*pmlzxl+temp(6,2,1)*pmlzyl+temp(6,3,1)*pmlzzl) +&
          fac2*(temp(3,1,2)*pmlxxl+temp(3,2,2)*pmlxyl+temp(3,3,2)*pmlxzl+&
                temp(5,1,2)*pmlyxl+temp(5,2,2)*pmlyyl+temp(5,3,2)*pmlyzl+&
                temp(6,1,2)*pmlzxl+temp(6,2,2)*pmlzyl+temp(6,3,2)*pmlzzl) +&
          fac3*(temp(3,1,3)*pmlxxl+temp(3,2,3)*pmlxyl+temp(3,3,3)*pmlxzl+&
                temp(5,1,3)*pmlyxl+temp(5,2,3)*pmlyyl+temp(5,3,3)*pmlyzl+&
                temp(6,1,3)*pmlzxl+temp(6,2,3)*pmlzyl+temp(6,3,3)*pmlzzl) -&
          wgll3_xyzl*jacobianl*(Qt1zl+Qt2zl+Qt3zl)
  pmlxxl = r1xl*ri1xl*d1l
  pmlxyl = r1xl*ri1yl*d1l
  pmlxzl = r1xl*ri1zl*d1l
  pmlyxl = r1yl*ri1xl*d1l
  pmlyyl = r1yl*ri1yl*d1l
  pmlyzl = r1yl*ri1zl*d1l
  pmlzxl = r1zl*ri1xl*d1l
  pmlzyl = r1zl*ri1yl*d1l
  pmlzzl = r1zl*ri1zl*d1l
  Qt_t(1,1,iglob) = Qt_t(1,1,iglob)+&
          fac1*(temp(1,1,1)*pmlxxl+temp(1,2,1)*pmlxyl+temp(1,3,1)*pmlxzl+&
                temp(2,1,1)*pmlyxl+temp(2,2,1)*pmlyyl+temp(2,3,1)*pmlyzl+&
                temp(3,1,1)*pmlzxl+temp(3,2,1)*pmlzyl+temp(3,3,1)*pmlzzl) +&
          fac2*(temp(1,1,2)*pmlxxl+temp(1,2,2)*pmlxyl+temp(1,3,2)*pmlxzl+&
                temp(2,1,2)*pmlyxl+temp(2,2,2)*pmlyyl+temp(2,3,2)*pmlyzl+&
                temp(3,1,2)*pmlzxl+temp(3,2,2)*pmlzyl+temp(3,3,2)*pmlzzl) +&
          fac3*(temp(1,1,3)*pmlxxl+temp(1,2,3)*pmlxyl+temp(1,3,3)*pmlxzl+&
                temp(2,1,3)*pmlyxl+temp(2,2,3)*pmlyyl+temp(2,3,3)*pmlyzl+&
                temp(3,1,3)*pmlzxl+temp(3,2,3)*pmlzyl+temp(3,3,3)*pmlzzl) 
  Qt_t(1,2,iglob) = Qt_t(1,2,iglob)+&
          fac1*(temp(2,1,1)*pmlxxl+temp(2,2,1)*pmlxyl+temp(2,3,1)*pmlxzl+&
                temp(4,1,1)*pmlyxl+temp(4,2,1)*pmlyyl+temp(4,3,1)*pmlyzl+&
                temp(5,1,1)*pmlzxl+temp(5,2,1)*pmlzyl+temp(5,3,1)*pmlzzl) +&
          fac2*(temp(2,1,2)*pmlxxl+temp(2,2,2)*pmlxyl+temp(2,3,2)*pmlxzl+&
                temp(4,1,2)*pmlyxl+temp(4,2,2)*pmlyyl+temp(4,3,2)*pmlyzl+&
                temp(5,1,2)*pmlzxl+temp(5,2,2)*pmlzyl+temp(5,3,2)*pmlzzl) +&
          fac3*(temp(2,1,3)*pmlxxl+temp(2,2,3)*pmlxyl+temp(2,3,3)*pmlxzl+&
                temp(4,1,3)*pmlyxl+temp(4,2,3)*pmlyyl+temp(4,3,3)*pmlyzl+&
                temp(5,1,3)*pmlzxl+temp(5,2,3)*pmlzyl+temp(5,3,3)*pmlzzl)
  Qt_t(1,3,iglob) = Qt_t(1,3,iglob)+&
          fac1*(temp(3,1,1)*pmlxxl+temp(3,2,1)*pmlxyl+temp(3,3,1)*pmlxzl+&
                temp(5,1,1)*pmlyxl+temp(5,2,1)*pmlyyl+temp(5,3,1)*pmlyzl+&
                temp(6,1,1)*pmlzxl+temp(6,2,1)*pmlzyl+temp(6,3,1)*pmlzzl) +&
          fac2*(temp(3,1,2)*pmlxxl+temp(3,2,2)*pmlxyl+temp(3,3,2)*pmlxzl+&
                temp(5,1,2)*pmlyxl+temp(5,2,2)*pmlyyl+temp(5,3,2)*pmlyzl+&
                temp(6,1,2)*pmlzxl+temp(6,2,2)*pmlzyl+temp(6,3,2)*pmlzzl) +&
          fac3*(temp(3,1,3)*pmlxxl+temp(3,2,3)*pmlxyl+temp(3,3,3)*pmlxzl+&
                temp(5,1,3)*pmlyxl+temp(5,2,3)*pmlyyl+temp(5,3,3)*pmlyzl+&
                temp(6,1,3)*pmlzxl+temp(6,2,3)*pmlzyl+temp(6,3,3)*pmlzzl)

  pmlxxl = r2xl*ri2xl*d2l
  pmlxyl = r2xl*ri2yl*d2l
  pmlxzl = r2xl*ri2zl*d2l
  pmlyxl = r2yl*ri2xl*d2l
  pmlyyl = r2yl*ri2yl*d2l
  pmlyzl = r2yl*ri2zl*d2l
  pmlzxl = r2zl*ri2xl*d2l
  pmlzyl = r2zl*ri2yl*d2l
  pmlzzl = r2zl*ri2zl*d2l
  Qt_t(2,1,iglob) = Qt_t(2,1,iglob)+&
          fac1*(temp(1,1,1)*pmlxxl+temp(1,2,1)*pmlxyl+temp(1,3,1)*pmlxzl+&
                temp(2,1,1)*pmlyxl+temp(2,2,1)*pmlyyl+temp(2,3,1)*pmlyzl+&
                temp(3,1,1)*pmlzxl+temp(3,2,1)*pmlzyl+temp(3,3,1)*pmlzzl) +&
          fac2*(temp(1,1,2)*pmlxxl+temp(1,2,2)*pmlxyl+temp(1,3,2)*pmlxzl+&
                temp(2,1,2)*pmlyxl+temp(2,2,2)*pmlyyl+temp(2,3,2)*pmlyzl+&
                temp(3,1,2)*pmlzxl+temp(3,2,2)*pmlzyl+temp(3,3,2)*pmlzzl) +&
          fac3*(temp(1,1,3)*pmlxxl+temp(1,2,3)*pmlxyl+temp(1,3,3)*pmlxzl+&
                temp(2,1,3)*pmlyxl+temp(2,2,3)*pmlyyl+temp(2,3,3)*pmlyzl+&
                temp(3,1,3)*pmlzxl+temp(3,2,3)*pmlzyl+temp(3,3,3)*pmlzzl)
  Qt_t(2,2,iglob) = Qt_t(2,2,iglob)+&
          fac1*(temp(2,1,1)*pmlxxl+temp(2,2,1)*pmlxyl+temp(2,3,1)*pmlxzl+&
                temp(4,1,1)*pmlyxl+temp(4,2,1)*pmlyyl+temp(4,3,1)*pmlyzl+&
                temp(5,1,1)*pmlzxl+temp(5,2,1)*pmlzyl+temp(5,3,1)*pmlzzl) +&
          fac2*(temp(2,1,2)*pmlxxl+temp(2,2,2)*pmlxyl+temp(2,3,2)*pmlxzl+&
                temp(4,1,2)*pmlyxl+temp(4,2,2)*pmlyyl+temp(4,3,2)*pmlyzl+&
                temp(5,1,2)*pmlzxl+temp(5,2,2)*pmlzyl+temp(5,3,2)*pmlzzl) +&
          fac3*(temp(2,1,3)*pmlxxl+temp(2,2,3)*pmlxyl+temp(2,3,3)*pmlxzl+&
                temp(4,1,3)*pmlyxl+temp(4,2,3)*pmlyyl+temp(4,3,3)*pmlyzl+&
                temp(5,1,3)*pmlzxl+temp(5,2,3)*pmlzyl+temp(5,3,3)*pmlzzl)
  Qt_t(2,3,iglob) = Qt_t(2,3,iglob)+&
          fac1*(temp(3,1,1)*pmlxxl+temp(3,2,1)*pmlxyl+temp(3,3,1)*pmlxzl+&
                temp(5,1,1)*pmlyxl+temp(5,2,1)*pmlyyl+temp(5,3,1)*pmlyzl+&
                temp(6,1,1)*pmlzxl+temp(6,2,1)*pmlzyl+temp(6,3,1)*pmlzzl) +&
          fac2*(temp(3,1,2)*pmlxxl+temp(3,2,2)*pmlxyl+temp(3,3,2)*pmlxzl+&
                temp(5,1,2)*pmlyxl+temp(5,2,2)*pmlyyl+temp(5,3,2)*pmlyzl+&
                temp(6,1,2)*pmlzxl+temp(6,2,2)*pmlzyl+temp(6,3,2)*pmlzzl) +&
          fac3*(temp(3,1,3)*pmlxxl+temp(3,2,3)*pmlxyl+temp(3,3,3)*pmlxzl+&
                temp(5,1,3)*pmlyxl+temp(5,2,3)*pmlyyl+temp(5,3,3)*pmlyzl+&
                temp(6,1,3)*pmlzxl+temp(6,2,3)*pmlzyl+temp(6,3,3)*pmlzzl)

  pmlxxl = r3xl*ri3xl*d3l
  pmlxyl = r3xl*ri3yl*d3l
  pmlxzl = r3xl*ri3zl*d3l
  pmlyxl = r3yl*ri3xl*d3l
  pmlyyl = r3yl*ri3yl*d3l
  pmlyzl = r3yl*ri3zl*d3l
  pmlzxl = r3zl*ri3xl*d3l
  pmlzyl = r3zl*ri3yl*d3l
  pmlzzl = r3zl*ri3zl*d3l
  Qt_t(3,1,iglob) = Qt_t(3,1,iglob)+&
          fac1*(temp(1,1,1)*pmlxxl+temp(1,2,1)*pmlxyl+temp(1,3,1)*pmlxzl+&
                temp(2,1,1)*pmlyxl+temp(2,2,1)*pmlyyl+temp(2,3,1)*pmlyzl+&
                temp(3,1,1)*pmlzxl+temp(3,2,1)*pmlzyl+temp(3,3,1)*pmlzzl) +&
          fac2*(temp(1,1,2)*pmlxxl+temp(1,2,2)*pmlxyl+temp(1,3,2)*pmlxzl+&
                temp(2,1,2)*pmlyxl+temp(2,2,2)*pmlyyl+temp(2,3,2)*pmlyzl+&
                temp(3,1,2)*pmlzxl+temp(3,2,2)*pmlzyl+temp(3,3,2)*pmlzzl) +&
          fac3*(temp(1,1,3)*pmlxxl+temp(1,2,3)*pmlxyl+temp(1,3,3)*pmlxzl+&
                temp(2,1,3)*pmlyxl+temp(2,2,3)*pmlyyl+temp(2,3,3)*pmlyzl+&
                temp(3,1,3)*pmlzxl+temp(3,2,3)*pmlzyl+temp(3,3,3)*pmlzzl)
  Qt_t(3,2,iglob) = Qt_t(3,2,iglob)+&
          fac1*(temp(2,1,1)*pmlxxl+temp(2,2,1)*pmlxyl+temp(2,3,1)*pmlxzl+&
                temp(4,1,1)*pmlyxl+temp(4,2,1)*pmlyyl+temp(4,3,1)*pmlyzl+&
                temp(5,1,1)*pmlzxl+temp(5,2,1)*pmlzyl+temp(5,3,1)*pmlzzl) +&
          fac2*(temp(2,1,2)*pmlxxl+temp(2,2,2)*pmlxyl+temp(2,3,2)*pmlxzl+&
                temp(4,1,2)*pmlyxl+temp(4,2,2)*pmlyyl+temp(4,3,2)*pmlyzl+&
                temp(5,1,2)*pmlzxl+temp(5,2,2)*pmlzyl+temp(5,3,2)*pmlzzl) +&
          fac3*(temp(2,1,3)*pmlxxl+temp(2,2,3)*pmlxyl+temp(2,3,3)*pmlxzl+&
                temp(4,1,3)*pmlyxl+temp(4,2,3)*pmlyyl+temp(4,3,3)*pmlyzl+&
                temp(5,1,3)*pmlzxl+temp(5,2,3)*pmlzyl+temp(5,3,3)*pmlzzl)
  Qt_t(3,3,iglob) = Qt_t(3,3,iglob)+&
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

subroutine compute_add_pml_physical()
  use constants, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, &
          NDIM, NGLLSQUARE, FOUR_THIRDS
  use pml_par, only: CPML_to_spec
  use pml_par, only: Qu, r_trans, r_trans_inv, pml_kappa, pml_d, &
           Qt_t, num_pml_physical, pml_physical_ispec, ibool_CPML, &
           pml_physical_ijk, pml_physical_normal, pml_physical_jacobian2Dw
  use specfem_par, only: ibool,xix,xiy,xiz,etax,etay,etaz,&
          gammax,gammay,gammaz,kappastore,mustore
  use specfem_par_elastic, only: displ 
  !!local variables
  implicit none
  integer :: iface, ispec, i, j, k, igll, iglob, ispec_glob
  real(kind=CUSTOM_REAL) :: nx, ny, nz, kappal, mul, lambdalplus2mul, lambdal,&
          duxdxl,duydxl,duzdxl,duxdyl,duydyl,duzdyl,duxdzl,duydzl,duzdzl, &
          sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yz,sigma_zz, &
          k1l,k2l,k3l,d1l,d2l,d3l,&
          r1xl,r1yl,r1zl,r2xl,r2yl,r2zl,r3xl,r3yl,r3zl,&
          ri1xl,ri1yl,ri1zl,ri2xl,ri2yl,ri2zl,ri3xl,ri3yl,ri3zl,&
          pmlxxl,pmlxyl,pmlxzl,pmlyxl,pmlyyl,pmlyzl,pmlzxl,pmlzyl,pmlzzl,&
          Qu1xl,Qu1yl,Qu1zl,Qu2xl,Qu2yl,Qu2zl,Qu3xl,Qu3yl,Qu3zl, &
          weight, tx, ty, tz
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: &
          dx_elem, dy_elem, dz_elem, val_elem
  do iface = 1, num_pml_physical
    ispec = pml_physical_ispec(iface)
    ispec_glob = CPML_to_spec(ispec)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec_glob)
          val_elem(1:NDIM,i,j,k) = displ(1:NDIM,iglob)
        enddo
      enddo
    enddo
    call get_gradient_element(val_elem, dx_elem, dy_elem, dz_elem, &
     xix(:,:,:,ispec_glob),xiy(:,:,:,ispec_glob),xiz(:,:,:,ispec_glob), &
     etax(:,:,:,ispec_glob),etay(:,:,:,ispec_glob),etaz(:,:,:,ispec_glob),&
     gammax(:,:,:,ispec_glob),gammay(:,:,:,ispec_glob),gammaz(:,:,:,ispec_glob))
    do igll = 1, NGLLSQUARE
      i = pml_physical_ijk(1,igll,iface)
      j = pml_physical_ijk(2,igll,iface)
      k = pml_physical_ijk(3,igll,iface)
      iglob = ibool_CPML(i,j,k,ispec)
      nx = pml_physical_normal(1,igll,iface)
      ny = pml_physical_normal(2,igll,iface)
      nz = pml_physical_normal(3,igll,iface)

      duxdxl = dx_elem(1,i,j,k)
      duydxl = dx_elem(2,i,j,k)
      duzdxl = dx_elem(3,i,j,k)
      duxdyl = dy_elem(1,i,j,k)
      duydyl = dy_elem(2,i,j,k)
      duzdyl = dy_elem(3,i,j,k)
      duxdzl = dz_elem(1,i,j,k)
      duydzl = dz_elem(2,i,j,k)
      duzdzl = dz_elem(3,i,j,k)

      kappal = kappastore(i,j,k,ispec_glob)
      mul = mustore(i,j,k,ispec_glob)
      lambdalplus2mul = kappal + FOUR_THIRDS * mul
      lambdal = lambdalplus2mul - 2.0 * mul
      
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

      Qu1xl = Qu(1,1,i,j,k,ispec)
      Qu1yl = Qu(1,2,i,j,k,ispec)
      Qu1zl = Qu(1,3,i,j,k,ispec)
      Qu2xl = Qu(2,1,i,j,k,ispec)
      Qu2yl = Qu(2,2,i,j,k,ispec)
      Qu2zl = Qu(2,3,i,j,k,ispec)
      Qu3xl = Qu(3,1,i,j,k,ispec)
      Qu3yl = Qu(3,2,i,j,k,ispec)
      Qu3zl = Qu(3,3,i,j,k,ispec)

      sigma_xx = lambdalplus2mul*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
               r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
                 lambdal*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
               r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
                 lambdal*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
               r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)
      sigma_xy = mul*(pmlxxl*duydxl+pmlxyl*duydyl+pmlxzl*duydzl+ &
               r1xl*Qu1yl+r2xl*Qu2yl+r3xl*Qu3yl) + &
                 mul*(pmlyxl*duxdxl+pmlyyl*duxdyl+pmlyzl*duxdzl+ &
               r1yl*Qu1xl+r2yl*Qu2xl+r3yl*Qu3xl)
      sigma_xz = mul*(pmlxxl*duzdxl+pmlxyl*duzdyl+pmlxzl*duzdzl+ &
               r1xl*Qu1zl+r2xl*Qu2zl+r3xl*Qu3zl) + &
                 mul*(pmlzxl*duxdxl+pmlzyl*duxdyl+pmlzzl*duxdzl+ &
               r1zl*Qu1xl+r2zl*Qu2xl+r3zl*Qu3xl)
      sigma_yy = lambdal*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
               r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
                 lambdalplus2mul*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
               r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
                 lambdal*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
               r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)
      sigma_yz = mul*(pmlyxl*duzdxl+pmlyyl*duzdyl+pmlyzl*duzdzl+ &
               r1yl*Qu1zl+r2yl*Qu2zl+r3yl*Qu3zl) + &
                 mul*(pmlzxl*duydxl+pmlzyl*duydyl+pmlzzl*duydzl+ &
               r1zl*Qu1yl+r2zl*Qu2yl+r3zl*Qu3yl)
      sigma_zz = lambdal*(pmlxxl*duxdxl+pmlxyl*duxdyl+pmlxzl*duxdzl+ &
               r1xl*Qu1xl+r2xl*Qu2xl+r3xl*Qu3xl) + &
                 lambdal*(pmlyxl*duydxl+pmlyyl*duydyl+pmlyzl*duydzl+ &
               r1yl*Qu1yl+r2yl*Qu2yl+r3yl*Qu3yl) + &
                 lambdalplus2mul*(pmlzxl*duzdxl+pmlzyl*duzdyl+pmlzzl*duzdzl+ &
               r1zl*Qu1zl+r2zl*Qu2zl+r3zl*Qu3zl)
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
  enddo ! iface = 1, num_pml_physical
end subroutine compute_add_pml_physical

subroutine get_gradient_element(s, dx_elem, dy_elem, dz_elem, &
     xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
  use constants, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NDIM
  use specfem_par, only: hprime_xxT,hprime_yyT,hprime_zzT
  implicit none
  real(kind=CUSTOM_REAL), dimension(NDIM, NGLLX, NGLLY, NGLLZ), intent(in) :: s
  real(kind=CUSTOM_REAL), dimension(NDIM, NGLLX, NGLLY, NGLLZ), intent(out) ::&
          dx_elem, dy_elem, dz_elem 
  integer :: i,j,k,l
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL), dimension(NDIM) :: temp1l, temp2l, temp3l
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix,xiy,xiz,etax,&
          etay,etaz,gammax,gammay,gammaz
  dx_elem = 0.0
  dy_elem = 0.0
  dz_elem = 0.0
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        temp1l = 0.0
        temp2l = 0.0
        temp3l = 0.0
        do l = 1, NGLLX
          hp1 = hprime_xxT(l,i)
          temp1l(:) = temp1l(:) + s(:,l,j,k) * hp1
        enddo
        do l = 1, NGLLY
          hp2 = hprime_yyT(l,j)
          temp2l(:) = temp2l(:) + s(:,i,l,k) * hp2
        enddo
        do l = 1, NGLLZ
          hp3 = hprime_zzT(l,k)
          temp3l(:) = temp3l(:) + s(:,i,j,l) * hp3
        enddo
        dx_elem(:,i,j,k)=temp1l(:)*xix(i,j,k)+&
                temp2l(:)*etax(i,j,k)+temp3l(:)*gammax(i,j,k)
        dy_elem(:,i,j,k)=temp1l(:)*xiy(i,j,k)+&
                temp2l(:)*etay(i,j,k)+temp3l(:)*gammay(i,j,k)
        dz_elem(:,i,j,k)=temp1l(:)*xiz(i,j,k)+&
                temp2l(:)*etaz(i,j,k)+temp3l(:)*gammaz(i,j,k)
      enddo
    enddo
  enddo
end subroutine get_gradient_element

subroutine update_pml_aux_conv()
  !use specfem_par, only: deltat, deltatover2
  use pml_par, only: Qt_t,Qu_t,Qt,Qu,Qt_t0,Qu_t0,nglob_CPML,rvolume,&
      coeff_exp1,coeff_exp2,coeff_exp3,&
      coeff_glob_exp1,coeff_glob_exp2,coeff_glob_exp3
  implicit none
  integer :: i
  do i = 1,nglob_CPML
    Qt_t(:,:,i) = Qt_t(:,:,i)*rvolume(i)
  enddo
  do i = 1, 3
    Qu(:,i,:,:,:,:)=Qu(:,i,:,:,:,:)*coeff_exp1(:,:,:,:,:)+&
                  coeff_exp2(:,:,:,:,:)*(Qu_t(:,i,:,:,:,:)+&
                       coeff_exp3(:,:,:,:,:)*Qu_t0(:,i,:,:,:,:))
    Qt(:,i,:)=Qt(:,i,:)*coeff_glob_exp1(:,:)+&
                  coeff_glob_exp2(:,:)*(Qt_t(:,i,:)+&
                       coeff_glob_exp3(:,:)*Qt_t0(:,i,:))
  enddo
  Qt_t0(:,:,:) = Qt_t(:,:,:)
  Qu_t0(:,:,:,:,:,:) = Qu_t(:,:,:,:,:,:)
  Qt_t(:,:,:) = 0.0
  Qu_t(:,:,:,:,:,:) = 0.0
end subroutine update_pml_aux_conv

subroutine update_pml_aux_corrector()
  use specfem_par, only: deltatover2
  use pml_par, only: Qt_t, Qu_t, Qt, Qu, Qt_t0, Qu_t0, nglob_CPML, rvolume
  implicit none
  !! local variables
  integer :: i
  do i = 1,nglob_CPML
    Qt_t(:,:,i) = Qt_t(:,:,i)*rvolume(i)
  enddo
  Qt(:,:,:)=Qt(:,:,:)+deltatover2*(Qt_t(:,:,:)-Qt_t0(:,:,:))
  Qu(:,:,:,:,:,:)=Qu(:,:,:,:,:,:)+deltatover2*(Qu_t(:,:,:,:,:,:)&
                      -Qu_t0(:,:,:,:,:,:))
end subroutine update_pml_aux_corrector

subroutine update_pml_aux_predictor()
  use specfem_par, only: deltatover2, deltat
  use pml_par, only: Qt_t, Qu_t, Qt, Qu, Qt_t0, Qu_t0
  Qt(:,:,:) = Qt(:,:,:) + deltat*Qt_t(:,:,:)
  Qu(:,:,:,:,:,:) = Qu(:,:,:,:,:,:) + deltat*Qu_t(:,:,:,:,:,:)
  Qt_t0(:,:,:) = Qt_t(:,:,:)
  Qu_t0(:,:,:,:,:,:) = Qu_t(:,:,:,:,:,:)
  Qt_t(:,:,:) = 0.0
  Qu_t(:,:,:,:,:,:) = 0.0
end subroutine update_pml_aux_predictor

subroutine prepare_adepml_update_coeff()
  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
  use specfem_par, only: deltat
  use pml_par, only: coeff_exp1,coeff_exp2,coeff_exp3,coeff_glob_exp1,&
    coeff_glob_exp2,coeff_glob_exp3,pml_beta,nspec_CPML,ibool_CPML
  implicit none
  integer :: i,j,k,iglob,ispec
  coeff_exp1(:,:,:,:,:)=exp(-pml_beta(:,:,:,:,:)*deltat)
  coeff_exp3(:,:,:,:,:)=sqrt(coeff_exp1(:,:,:,:,:))
  coeff_exp2(:,:,:,:,:)=(1._CUSTOM_REAL-coeff_exp3(:,:,:,:,:))/pml_beta(:,:,:,:,:)
  do ispec = 1, nspec_CPML
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_CPML(i,j,k,ispec)
          coeff_glob_exp1(:,iglob)=coeff_exp1(:,i,j,k,ispec)
          coeff_glob_exp2(:,iglob)=coeff_exp2(:,i,j,k,ispec)
          coeff_glob_exp3(:,iglob)=coeff_exp3(:,i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo
  if (maxval(coeff_exp1) > 1.0) print *, 'bad value coeff_exp1'
  if (maxval(coeff_exp3) > 1.0) print *, 'bad value coeff_exp3'
  if (maxval(coeff_glob_exp1) > 1.0) print *, 'bad value coeff_glob_exp1'
  if (maxval(coeff_glob_exp3) > 1.0) print *, 'bad value coeff_glob_exp3'
end subroutine prepare_adepml_update_coeff
