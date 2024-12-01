subroutine compute_forces_viscoelastic_adepml_element(ispec,displ,accel)
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD,FOUR_THIRDS
  use specfem_par, only: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
               NGLOB_AB,hprime_xxT,hprime_yyT,hprime_zzT, &
               hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
               wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
               kappastore,mustore,jacobian,ibool,&
               irregular_element_number,xix_regular,jacobian_regular, &
               t_force_pml, hprime_xx ! nqdu
  use pml_par, only: spec_to_CPML, accel_elastic_CPML
  implicit none
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ,accel
  ! local parameters
  integer :: ispec, iglob, ispec_irreg
  integer :: i,j,k,l
  integer :: ispec_CPML

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  !real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  !real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  
  !real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal, jacobianl

  real(kind=CUSTOM_REAL), dimension(6,NGLLX,NGLLY,NGLLZ) :: sigma_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            dummyx_loc,dummyy_loc,dummyz_loc, &
            tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
          
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: zero_array
  real(kind=CUSTOM_REAL), dimension(6,3,3,NGLLX,NGLLY,NGLLZ) :: temp_adepml
  real(kind=CUSTOM_REAL), dimension(6,3,3) :: newtemp_adepml
  double precision, external :: wtime
  double precision :: t_clock

  ! TL: timing
  t_clock = wtime()

  !nqdu
  ! generate an array equal to zero
  !zero_array(:,:,:) = 0._CUSTOM_REAL

  ispec_CPML = spec_to_CPML(ispec)

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        iglob = ibool(i,j,k,ispec)
        dummyx_loc(i,j,k) = displ(1,iglob)
        dummyy_loc(i,j,k) = displ(2,iglob)
        dummyz_loc(i,j,k) = displ(3,iglob)
      enddo
    enddo
  enddo
  ! dv^{n+1/2}/dx
  !nqdu
  ! call compute_strain_in_element( &
  !                tempx1,tempx2,tempx3,zero_array,zero_array,zero_array, &
  !                tempy1,tempy2,tempy3,zero_array,zero_array,zero_array, &
  !                tempz1,tempz2,tempz3,zero_array,zero_array,zero_array, &
  !                dummyx_loc,dummyy_loc,dummyz_loc, &
  !                hprime_xxT,hprime_yyT,hprime_zzT)
  call compute_strain_in_element( &
                 tempx1,tempx2,tempx3, &
                 tempy1,tempy2,tempy3, &
                 tempz1,tempz2,tempz3, &
                 dummyx_loc,dummyy_loc,dummyz_loc, &
                 hprime_xxT,hprime_xx)
  ispec_irreg = irregular_element_number(ispec)
  if (ispec_irreg == 0) jacobianl = jacobian_regular

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
          jacobianl = jacobian(i,j,k,ispec_irreg)

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
        ! precompute some sums to save CPU time
        !duxdxl_plus_duydyl = duxdxl + duydyl
        !duxdxl_plus_duzdzl = duxdxl + duzdzl
        !duydyl_plus_duzdzl = duydyl + duzdzl
        !duxdyl_plus_duydxl = duxdyl + duydxl
        !duzdxl_plus_duxdzl = duzdxl + duxdzl
        !duzdyl_plus_duydzl = duzdyl + duydzl

        kappal = kappastore(i,j,k,ispec)
        mul = mustore(i,j,k,ispec)

        ! isotropic case
        lambdalplus2mul = kappal + FOUR_THIRDS * mul
        lambdal = lambdalplus2mul - 2._CUSTOM_REAL * mul
        ! d\sigma^{n+1/2}/dt = ... dv^{n+1/2}/dx+r*Qu^{n+1/2}
        call compute_sigma_adepml(ispec_CPML,i,j,k,&
                                lambdalplus2mul,lambdal,mul, &
                                sigma_elem(:,i,j,k), &
                                duxdxl, duydxl, duzdxl, &
                                duxdyl, duydyl, duzdyl, &
                                duxdzl, duydzl, duzdzl)
        call compute_Qu_t_point(ispec_CPML,i,j,k,&
                              duxdxl, duydxl, duzdxl, &
                              duxdyl, duydyl, duzdyl, &
                              duxdzl, duydzl, duzdzl)
        ! \sigma^{n+1/2}=\sigma^{n}+1/2*\delta t*d\sigma^{n+1/2}/dt
        !sigma_mid_elem(:,i,j,k) = pml_sigma(:,i,j,k,ispec_CPML) + &
        !           dsigma_elem(:,i,j,k) * deltatover2
        ! \sigma^{n+1}=\sigma^{n+1/2}+1/2*\delta t*d\sigma^{n+1/2}/dt
        !pml_sigma(:,i,j,k,ispec_CPML) = pml_sigma(:,i,j,k,ispec_CPML) + &
        !           dsigma_elem(:,i,j,k) * deltat
        !pml_sigma(2,i,j,k,ispec_CPML) = pml_sigma(2,i,j,k,ispec_CPML) + &
        !           dsigma_elem(2,i,j,k) * deltatover2 
        !pml_sigma(3,i,j,k,ispec_CPML) = pml_sigma(3,i,j,k,ispec_CPML) + &
        !           dsigma_elem(3,i,j,k) * deltatover2
        !pml_sigma(4,i,j,k,ispec_CPML) = pml_sigma(4,i,j,k,ispec_CPML) + &
        !           dsigma_elem(4,i,j,k) * deltatover2
        !pml_sigma(5,i,j,k,ispec_CPML) = pml_sigma(5,i,j,k,ispec_CPML) + &
        !           dsigma_elem(5,i,j,k) * deltatover2
        !pml_sigma(6,i,j,k,ispec_CPML) = pml_sigma(6,i,j,k,ispec_CPML) + &
        !           dsigma_elem(6,i,j,k) * deltatover2
      enddo
    enddo
  enddo
  !do i = 1, 6
  !  temp_adepml(i,1,1,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*xix(:,:,:,ispec_irreg)
  !  temp_adepml(i,2,1,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*xiy(:,:,:,ispec_irreg)
  !  temp_adepml(i,3,1,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*xiz(:,:,:,ispec_irreg)
  !  temp_adepml(i,1,2,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*etax(:,:,:,ispec_irreg)
  !  temp_adepml(i,2,2,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*etay(:,:,:,ispec_irreg)
  !  temp_adepml(i,3,2,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*etaz(:,:,:,ispec_irreg)
  !  temp_adepml(i,1,3,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*gammax(:,:,:,ispec_irreg)
  !  temp_adepml(i,2,3,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*gammay(:,:,:,ispec_irreg)
  !  temp_adepml(i,3,3,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*gammaz(:,:,:,ispec_irreg)
  !enddo
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        jacobianl=jacobian(i,j,k,ispec_irreg)
        xixl=xix(i,j,k,ispec_irreg)*jacobianl
        xiyl = xiy(i,j,k,ispec_irreg)*jacobianl
        xizl = xiz(i,j,k,ispec_irreg)*jacobianl
        etaxl = etax(i,j,k,ispec_irreg)*jacobianl
        etayl = etay(i,j,k,ispec_irreg)*jacobianl
        etazl = etaz(i,j,k,ispec_irreg)*jacobianl
        gammaxl = gammax(i,j,k,ispec_irreg)*jacobianl
        gammayl = gammay(i,j,k,ispec_irreg)*jacobianl
        gammazl = gammaz(i,j,k,ispec_irreg)*jacobianl
       
        temp_adepml(:,1,1,i,j,k)=sigma_elem(:,i,j,k)*xixl
        !temp_adepml_mid(:,1,1,i,j,k)=sigma_mid_elem(:,i,j,k)*xixl

        temp_adepml(:,2,1,i,j,k)=sigma_elem(:,i,j,k)*xiyl
        !temp_adepml_mid(:,2,1,i,j,k)=sigma_mid_elem(:,i,j,k)*xiyl
 
        temp_adepml(:,3,1,i,j,k)=sigma_elem(:,i,j,k)*xizl
        !temp_adepml_mid(:,3,1,i,j,k)=sigma_mid_elem(:,i,j,k)*xizl

        temp_adepml(:,1,2,i,j,k)=sigma_elem(:,i,j,k)*etaxl
        !temp_adepml_mid(:,1,2,i,j,k)=sigma_mid_elem(:,i,j,k)*etaxl

        temp_adepml(:,2,2,i,j,k)=sigma_elem(:,i,j,k)*etayl
        !temp_adepml_mid(:,2,2,i,j,k)=sigma_mid_elem(:,i,j,k)*etayl

        temp_adepml(:,3,2,i,j,k)=sigma_elem(:,i,j,k)*etazl
        !temp_adepml_mid(:,3,2,i,j,k)=sigma_mid_elem(:,i,j,k)*etazl

        temp_adepml(:,1,3,i,j,k)=sigma_elem(:,i,j,k)*gammaxl
        !temp_adepml_mid(:,1,3,i,j,k)=sigma_mid_elem(:,i,j,k)*gammaxl

        temp_adepml(:,2,3,i,j,k)=sigma_elem(:,i,j,k)*gammayl
        !temp_adepml_mid(:,2,3,i,j,k)=sigma_mid_elem(:,i,j,k)*gammayl

        temp_adepml(:,3,3,i,j,k)=sigma_elem(:,i,j,k)*gammazl
        !temp_adepml_mid(:,3,3,i,j,k)=sigma_mid_elem(:,i,j,k)*gammazl
        
      enddo
    enddo
  enddo

  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        !jacobianl = jacobian(i,j,k,ispec_irreg)
        newtemp_adepml(:,:,:) = 0._CUSTOM_REAL
        !newtemp_adepml_mid(:,:,:) = 0._CUSTOM_REAL
        do l=1,NGLLX
          fac1 = hprimewgll_xx(l,i)
          newtemp_adepml(:,:,1) = newtemp_adepml(:,:,1) + &
                                temp_adepml(:,:,1,l,j,k)*fac1
          !newtemp_adepml_mid(:,:,1) = newtemp_adepml_mid(:,:,1) + &
          !                      temp_adepml_mid(:,:,1,l,j,k)*fac1
          fac2 = hprimewgll_yy(l,j)
          newtemp_adepml(:,:,2) = newtemp_adepml(:,:,2) + &
                                temp_adepml(:,:,2,i,l,k)*fac2
          !newtemp_adepml_mid(:,:,2) = newtemp_adepml_mid(:,:,2) + &
          !                      temp_adepml_mid(:,:,2,i,l,k)*fac2
          fac3 = hprimewgll_zz(l,k)
          newtemp_adepml(:,:,3) = newtemp_adepml(:,:,3) + &
                                temp_adepml(:,:,3,i,j,l)*fac3
          !newtemp_adepml_mid(:,:,3) = newtemp_adepml_mid(:,:,3) + &
          !                      temp_adepml_mid(:,:,3,i,j,l)*fac3
        enddo
        fac1 = wgllwgll_yz(j,k)
        fac2 = wgllwgll_xz(i,k)
        fac3 = wgllwgll_xy(i,j)
        !wgll3_xyzl = wgll3_xyz(i,j,k)
        ! dQt^{n+1/2}/dt=...\sigma^{n+1/2} (volume integral)
        call compute_Qt_t(ispec_CPML,i,j,k,newtemp_adepml,fac1,fac2,fac3)
        ! a^{n+1}=-...\sigma^{n+1}
        call compute_accel_adepml(ispec_CPML,i,j,k,newtemp_adepml,fac1,fac2,fac3)
      enddo
    enddo
  enddo
  ! dQt^{n+1/2}/dt=...\sigma^{n+1/2} (surface integral)
  call compute_add_pml_physical_element(ispec_CPML,sigma_elem)  

  !! \sigma^{n+1}=\sigma^{n+1/2}+1/2*\delta t*d\sigma^{n+1/2}/dt
  !pml_sigma(:,:,:,:,ispec_CPML) = pml_sigma(:,:,:,:,ispec_CPML) + &
  !                dsigma_elem(:,:,:,:) * deltatover2
  !do i = 1, 6
  !  temp_adepml(i,1,1,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*xix(:,:,:,ispec_irreg)
  !  temp_adepml(i,2,1,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*xiy(:,:,:,ispec_irreg)
  !  temp_adepml(i,3,1,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*xiz(:,:,:,ispec_irreg)
  !  temp_adepml(i,1,2,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*etax(:,:,:,ispec_irreg)
  !  temp_adepml(i,2,2,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*etay(:,:,:,ispec_irreg)
  !  temp_adepml(i,3,2,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*etaz(:,:,:,ispec_irreg)
  !  temp_adepml(i,1,3,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*gammax(:,:,:,ispec_irreg)
  !  temp_adepml(i,2,3,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*gammay(:,:,:,ispec_irreg)
  !  temp_adepml(i,3,3,:,:,:) = jacobian(:,:,:,ispec_irreg)*&
  !                pml_sigma(i,:,:,:,ispec_CPML)*gammaz(:,:,:,ispec_irreg)
  !enddo
!  do k = 1, NGLLZ
!    do j = 1, NGLLY
!      do i = 1, NGLLX
!        do l = 1, 6
!          temp_adepml(l,1,1,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*xix(i,j,k,ispec_irreg)
!        enddo
!        do l = 1, 6
!          temp_adepml(l,2,1,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*xiy(i,j,k,ispec_irreg)
!        enddo
!        do l = 1, 6
!          temp_adepml(l,3,1,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*xiz(i,j,k,ispec_irreg)
!        enddo
!        do l = 1, 6
!          temp_adepml(l,1,2,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*etax(i,j,k,ispec_irreg)
!        enddo
!        do l = 1, 6
!          temp_adepml(l,2,2,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*etay(i,j,k,ispec_irreg)
!        enddo
!        do l = 1, 6
!          temp_adepml(l,3,2,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*etaz(i,j,k,ispec_irreg)
!        enddo
!        do l = 1, 6
!          temp_adepml(l,1,3,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*gammax(i,j,k,ispec_irreg)
!        enddo
!        do l = 1, 6
!          temp_adepml(l,2,3,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*gammay(i,j,k,ispec_irreg)
!        enddo
!        do l = 1, 6
!          temp_adepml(l,3,3,i,j,k) = jacobian(i,j,k,ispec_irreg)*&
!                  pml_sigma(l,i,j,k,ispec_CPML)*gammaz(i,j,k,ispec_irreg)
!        enddo
!      enddo
!    enddo
!  enddo
!  
!  do k = 1, NGLLZ
!    do j = 1, NGLLY
!      do i = 1, NGLLX
!        !jacobianl = jacobian(i,j,k,ispec_irreg)
!        newtemp_adepml(:,:,:) = 0._CUSTOM_REAL
!        do l=1,NGLLX
!          fac1 = hprimewgll_xx(l,i)
!          newtemp_adepml(:,:,1) = newtemp_adepml(:,:,1) + &
!                                temp_adepml(:,:,1,l,j,k)*fac1
!          fac2 = hprimewgll_yy(l,j)
!          newtemp_adepml(:,:,2) = newtemp_adepml(:,:,2) + &
!                                temp_adepml(:,:,2,i,l,k)*fac2
!          fac3 = hprimewgll_zz(l,k)
!          newtemp_adepml(:,:,3) = newtemp_adepml(:,:,3) + &
!                                temp_adepml(:,:,3,i,j,l)*fac3
!        enddo
!        fac1 = wgllwgll_yz(j,k)
!        fac2 = wgllwgll_xz(i,k)
!        fac3 = wgllwgll_xy(i,j)
!        !wgll3_xyzl = wgll3_xyz(i,j,k)
!        ! a^{n+1}=-...\sigma^{n+1}
!        call compute_accel_adepml(ispec_CPML,i,j,k,newtemp_adepml,fac1,fac2,fac3)
!      enddo
!    enddo
!  enddo

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)
        accel(1,iglob) = accel(1,iglob) - accel_elastic_CPML(1,i,j,k)
        accel(2,iglob) = accel(2,iglob) - accel_elastic_CPML(2,i,j,k)
        accel(3,iglob) = accel(3,iglob) - accel_elastic_CPML(3,i,j,k)
      enddo
    enddo
  enddo  
  ! TL: timing
  t_force_pml = t_force_pml + (wtime() - t_clock)
  
  !nqdu
  ! contains
  !  subroutine compute_strain_in_element(tempx1_att,tempx2_att,tempx3_att,tempx1,tempx2,tempx3, &
  !                                           tempy1_att,tempy2_att,tempy3_att,tempy1,tempy2,tempy3, &
  !                                           tempz1_att,tempz2_att,tempz3_att,tempz1,tempz2,tempz3, &
  !                                           dummyx_loc,dummyy_loc,dummyz_loc,hprime_xxT,hprime_yyT,hprime_zzT)

  ! use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_att,tempx2_att,tempx3_att, &
  !                                                         tempy1_att,tempy2_att,tempy3_att, &
  !                                                         tempz1_att,tempz2_att,tempz3_att

  ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3,&
  !                                                         tempy1,tempy2,tempy3,&
  !                                                         tempz1,tempz2,tempz3

  ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT
  ! real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yyT
  ! real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zzT

  ! ! local variables
  ! integer :: i,j,k,l
  ! real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  ! tempx1_att(:,:,:) = tempx1(:,:,:)
  ! tempx2_att(:,:,:) = tempx2(:,:,:)
  ! tempx3_att(:,:,:) = tempx3(:,:,:)

  ! tempy1_att(:,:,:) = tempy1(:,:,:)
  ! tempy2_att(:,:,:) = tempy2(:,:,:)
  ! tempy3_att(:,:,:) = tempy3(:,:,:)

  ! tempz1_att(:,:,:) = tempz1(:,:,:)
  ! tempz2_att(:,:,:) = tempz2(:,:,:)
  ! tempz3_att(:,:,:) = tempz3(:,:,:)
  
  ! ! use first order Taylor expansion of displacement for local storage of
  ! ! stresses
  ! ! at this current time step, to fix attenuation in a consistent way
  ! do k=1,NGLLZ
  !   do j=1,NGLLY
  !     do i=1,NGLLX

  !       ! we can merge these loops because NGLLX = NGLLY = NGLLZ
  !       do l=1,NGLLX
  !         hp1 = hprime_xxT(l,i)
  !         tempx1_att(i,j,k) = tempx1_att(i,j,k) + dummyx_loc(l,j,k) * hp1
  !         tempy1_att(i,j,k) = tempy1_att(i,j,k) + dummyy_loc(l,j,k) * hp1
  !         tempz1_att(i,j,k) = tempz1_att(i,j,k) + dummyz_loc(l,j,k) * hp1

  !         hp2 = hprime_yyT(l,j)
  !         tempx2_att(i,j,k) = tempx2_att(i,j,k) + dummyx_loc(i,l,k) * hp2
  !         tempy2_att(i,j,k) = tempy2_att(i,j,k) + dummyy_loc(i,l,k) * hp2
  !         tempz2_att(i,j,k) = tempz2_att(i,j,k) + dummyz_loc(i,l,k) * hp2

  !         hp3 = hprime_zzT(l,k)
  !         tempx3_att(i,j,k) = tempx3_att(i,j,k) + dummyx_loc(i,j,l) * hp3
  !         tempy3_att(i,j,k) = tempy3_att(i,j,k) + dummyy_loc(i,j,l) * hp3
  !         tempz3_att(i,j,k) = tempz3_att(i,j,k) + dummyz_loc(i,j,l) * hp3
  !       enddo

  !     enddo
  !   enddo
  ! enddo

  ! end subroutine compute_strain_in_element
end subroutine compute_forces_viscoelastic_adepml_element

