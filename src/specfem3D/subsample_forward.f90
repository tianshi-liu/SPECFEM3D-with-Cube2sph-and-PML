subroutine open_forward_wavefield_write()
  use specfem_par, only: NGLOB_AB,NDIM,IFILE_FORWARD_WAVEFIELD,&
                         prname,CUSTOM_REAL,MAX_STRING_LEN
  implicit none
  integer :: ier
  character (len=MAX_STRING_LEN) :: file_name
  file_name = prname(1:len_trim(prname)) // 'forward_wavefield.bin'
  open(unit=IFILE_FORWARD_WAVEFIELD, file=file_name(1:len_trim(file_name)), &
       access='direct', recl=NDIM*NGLOB_AB*CUSTOM_REAL, &
       action='write',form='unformatted', iostat=ier) 
  if (ier /= 0) &
       call exit_MPI_without_rank('error opening to write:'// file_name)
end subroutine open_forward_wavefield_write

subroutine open_forward_wavefield_read()
  use specfem_par, only: NGLOB_AB,NDIM,IFILE_FORWARD_WAVEFIELD,&
                         prname,CUSTOM_REAL,MAX_STRING_LEN
  implicit none
  integer :: ier
  character (len=MAX_STRING_LEN) :: file_name
  file_name = prname(1:len_trim(prname)) // 'forward_wavefield.bin'
  open(unit=IFILE_FORWARD_WAVEFIELD, file=file_name(1:len_trim(file_name)), &
       access='direct', recl=NDIM*NGLOB_AB*CUSTOM_REAL, &
       action='read',form='unformatted', iostat=ier) 
  if (ier /= 0) &
       call exit_MPI_without_rank('error opening to read:'// file_name)
end subroutine open_forward_wavefield_read

subroutine close_forward_wavefield()
  use specfem_par, only: IFILE_FORWARD_WAVEFIELD
  close(IFILE_FORWARD_WAVEFIELD)
end subroutine close_forward_wavefield

subroutine write_subsampled_forward_wavefield(it)
  use specfem_par, only: IFILE_FORWARD_WAVEFIELD, NSTEP_PER_FORWARD_OUTPUT
  use specfem_par_elastic, only: displ

  implicit none
  integer :: i_save, it

  i_save = int(it / NSTEP_PER_FORWARD_OUTPUT)

  write(IFILE_FORWARD_WAVEFIELD, rec=i_save) displ
                         
end subroutine write_subsampled_forward_wavefield

subroutine read_subsampled_forward_wavefield(b_it)
  use specfem_par, only: IFILE_FORWARD_WAVEFIELD, NSTEP_PER_FORWARD_OUTPUT
  use specfem_par_elastic, only: b_displ
  implicit none
  integer :: i_save, b_it

  i_save = int(b_it / NSTEP_PER_FORWARD_OUTPUT)

  read(IFILE_FORWARD_WAVEFIELD, rec=i_save) b_displ
  
end subroutine read_subsampled_forward_wavefield

subroutine compute_kernels_from_subsampled_wavefield()
  use specfem_par
  use specfem_par_elastic

  implicit none
  ! local parameters
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL),dimension(21) :: prod
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc,b_epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    epsilon_trace_over_3_elem,epsilondev_xx_elem, epsilondev_yy_elem, &
    epsilondev_xy_elem, epsilondev_xz_elem, epsilondev_yz_elem, &
    b_epsilon_trace_over_3_elem,b_epsilondev_xx_elem, b_epsilondev_yy_elem, &
    b_epsilondev_xy_elem, b_epsilondev_xz_elem, b_epsilondev_yz_elem
  !! cannot do GPU_MODE or ANISOTROPIC_VELOCITY_KL yet
  do ispec = 1, NSPEC_AB
    ! compute strain in each element for forward and adjoint field
    call compute_epsilon_element(displ,ispec, &
       epsilon_trace_over_3_elem,epsilondev_xx_elem, epsilondev_yy_elem, &
       epsilondev_xy_elem, epsilondev_xz_elem, epsilondev_yz_elem)
    call compute_epsilon_element(b_displ,ispec, &
       b_epsilon_trace_over_3_elem,b_epsilondev_xx_elem, b_epsilondev_yy_elem,&
       b_epsilondev_xy_elem, b_epsilondev_xz_elem, b_epsilondev_yz_elem)
    do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
      iglob = ibool(i,j,k,ispec)
      epsilondev_loc(1) = epsilondev_xx_elem(i,j,k)
      epsilondev_loc(2) = epsilondev_yy_elem(i,j,k)
      epsilondev_loc(3) = epsilondev_xy_elem(i,j,k)
      epsilondev_loc(4) = epsilondev_xz_elem(i,j,k)
      epsilondev_loc(5) = epsilondev_yz_elem(i,j,k)

      b_epsilondev_loc(1) = b_epsilondev_xx_elem(i,j,k)
      b_epsilondev_loc(2) = b_epsilondev_yy_elem(i,j,k)
      b_epsilondev_loc(3) = b_epsilondev_xy_elem(i,j,k)
      b_epsilondev_loc(4) = b_epsilondev_xz_elem(i,j,k)
      b_epsilondev_loc(5) = b_epsilondev_yz_elem(i,j,k)

      rho_kl(i,j,k,ispec) =  rho_kl(i,j,k,ispec)+deltat* &
        NSTEP_PER_FORWARD_OUTPUT*dot_product(accel(:,iglob),b_displ(:,iglob))
      ! For anisotropic kernels
      if (ANISOTROPIC_KL) then
        call compute_strain_product(prod,epsilon_trace_over_3_elem(i,j,k),&
             epsilondev_loc,b_epsilon_trace_over_3_elem(i,j,k),b_epsilondev_loc)
        cijkl_kl(:,i,j,k,ispec) = cijkl_kl(:,i,j,k,ispec)+deltat *&
          NSTEP_PER_FORWARD_OUTPUT*prod(:)
      else
        mu_kl(i,j,k,ispec)=mu_kl(i,j,k,ispec)+deltat*NSTEP_PER_FORWARD_OUTPUT*&
         (epsilondev_loc(1)*b_epsilondev_loc(1)+&
            epsilondev_loc(2)*b_epsilondev_loc(2)+&
            (epsilondev_loc(1)+epsilondev_loc(2))*&
            (b_epsilondev_loc(1)+b_epsilondev_loc(2))+&
            2 * (epsilondev_loc(3)*b_epsilondev_loc(3)+&
            epsilondev_loc(4)*b_epsilondev_loc(4) + &
            epsilondev_loc(5)*b_epsilondev_loc(5)) )
        kappa_kl(i,j,k,ispec)=kappa_kl(i,j,k,ispec)+deltat*&
          NSTEP_PER_FORWARD_OUTPUT *&
            (9 * epsilon_trace_over_3_elem(i,j,k) &
            * b_epsilon_trace_over_3_elem(i,j,k))
      endif
    enddo;enddo;enddo
  enddo ! ispec = 1, NSPEC_AB
end subroutine compute_kernels_from_subsampled_wavefield

subroutine compute_epsilon_element(u,ispec,e,exx,eyy,exy,exz,eyz)
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD,FOUR_THIRDS
  use specfem_par, only: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                         NGLOB_AB,hprime_xxT,hprime_yyT,hprime_zzT, &
                         ibool,irregular_element_number,xix_regular
  implicit none
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: u
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: e,exx,eyy,exy,exz,eyz
  integer :: ispec, iglob, ispec_irreg
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,&
                            duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,&
                            dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,&
                            gammaxl,gammayl,gammazl,templ
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: zero_array
  zero_array(:,:,:) = 0._CUSTOM_REAL
  do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
    iglob = ibool(i,j,k,ispec)
    dummyx_loc(i,j,k) = u(1,iglob)
    dummyy_loc(i,j,k) = u(2,iglob)
    dummyz_loc(i,j,k) = u(3,iglob)
  enddo;enddo;enddo
  call compute_strain_in_element( &
                 tempx1,tempx2,tempx3,zero_array,zero_array,zero_array, &
                 tempy1,tempy2,tempy3,zero_array,zero_array,zero_array, &
                 tempz1,tempz2,tempz3,zero_array,zero_array,zero_array, &
                 dummyx_loc,dummyy_loc,dummyz_loc, &
                 hprime_xxT,hprime_yyT,hprime_zzT)
  ispec_irreg = irregular_element_number(ispec)
  do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
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
    templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
    e(i,j,k) = templ
    exx(i,j,k) = duxdxl - templ
    eyy(i,j,k) = duydyl - templ
    exy(i,j,k) = 0.5_CUSTOM_REAL * (duxdyl + duydxl)
    exz(i,j,k) = 0.5_CUSTOM_REAL * (duzdxl + duxdzl)
    eyz(i,j,k) = 0.5_CUSTOM_REAL * (duzdyl + duydzl)
  enddo;enddo;enddo

  contains

  subroutine compute_strain_in_element(tempx1_att,tempx2_att,tempx3_att,tempx1,tempx2,tempx3, &
                                            tempy1_att,tempy2_att,tempy3_att,tempy1,tempy2,tempy3, &
                                            tempz1_att,tempz2_att,tempz3_att,tempz1,tempz2,tempz3, &
                                            dummyx_loc,dummyy_loc,dummyz_loc,hprime_xxT,hprime_yyT,hprime_zzT)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_att,tempx2_att,tempx3_att, &
                                                          tempy1_att,tempy2_att,tempy3_att, &
                                                          tempz1_att,tempz2_att,tempz3_att

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3, &
                                                          tempy1,tempy2,tempy3, &
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
end subroutine compute_epsilon_element
