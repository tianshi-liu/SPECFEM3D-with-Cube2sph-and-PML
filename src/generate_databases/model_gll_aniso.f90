!--------------------------------------------------------------------------------------------------
!! nqdu ANISO GLL Model C21
!
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!
!-----------

subroutine model_gll_aniso(myrank,nspec,LOCAL_PATH)
  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,&
    MAX_STRING_LEN,ATTENUATION,CUSTOM_REAL,ibool

  use create_regions_mesh_ext_par, only: rhostore,kappastore,mustore,rho_vp,rho_vs,&
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store,eta_anistore, &
            xstore_dummy, ystore_dummy, zstore_dummy
  use create_regions_mesh_ext_par, only: CUBE2SPH_MESH, &
               kappavstore,muvstore

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:), allocatable :: vp_read,vs_read,rho_read
  
  integer :: ier,ispec,i,j,k,iglob
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26,&
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  character(len=MAX_STRING_LEN) :: prname_lp,filename
  real(kind=CUSTOM_REAL) :: xp, yp, zp, r_dummy, theta, phi

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using GLL model from: ',trim(LOCAL_PATH)
  endif

  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)// '/' //'proc',myrank,'_'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! if only vp structure is available (as is often the case in exploration seismology),
  !!! use lines for vp only

  ! density
  allocate(rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 647')
  if (ier /= 0) stop 'error allocating array rho_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: rho.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'rho.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading rho.bin file'
  endif

  read(28) rho_read
  close(28)

  ! vp
  allocate(vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
  if (ier /= 0) stop 'error allocating array vp_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: vp.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vp.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vp.bin file'
  endif

  read(28) vp_read
  close(28)

  ! vs
  allocate(vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 649')
  if (ier /= 0) stop 'error allocating array vs_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: vs.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vs.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vs.bin file'
  endif

  read(28) vs_read
  close(28)

  ! C21
  if (myrank == 0) write(IMAIN,*) '     reading in: radial cijkl.bin'

  filename = prname_lp(1:len_trim(prname_lp))// 'c11.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c11.bin'
  read(28) c11store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c12.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c12.bin'
  read(28) c12store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c13.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c13.bin'
  read(28) c13store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c14.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c14.bin'
  read(28) c14store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c15.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c15.bin'
  read(28) c15store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c16.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c16.bin'
  read(28) c16store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c22.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c22.bin'
  read(28) c22store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c23.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c23.bin'
  read(28) c23store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c24.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c24.bin'
  read(28) c24store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c25.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c25.bin'
  read(28) c25store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c26.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c26.bin'
  read(28) c26store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c33.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c33.bin'
  read(28) c33store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c34.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c34.bin'
  read(28) c34store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c35.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c35.bin'
  read(28) c35store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c36.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c36.bin'
  read(28) c36store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c44.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c44.bin'
  read(28) c44store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c45.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c45.bin'
  read(28) c45store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c46.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c46.bin'
  read(28) c46store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c55.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c55.bin'
  read(28) c55store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c56.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c56.bin'
  read(28) c56store
  close(28)

  filename = prname_lp(1:len_trim(prname_lp))// 'c66.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file c66.bin'
  read(28) c66store
  close(28)

  eta_anistore = 1.

  if (CUBE2SPH_MESH) then
    do ispec = 1, nspec
      do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX
        d11 = c11store(i,j,k,ispec)
        d12 = c12store(i,j,k,ispec)
        d13 = c13store(i,j,k,ispec)
        d14 = c14store(i,j,k,ispec)
        d15 = c15store(i,j,k,ispec)
        d16 = c16store(i,j,k,ispec)
        d22 = c22store(i,j,k,ispec)
        d23 = c23store(i,j,k,ispec)
        d24 = c24store(i,j,k,ispec)
        d25 = c25store(i,j,k,ispec)
        d26 = c26store(i,j,k,ispec)
        d33 = c33store(i,j,k,ispec)
        d34 = c34store(i,j,k,ispec)
        d35 = c35store(i,j,k,ispec)
        d36 = c36store(i,j,k,ispec)
        d44 = c44store(i,j,k,ispec)
        d45 = c45store(i,j,k,ispec)
        d46 = c46store(i,j,k,ispec)
        d55 = c55store(i,j,k,ispec)
        d56 = c56store(i,j,k,ispec)
        d66 = c66store(i,j,k,ispec)

        iglob = ibool(i,j,k,ispec)
        xp = xstore_dummy(iglob)
        yp = ystore_dummy(iglob)
        zp = zstore_dummy(iglob)

        ! now z axis is in r direction
        call xyz_2_rthetaphi(xp,yp,zp,r_dummy,theta,phi)
        
        ! nqdu add
        call rotate_tensor_radial_to_global(&
            dble(theta),dble(phi),d11,d12,d13,d14,d15,d16, &
            d22,d23,d24,d25,d26, &
            d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
            c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66 &
        )

        ! call rotate_aniso_tensor(dble(theta),dble(phi),d11,d12,d13,d14,d15,d16, &
        !                     d22,d23,d24,d25,d26, &
        !                     d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
        !                     c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
        !                     c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

        c11store(i,j,k,ispec) = real(c11, kind=CUSTOM_REAL)
        c12store(i,j,k,ispec) = real(c12, kind=CUSTOM_REAL)
        c13store(i,j,k,ispec) = real(c13, kind=CUSTOM_REAL)
        c14store(i,j,k,ispec) = real(c14, kind=CUSTOM_REAL)
        c15store(i,j,k,ispec) = real(c15, kind=CUSTOM_REAL)
        c16store(i,j,k,ispec) = real(c16, kind=CUSTOM_REAL)
        c22store(i,j,k,ispec) = real(c22, kind=CUSTOM_REAL)
        c23store(i,j,k,ispec) = real(c23, kind=CUSTOM_REAL)
        c24store(i,j,k,ispec) = real(c24, kind=CUSTOM_REAL)
        c25store(i,j,k,ispec) = real(c25, kind=CUSTOM_REAL)
        c26store(i,j,k,ispec) = real(c26, kind=CUSTOM_REAL)
        c33store(i,j,k,ispec) = real(c33, kind=CUSTOM_REAL)
        c34store(i,j,k,ispec) = real(c34, kind=CUSTOM_REAL)
        c35store(i,j,k,ispec) = real(c35, kind=CUSTOM_REAL)
        c36store(i,j,k,ispec) = real(c36, kind=CUSTOM_REAL)
        c44store(i,j,k,ispec) = real(c44, kind=CUSTOM_REAL)
        c45store(i,j,k,ispec) = real(c45, kind=CUSTOM_REAL)
        c46store(i,j,k,ispec) = real(c46, kind=CUSTOM_REAL)
        c55store(i,j,k,ispec) = real(c55, kind=CUSTOM_REAL)
        c56store(i,j,k,ispec) = real(c56, kind=CUSTOM_REAL)
        c66store(i,j,k,ispec) = real(c66, kind=CUSTOM_REAL)
      enddo;enddo;enddo;
    enddo
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where shear wavespeed structure is not given
  !!! modify according to your desire

  !   vs_read = 0.0
  !   where ( mustore > 100.0 )       vs_read = vp_read / sqrt(3.0)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! update arrays that will be saved and used in the solver xspecfem3D
  !!! the following part is neccessary if you uncommented something above

  rhostore(:,:,:,:) = rho_read(:,:,:,:)
  kappastore(:,:,:,:) = rhostore(:,:,:,:) * ( vp_read(:,:,:,:) * vp_read(:,:,:,:) &
                                              - FOUR_THIRDS * vs_read(:,:,:,:) * vs_read(:,:,:,:) )
  mustore(:,:,:,:) = rhostore(:,:,:,:) * vs_read(:,:,:,:) * vs_read(:,:,:,:)
  kappavstore = kappastore
  muvstore = mustore
  rho_vp(:,:,:,:) = rhostore(:,:,:,:) * vp_read(:,:,:,:)
  rho_vs(:,:,:,:) = rhostore(:,:,:,:) * vs_read(:,:,:,:)

  ! free memory
  deallocate(rho_read,vp_read,vs_read)

  ! gets attenuation arrays from files
  if (ATTENUATION) then
    stop 'Anisotropic GLL model with attenuation is not supported yet.'
  endif

end subroutine 