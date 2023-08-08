!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
!! TL: anisotropic GLL model
!
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!
!--------------------------------------------------------------------------------------------------

  subroutine model_gll_aniso(myrank,nspec,LOCAL_PATH)

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,&
    MAX_STRING_LEN,ATTENUATION,CUSTOM_REAL,ibool

  use create_regions_mesh_ext_par, only: rhostore,kappastore,mustore,rho_vp,rho_vs,&
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store, &
            xstore_dummy, ystore_dummy, zstore_dummy
  use create_regions_mesh_ext_par, only: AZIMUTHAL_ANISOTROPY, CUBE2SPH_MESH, &
               kappavstore,muvstore,eta_anistore,Gc_nondimstore,Gs_nondimstore

  implicit none

  integer, intent(in) :: myrank,nspec
  logical, parameter :: ISOTROPIC_BULK=.true. ! bulk velocity is isotropic
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:), allocatable :: vph_read,vpv_read,&
                                           vsh_read,vsv_read,rho_read,&
                                           vbulk_read,&
                                           Gc_nondim_read,Gs_nondim_read
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26,&
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision :: rho,vpv,vph,vsv,vsh,eta_aniso,Gc_nondim,Gs_nondim
  double precision :: aa,cc,nn,ll,ff
  double precision :: A,C,F,AL,AN,Gc,Gs
  real(kind=CUSTOM_REAL) :: xp, yp, zp, r_dummy, theta, phi
  
  integer :: ier
  integer :: i,j,k,ispec, iglob
  character(len=MAX_STRING_LEN) :: prname_lp,filename

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

  if (.not. ISOTROPIC_BULK) then
    ! vpv
    allocate(vpv_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
    if (ier /= 0) stop 'error allocating array vpv_read'

    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: vpv.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'vpv.bin'
    open(unit=28,file=trim(filename),status='old',action='read',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(filename)
      stop 'error reading vpv.bin file'
    endif

    read(28) vpv_read
    close(28)

    ! vph
    allocate(vph_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
    if (ier /= 0) stop 'error allocating array vph_read'

    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: vph.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'vph.bin'
    open(unit=28,file=trim(filename),status='old',action='read',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(filename)
      stop 'error reading vph.bin file'
    endif

    read(28) vph_read
    close(28)
  else ! bulk velocity is isotropic
    !vbulk
    allocate(vbulk_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
    if (ier /= 0) stop 'error allocating array vbulk_read'

    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: vbulk.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'vbulk.bin'
    open(unit=28,file=trim(filename),status='old',action='read',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(filename)
      stop 'error reading vbulk.bin file'
    endif

    read(28) vbulk_read
    close(28)
  endif

  ! vsv
  allocate(vsv_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 649')
  if (ier /= 0) stop 'error allocating array vsv_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: vsv.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vsv.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vsv.bin file'
  endif

  read(28) vsv_read
  close(28)

  ! vsh
  allocate(vsh_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 649')
  if (ier /= 0) stop 'error allocating array vsh_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: vsh.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vsh.bin'
  open(unit=28,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vsh.bin file'
  endif

  read(28) vsh_read
  close(28)

  if (AZIMUTHAL_ANISOTROPY) then
    ! Non-dimensionalized Gc
    allocate(Gc_nondim_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
    if (ier /= 0) stop 'error allocating array Gc_nondim_read'

    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: Gc_nondim.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'Gc_nondim.bin'
    open(unit=28,file=trim(filename),status='old',action='read',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(filename)
      stop 'error reading Gc_nondim.bin file'
    endif

    read(28) Gc_nondim_read
    close(28)

    ! Non-dimonsionalized Gs
    allocate(Gs_nondim_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
    if (ier /= 0) stop 'error allocating array Gs_nondim_read'

    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: Gs_nondim.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'Gs_nondim.bin'
    open(unit=28,file=trim(filename),status='old',action='read',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(filename)
      stop 'error reading Gs_nondim.bin file'
    endif

    read(28) Gs_nondim_read
    close(28)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where density structure is not given
  !!! modify according to your desire

  !  rho_read = 1000.0
  !  where ( mustore > 100.0 )  &
  !           rho_read = (1.6612 * (vp_read / 1000.0)     &
  !                      -0.4720 * (vp_read / 1000.0)**2  &
  !                      +0.0671 * (vp_read / 1000.0)**3  &
  !                      -0.0043 * (vp_read / 1000.0)**4  &
  !                      +0.000106*(vp_read / 1000.0)**5)*1000.0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where shear wavespeed structure is not given
  !!! modify according to your desire

  !   vs_read = 0.0
  !   where ( mustore > 100.0 )       vs_read = vp_read / sqrt(3.0)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! update arrays that will be saved and used in the solver xspecfem3D
  !!! the following part is neccessary if you uncommented something above

  rhostore(:,:,:,:) = rho_read(:,:,:,:)
  !! TL: set rho_vp, rho_vs according to vph and vsh, 
  !! same as SPECFEM3D_GLOBE
  if (.not. ISOTROPIC_BULK) then
    kappastore(:,:,:,:) = rhostore(:,:,:,:) * (vph_read(:,:,:,:) * vph_read(:,:,:,:) &
                        - FOUR_THIRDS * vsh_read(:,:,:,:) * vsh_read(:,:,:,:))
    kappavstore(:,:,:,:) = rhostore(:,:,:,:) * (vpv_read(:,:,:,:) * vpv_read(:,:,:,:) &
                        - FOUR_THIRDS * vsv_read(:,:,:,:) * vsv_read(:,:,:,:))
  else
    kappastore(:,:,:,:) = rhostore(:,:,:,:) * vbulk_read(:,:,:,:) * vbulk_read(:,:,:,:)
    kappavstore(:,:,:,:) = rhostore(:,:,:,:) * vbulk_read(:,:,:,:) * vbulk_read(:,:,:,:)
    allocate(vpv_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    allocate(vph_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    vpv_read(:,:,:,:) = sqrt(vbulk_read(:,:,:,:) * vbulk_read(:,:,:,:) + &
                  FOUR_THIRDS * vsv_read(:,:,:,:) * vsv_read(:,:,:,:))
    vph_read(:,:,:,:) = sqrt(vbulk_read(:,:,:,:) * vbulk_read(:,:,:,:) + &
                  FOUR_THIRDS * vsh_read(:,:,:,:) * vsh_read(:,:,:,:))
  endif
  mustore(:,:,:,:) = rhostore(:,:,:,:) * vsh_read(:,:,:,:) * vsh_read(:,:,:,:)

  muvstore(:,:,:,:) = rhostore(:,:,:,:) * vsv_read(:,:,:,:) * vsv_read(:,:,:,:)
  eta_anistore(:,:,:,:) = 1.0
  if (AZIMUTHAL_ANISOTROPY) then
    Gc_nondimstore(:,:,:,:) = Gc_nondim_read(:,:,:,:)
    Gs_nondimstore(:,:,:,:) = Gs_nondim_read(:,:,:,:)
  endif
  rho_vp(:,:,:,:) = rhostore(:,:,:,:) * vph_read(:,:,:,:)
  rho_vs(:,:,:,:) = rhostore(:,:,:,:) * vsh_read(:,:,:,:)
  eta_aniso = 1.0
  do ispec = 1, nspec
    do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX
      vph = dble(vph_read(i,j,k,ispec))
      vpv = dble(vpv_read(i,j,k,ispec))
      vsh = dble(vsh_read(i,j,k,ispec))
      vsv = dble(vsv_read(i,j,k,ispec))
      rho = dble(rho_read(i,j,k,ispec))
      if (AZIMUTHAL_ANISOTROPY) then
        Gc_nondim = dble(Gc_nondim_read(i,j,k,ispec))
        Gs_nondim = dble(Gs_nondim_read(i,j,k,ispec))
      endif
  
      aa = rho*vph*vph
      cc = rho*vpv*vpv
      nn = rho*vsh*vsh
      ll = rho*vsv*vsv
      ff = eta_aniso*(aa - 2.*ll)
      if (AZIMUTHAL_ANISOTROPY) then
        Gc = Gc_nondim * rho * (2.0*vsv*vsv+vsh*vsh) / 3.0
        Gs = Gs_nondim * rho * (2.0*vsv*vsv+vsh*vsh) / 3.0
      endif
  
      A = aa
      C = cc
      AN = nn
      AL = ll
      F = ff

! The mapping from the elastic coefficients to the elastic tensor elements
! in the local Cartesian coordinate system (classical geographic) used in the
! global code (1---South, 2---East, 3---up)
! Always keep the following part when you modify this subroutine
      d11 = A
      d12 = A - 2.0 * AN
      d13 = F
      d14 = 0.0
      d15 = 0.0
      d16 = 0.0
      d22 = A
      d23 = F
      d24 = 0.0
      d25 = 0.0
      d26 = 0.0
      d33 = C
      d34 = 0.0
      d35 = 0.0
      d36 = 0.0
      d44 = AL
      d45 = 0.0
      d46 = 0.0
      d55 = AL
      d56 = 0.0
      d66 = AN
      if (AZIMUTHAL_ANISOTROPY) then
        d44 = d44 - Gc
        d45 = d45 - Gs
        d55 = d55 + Gc
      endif

      iglob = ibool(i,j,k,ispec)
      xp = xstore_dummy(iglob)
      yp = ystore_dummy(iglob)
      zp = zstore_dummy(iglob)
      if (CUBE2SPH_MESH) then
        call xyz_2_rthetaphi(xp,yp,zp,r_dummy,theta,phi)

        call rotate_aniso_tensor(dble(theta),dble(phi),d11,d12,d13,d14,d15,d16, &
                           d22,d23,d24,d25,d26, &
                           d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                           c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
      else
! The mapping to the global Cartesian coordinate system used in the code
! (1---East, 2---North, 3---up)
        c11 = d22
        c12 = d12
        c13 = d23
        c14 = - d25
        c15 = d24
        c16 = - d26
        c22 = d11
        c23 = d13
        c24 = - d15
        c25 = d14
        c26 = - d16
        c33 = d33
        c34 = - d35
        c35 = d34
        c36 = - d36
        c44 = d55
        c45 = - d45
        c46 = d56
        c55 = d44
        c56 = - d46
        c66 = d66        
      endif
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

    enddo;enddo;enddo
  enddo  

  ! gets attenuation arrays from files
  if (ATTENUATION) then
    stop 'Anisotropic GLL model with attenuation is not supported yet.'
  endif

  ! free memory
  deallocate(rho_read,vpv_read,vph_read,vsv_read,vsh_read)
  if (ISOTROPIC_BULK) deallocate(vbulk_read)
  if (AZIMUTHAL_ANISOTROPY) deallocate(Gc_nondim_read,Gs_nondim_read)
  end subroutine model_gll_aniso

