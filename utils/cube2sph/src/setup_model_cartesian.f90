program setup_model_cartesian
  use mpi
  use meshfem3D_models_par, only: myrank
  implicit none
  integer :: ier
  integer, dimension(:,:,:,:), allocatable :: ibool
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read,x,y,z,vp_new,vs_new,rho_new
  real, dimension(:), allocatable :: xstore, ystore, zstore
  integer, parameter :: MAX_STRING_LEN = 512  ! constants.h
  character(len=MAX_STRING_LEN), parameter :: LOCAL_PATH = 'DATABASES_MPI'
  character(len=MAX_STRING_LEN), parameter :: LOCAL_PATH_REF = 'DATABASES_MPI_REF'
  character(len=MAX_STRING_LEN) :: prname_lp,filename
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX
  integer, parameter :: IIN = 40,IOUT = 41
  integer :: nspec, i, j, k, ispec, nspec_irregular, nglob, iglob
  call MPI_Init(ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  ! user output
  if (myrank == 0) then
    write(*,*) '     reading mesh from: ',trim(LOCAL_PATH_REF)
  endif

  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH_REF)// '/' //'proc',myrank,'_'
  ! reading in number of elements
  if (myrank == 0) write(*,*) '     reading in: external_mesh.bin'
  filename = prname_lp(1:len_trim(prname_lp))//'external_mesh.bin'
  open(unit=IIN,file=trim(filename),status='unknown',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening database proc######_external_mesh.bin'
  read(IIN) nspec
  allocate(ibool(NGLLX, NGLLY, NGLLZ, nspec), &
           x(NGLLX, NGLLY, NGLLZ, nspec), &
           y(NGLLX, NGLLY, NGLLZ, nspec), &
           z(NGLLX, NGLLY, NGLLZ, nspec), stat=ier)
  if (ier/=0) stop 'error allocating xyz array'
  read(IIN) nglob
  allocate(xstore(nglob), ystore(nglob), zstore(nglob), stat=ier)
  if (ier/=0) stop 'error allocating xyz store array'
  read(IIN) nspec_irregular
  read(IIN) ibool
  read(IIN) xstore
  read(IIN) ystore
  read(IIN) zstore
  close(IIN)
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          x(i,j,k,ispec) = xstore(iglob)
          y(i,j,k,ispec) = ystore(iglob)
          z(i,j,k,ispec) = zstore(iglob)
        enddo
      enddo
    enddo
  enddo
  deallocate(xstore, ystore, zstore, ibool)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! if only vp structure is available (as is often the case in exploration seismology),
  !!! use lines for vp only

  ! density
  allocate(rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array rho_read'

    ! user output
  !if (myrank == 0) write(*,*) '     reading in: rho.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'rho.bin'
  if (myrank == 0) print *, 'reading ', filename
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading rho.bin file'
  endif

  read(IIN) rho_read
  close(IIN)

  ! vp
  allocate(vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array vp_read'

  ! user output
  !if (myrank == 0) write(*,*) '     reading in: vp.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vp.bin'
  if (myrank == 0) print *, 'reading ', filename
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vp.bin file'
  endif

  read(IIN) vp_read
  close(IIN)

  ! vs
  allocate(vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array vs_read'

  ! user output
  !if (myrank == 0) write(*,*) '     reading in: vs.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vs.bin'
  if (myrank == 0) print *, 'reading ', filename
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vs.bin file'
  endif

  read(IIN) vs_read
  close(IIN)

  allocate(vp_new(NGLLX,NGLLY,NGLLZ,nspec), &
           vs_new(NGLLX,NGLLY,NGLLZ,nspec), &
           rho_new(NGLLX,NGLLY,NGLLZ,nspec), stat=ier)
  if (ier/=0) then
    stop 'error allocating array for new structural parameters'
  endif
  call update_parameters_ref(x, y, z, rho_read, vp_read, vs_read, &
                         rho_new, vp_new, vs_new, &
                         NGLLX, NGLLY, NGLLZ, nspec)
  !!!!!!!!!!!!!!!
  !rho_read(:,:,:,:) = rho_new(:,:,:,:)
  !vs_read(:,:,:,:) = vs_new(:,:,:,:)
  !vp_read(:,:,:,:) = vp_new(:,:,:,:)
  !!!!!!!!!!!!!!!
  !call update_parameters_from_netcdf(x, y, z, rho_read, vp_read, vs_read, &
  !                                   rho_new, vp_new, vs_new, &
  !                                   NGLLX, NGLLY, NGLLZ, nspec)
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)// '/' //'proc',myrank,'_'
  if (myrank == 0) print *, 'writing to ', prname_lp(1:len_trim(prname_lp))//'vp.bin'
  ! write new vp
  open(unit=IOUT,file=prname_lp(1:len_trim(prname_lp))//'vp.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file vp.bin'
  write(IOUT) vp_new
  close(IOUT)

  ! write new vs
  if (myrank == 0) print *, 'writing to ', prname_lp(1:len_trim(prname_lp))//'vs.bin'
  open(unit=IOUT,file=prname_lp(1:len_trim(prname_lp))//'vs.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file vs.bin'
  write(IOUT) vs_new
  close(IOUT)

  ! write new rho
  if (myrank == 0) print *, 'writing to ', prname_lp(1:len_trim(prname_lp))//'rho.bin'
  open(unit=IOUT,file=prname_lp(1:len_trim(prname_lp))//'rho.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening file rho.bin'
  write(IOUT) rho_new
  close(IOUT)

  deallocate(vp_read, vs_read, rho_read, x, y, z, vp_new, vs_new, rho_new)
  call MPI_Finalize(ier)
end program setup_model_cartesian

subroutine update_parameters_ref(x, y, z, rho_old, vp_old, vs_old, rho_new, &
                     vp_new, vs_new, nx, ny, nz, nspec)
  use constants
  use meshfem3D_par, only: &
    iregion_code, idoubling, &
    RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER, &
    ISOTROPIC_3D_MANTLE, USE_FULL_TISO_MANTLE, LOCAL_PATH, &
    RCMB,RICB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220, &
    R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
    ABSORBING_CONDITIONS

  use meshfem3D_models_par

  use create_regions_mesh_par2, only: &
    Qmu_store,tau_e_store,tau_s,T_c_source
  implicit none
  integer :: ispec,nspec,inode
  logical :: elem_in_crust,elem_in_mantle,found_crust, elem_is_tiso
  logical, dimension(:), allocatable :: ispec_is_tiso
  ! local parameters
  double precision, dimension(NGNOD) :: xelm, yelm, zelm
  double precision :: xmesh,ymesh,zmesh, rmin,rmax
  ! the 21 coefficients for an anisotropic medium in reduced notation
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                   c34,c35,c36,c44,c45,c46,c55,c56,c66

  double precision :: Qkappa,Qmu
  double precision, dimension(N_SLS) :: tau_e

  double precision :: rho,dvp,scale_vel,scale_rho
  double precision :: vpv,vph,vsv,vsh,eta_aniso
  double precision :: lat,lon,r_dummy,theta,phi,vpvc,vphc,vsvc,vshc,etac,rhoc
  double precision :: r,r_prem,moho,rmid
  integer :: i,j,k,i_sls
  integer :: nx,ny,nz,imid,jmid,kmid
  integer, dimension(NGNOD) :: iaddx, iaddy, iaddz
  real, dimension(nx,ny,nz,nspec) :: x,y,z,rho_old,vs_old,vp_old, &
                                        rho_new,vs_new,vp_new
  double precision, dimension(nx,ny,nz,nspec) :: xstore,ystore,zstore

  double precision xigll(nx),yigll(ny),zigll(nz)
  double precision wxgll(nx),wygll(ny),wzgll(nz)

  double precision shape3D(NGNOD,nx,ny,nz),dershape3D(NDIM,NGNOD,nx,ny,nz)
  double precision :: vs_min

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !set parameters
  !myrank = 0
  scale_vel = R_EARTH * dsqrt(PI*GRAV*RHOAV) ! m/s
  scale_rho = RHOAV ! kg/m^3
  imid = 3
  jmid = 3
  kmid = 3
  LOCAL_PATH = 'DATABASES_MPI'
  ELLIPTICITY = .true.
  TOPOGRAPHY = .true.
  CASE_3D = .true.
  CRUSTAL = .true.
  ISOTROPIC_3D_MANTLE = .true.
  ONE_CRUST = .true.
  !REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
  !THREE_D_MODEL = THREE_D_MODEL_GLL
  TRANSVERSE_ISOTROPY = .false.
  !CRUSTAL = .true.
  !!!!!!! change mantle models here
  !REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
  !THREE_D_MODEL = THREE_D_MODEL_S362ANI
  REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
  THREE_D_MODEL = THREE_D_MODEL_S40RTS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RMIDDLE_CRUST = 6356000.d0
  RMOHO = 6346600.d0 ! at 24.4km depth
  R80  = 6291000.d0
  R220 = 6151000.d0
  R400 = 5961000.d0 ! 410km discontinuity
  R600 = 5771000.d0
  R670 = 5721000.d0 ! 650km discontinuity
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  ! 3D models do not honor PREM moho but a fictitious moho at 40km depth:
  ! either to make simulation cheaper or to have a 3D crustal structure
  RMOHO_FICTITIOUS_IN_MESHER = (R80 + R_EARTH) / 2.0d0
  iregion_code = IREGION_CRUST_MANTLE
  R80_FICTITIOUS_IN_MESHER = R80
  vs_min = 1000.0  ! m/s
  !vs_min = 0.0
  if (CRUSTAL .and. CASE_3D) then
    !> Hejun
    ! mesh will honor 3D crustal moho topography
    ! moves MOHO up 5km to honor moho topography deeper than 35 km
    ! moves R80 down to 120km depth in order to have less squeezing for elements below moho
    RMOHO_FICTITIOUS_IN_MESHER = RMOHO_FICTITIOUS_IN_MESHER + RMOHO_STRETCH_ADJUSTMENT
    R80_FICTITIOUS_IN_MESHER = R80_FICTITIOUS_IN_MESHER + R80_STRETCH_ADJUSTMENT
  endif
  allocate(idoubling(nspec), ispec_is_tiso(nspec))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !iregion_code = IREGION_CRUST_MANTLE
  !CRUSTAL = .true.
  !! turn on this line if REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
    call model_1dref_broadcast(CRUSTAL)
  endif
  !THREE_D_MODEL = THREE_D_MODEL_S29EA
  select case (THREE_D_MODEL)
    case (THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
            THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
      call model_s362ani_broadcast(THREE_D_MODEL)
    case (THREE_D_MODEL_S40RTS)
      call model_s40rts_broadcast()
  end select
  call meshfem3D_crust_broadcast()
  call hex_nodes(iaddx, iaddy, iaddz)
  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,nx,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,ny,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,nz,GAUSSALPHA,GAUSSBETA)

  ! get the 3-D shape functions
  call get_shape3D(shape3D,dershape3D,xigll,yigll,zigll)
  ! loops over all GLL points for this spectral element
  do ispec = 1,nspec
    elem_in_crust = .false.
    elem_in_mantle = .false.
    elem_is_tiso = .false.
    do inode = 1, NGNOD
      i = iaddx(inode) * (imid - 1) + 1
      j = iaddy(inode) * (jmid - 1) + 1
      k = iaddz(inode) * (kmid - 1) + 1
      xelm(inode) = x(i,j,k,ispec) / R_EARTH
      yelm(inode) = y(i,j,k,ispec) / R_EARTH
      zelm(inode) = z(i,j,k,ispec) / R_EARTH
    enddo
    rmid = sqrt(xelm(NGNOD) * xelm(NGNOD) + yelm(NGNOD) * yelm(NGNOD) &
                  + zelm(NGNOD) * zelm(NGNOD))
    if (rmid > RMOHO_FICTITIOUS_IN_MESHER / R_EARTH) then
      idoubling(ispec) = IFLAG_CRUST
      rmax = ONE
      rmin = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
    else if (rmid > R80_FICTITIOUS_IN_MESHER / R_EARTH &
            .and. rmid <= RMOHO_FICTITIOUS_IN_MESHER / R_EARTH) then
      idoubling(ispec) = IFLAG_80_MOHO
      rmax = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmin = R80_FICTITIOUS_IN_MESHER / R_EARTH
    else if (rmid > R220 / R_EARTH .and. &
            rmid <= R80_FICTITIOUS_IN_MESHER / R_EARTH) then
      idoubling(ispec) = IFLAG_220_80
      rmax = R80_FICTITIOUS_IN_MESHER / R_EARTH
      rmin = R220 / R_EARTH
    else if (rmid > R670 / R_EARTH .and. rmid <= R220 / R_EARTH) then
      idoubling(ispec) = IFLAG_670_220
      if (rmid > R400 / R_EARTH) then
        rmax = R220 / R_EARTH
        rmin = R400 / R_EARTH
      else if (rmid > R600 / R_EARTH .and. rmid <= R400 / R_EARTH) then
        rmax = R400 / R_EARTH
        rmin = R600 / R_EARTH
      else
        rmax = R600 / R_EARTH
        rmin = R670 / R_EARTH
      end if
    else
      idoubling(ispec) = IFLAG_MANTLE_NORMAL
      if (rmid > R771 / R_EARTH) then
        rmax = R670 / R_EARTH
        rmin = R771 / R_EARTH
      else
        rmax = R771 / R_EARTH
        rmin = RTOPDDOUBLEPRIME / R_EARTH
      end if
    end if
    if (iregion_code == IREGION_CRUST_MANTLE) then
      if (CRUSTAL .and. CASE_3D) then
        ! 3D crustal models
        if (idoubling(ispec) == IFLAG_CRUST &
          .or. idoubling(ispec) == IFLAG_220_80 &
          .or. idoubling(ispec) == IFLAG_80_MOHO) then
          ! Stretch mesh to honor smoothed moho thickness from crust2.0
          ! mesh is stretched between surface and 220
          !
          ! differentiate between regional and global meshing
          if (REGIONAL_MOHO_MESH) then
            call moho_stretching_honor_crust_reg(myrank,xelm,yelm,zelm, &
                                                 elem_in_crust,elem_in_mantle)
          else
            call moho_stretching_honor_crust(myrank,xelm,yelm,zelm, &
                                             elem_in_crust,elem_in_mantle)
          endif
        else
          ! element below 220km
          ! sets element flag for mantle
          elem_in_mantle = .true.
        endif
      else
        ! 1D crust, no stretching
        ! sets element flags
        if (idoubling(ispec) == IFLAG_CRUST) then
          elem_in_crust = .true.
        else
          elem_in_mantle = .true.
        endif
      endif
      ! sets transverse isotropic flag for elements in mantle
      if (TRANSVERSE_ISOTROPY) then
        ! modifies tiso to have it for all mantle elements
        ! preferred for example, when using 1Dref (STW model)
        if (USE_FULL_TISO_MANTLE) then
          ! all elements below the actual moho will be used for transverse isotropy
          ! note: this will increase the computation time by ~ 45 %
          if (elem_in_mantle) then
            elem_is_tiso = .true.
          endif
        else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
          ! transverse isotropic mantle between fictitious moho to 670km depth
          ! preferred for Harvard (Kustowski's) models using STW 1D reference, i.e.
          ! THREE_D_MODEL_S362ANI
          ! THREE_D_MODEL_S362WMANI
          ! THREE_D_MODEL_S29EA
          ! THREE_D_MODEL_GLL
          ! which show significant transverse isotropy also below 220km depth
          if (USE_OLD_VERSION_5_1_5_FORMAT) then
            ! assigns TI only to elements below (2-layer) fictitious moho down to 670
            if (idoubling(ispec) == IFLAG_220_80 &
              .or. idoubling(ispec) == IFLAG_80_MOHO &
              .or. idoubling(ispec) == IFLAG_670_220) then
              elem_is_tiso = .true.
            endif
          else
            ! assigns TI to elements in mantle elements just below actual moho down to 670
            if (idoubling(ispec) == IFLAG_220_80 &
              .or. idoubling(ispec) == IFLAG_80_MOHO &
              .or. idoubling(ispec) == IFLAG_670_220 &
              .or. (idoubling(ispec) == IFLAG_CRUST .and. elem_in_mantle )) then
              elem_is_tiso = .true.
            endif
          endif
        else if (idoubling(ispec) == IFLAG_220_80 .or. idoubling(ispec) == IFLAG_80_MOHO) then
          ! default case for PREM reference models:
          ! models use only transverse isotropy between moho and 220 km depth
          elem_is_tiso = .true.
          ! checks mantle flag to be sure
          !if (elem_in_mantle .eqv. .false. ) stop 'Error mantle flag confused between moho and 220'
        endif
      endif
    endif
    ! sets element tiso flag
    ispec_is_tiso(ispec) = elem_is_tiso

    ! interpolates and stores GLL point locations
    call compute_element_GLL_locations(xelm,yelm,zelm,ispec,nspec, &
                                     xstore,ystore,zstore,shape3D)

    ! computes velocity/density/... values for the chosen Earth model
    ! (only needed for second meshing phase)
    ! it is *CRUCIAL* to leave this initialization here, this was the cause of the "s362ani + attenuation" bug in 2013 and 2014
    ! thus please never remove the line below
    moho = 0.d0
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx

          ! initializes values
          rho = 0.d0
          vpv = 0.d0
          vph = 0.d0
          vsv = 0.d0
          vsh = 0.d0
          eta_aniso = 0.d0
          c11 = 0.d0
          c12 = 0.d0
          c13 = 0.d0
          c14 = 0.d0
          c15 = 0.d0
          c16 = 0.d0
          c22 = 0.d0
          c23 = 0.d0
          c24 = 0.d0
          c25 = 0.d0
          c26 = 0.d0
          c33 = 0.d0
          c34 = 0.d0
          c35 = 0.d0
          c36 = 0.d0
          c44 = 0.d0
          c45 = 0.d0
          c46 = 0.d0
          c55 = 0.d0
          c56 = 0.d0
          c66 = 0.d0
          Qmu = 0.d0
          Qkappa = 0.d0 ! not used, not stored so far...
          tau_e(:) = 0.d0
          dvp = 0.d0

          xmesh = xstore(i,j,k,ispec)
          ymesh = ystore(i,j,k,ispec)
          zmesh = zstore(i,j,k,ispec)

          r = dsqrt(xmesh*xmesh + ymesh*ymesh + zmesh*zmesh)
          r_prem = r
          if (r <= rmin*1.000001d0) r_prem = rmin*1.000001d0
          if (r >= rmax*0.999999d0) r_prem = rmax*0.999999d0
          call meshfem3D_models_get1D_val(iregion_code,idoubling, &
                          r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                          Qkappa,Qmu,RICB,RCMB, &
                          RTOPDDOUBLEPRIME,R80,R120,R220,R400,R600,R670,R771, &
                          RMOHO,RMIDDLE_CRUST,ROCEAN)
          ! gets the 3-D model parameters for the mantle
          call meshfem3D_models_get3Dmntl_val(iregion_code,r_prem,rho,dvp, &
                              vpv,vph,vsv,vsh,eta_aniso, &
                              RCMB,R670,RMOHO, &
                              xmesh,ymesh,zmesh,r, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                              c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
          if (CRUSTAL .and. .not. elem_in_mantle) &
            call meshfem3D_models_get3Dcrust_val(iregion_code,xmesh,ymesh, &
                              zmesh,r, &
                              vpv,vph,vsv,vsh,rho,eta_aniso,dvp, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                              c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                              elem_in_crust,moho)
          rho_new(i,j,k,ispec) = real(rho * scale_rho, kind=CUSTOM_REAL)
          vp_new(i,j,k,ispec) = real(vph * scale_vel, kind=CUSTOM_REAL)
          vs_new(i,j,k,ispec) = real(vsh * scale_vel, kind=CUSTOM_REAL)
          if (vs_new(i,j,k,ispec) < vs_min) then
            !vp_new(i,j,k,ispec) = sqrt(vp_new(i,j,k,ispec) * &
            !vp_new(i,j,k,ispec) + 4.d0 * (vs_min * vs_min &
            !- vs_new(i,j,k,ispec) * vs_new(i,j,k,ispec)) / 3.d0)
            vp_new(i,j,k,ispec) = vp_new(i,j,k,ispec) * &
                             (vs_min / vs_new(i,j,k,ispec))
            vs_new(i,j,k,ispec) = vs_min
            if (elem_in_mantle) print *, 'bad vs: ', myrank
          endif
        enddo
      enddo
    enddo
  enddo

  deallocate(idoubling, ispec_is_tiso)


end subroutine update_parameters_ref

subroutine update_parameters_from_netcdf(x, y, z, rho, vp, vs, rho_new, vp_new, vs_new, &
                nx, ny, nz, nspec)
  implicit none
  integer :: nx, ny, nz, nspec, i, j, k, ispec
  integer, parameter :: SIZE_REAL = 4, SIZE_DOUBLE = 8
  integer, parameter :: CUSTOM_REAL = SIZE_REAL
  real(kind=CUSTOM_REAL),dimension(nx,ny,nz,nspec) :: x, y, z, rho, vp, vs, &
          rho_new, vp_new, vs_new
  !real(kind=CUSTOM_REAL) :: r, theta, phi, lat, lon, depth, xmesh, ymesh, zmesh
  double precision :: r, theta, phi, lat, lon, depth, xmesh, ymesh, zmesh
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: PI_OVER_TWO = PI / 2.0d0
  double precision, parameter :: RADIANS_TO_DEGREES = 180.d0 / PI
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0
  double precision, parameter :: R_EARTH = 6371000.d0
  double precision, parameter :: R_EARTH_KM = R_EARTH / 1000.d0
  double precision, parameter :: LARGE_VALUE = 9990000.d0 ! m/s
  !!!!!! put model metadata here !!!!!!!
  character(len = *), parameter :: model_fn = 'model_Ward18/Alaska.ANT+RF.Ward.2018_kmps.nc'
  double precision, parameter :: depth_start = 1.0, lat_start = 52.0, lon_start = -113.0, &
                            ddepth = 1.0, dlat = 0.1, dlon = -0.1
  integer, parameter :: ndepth = 70, nlat = 211, nlon = 601
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real, dimension(nlon,nlat,ndepth) :: vs_rec
  real, dimension(ndepth) :: depth_grid
  real, dimension(nlat) :: lat_grid
  real, dimension(nlon) :: lon_grid
  real, dimension(2, 2, 2) :: vs_square
  integer :: idepth, ilat, ilon
  double precision :: ratio_depth, ratio_lat, ratio_lon

  call read_netcdf_model(model_fn, ndepth, nlat, nlon, depth_grid, lat_grid, lon_grid, vs_rec)
  vs_rec(:,:,:) = vs_rec(:,:,:) * 1000.d0 ! m/s 
  do ispec = 1,nspec
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          xmesh = x(i,j,k,ispec)
          ymesh = y(i,j,k,ispec)
          zmesh = z(i,j,k,ispec)

          !call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r,theta,phi)
          !call reduce(theta,phi)

          ! lat/lon in degrees (range lat/lon = [-90,90] / [-180,180]
          !lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
          !lon = phi * RADIANS_TO_DEGREES
          call xyz_2_rlatlon_dble(xmesh,ymesh,zmesh,r,lat,lon)
          if (lon > 180.0d0 ) lon = lon - 360.0d0
          !depth = (ONE - r) * R_EARTH_KM
          depth = R_EARTH_KM - r / 1000.0d0
          !!!! interpolation !!!!
          idepth = floor((depth - depth_start) / ddepth) + 1
          ilat = floor((lat - lat_start) / dlat) + 1
          ilon = floor((lon - lon_start) / dlon) + 1
          if ((idepth .ge. ndepth) .or. &
              (ilat .le. 0) .or. (ilat .ge. nlat) .or. &
              (ilon .le. 0) .or. (ilon .ge. nlon)) then  !! out of grid
            vs_new(i,j,k,ispec) = vs(i,j,k,ispec)
          else
            if (idepth .le. 0) then
              vs_square(:,:,1) = vs_rec(ilon:ilon + 1, ilat:ilat + 1, 1)
              vs_square(:,:,2) = vs_rec(ilon:ilon + 1, ilat:ilat + 1, 1)
            else
              vs_square(:,:,:) = vs_rec(ilon:ilon + 1, ilat:ilat + 1, idepth:idepth + 1)
            endif
            if ((vs_square(1,1,1) .gt. LARGE_VALUE) .or. &
                (vs_square(2,1,1) .gt. LARGE_VALUE) .or. &
                (vs_square(1,2,1) .gt. LARGE_VALUE) .or. &
                (vs_square(2,2,1) .gt. LARGE_VALUE) .or. &
                (vs_square(1,1,2) .gt. LARGE_VALUE) .or. &
                (vs_square(2,1,2) .gt. LARGE_VALUE) .or. &
                (vs_square(1,2,2) .gt. LARGE_VALUE) .or. &
                (vs_square(2,2,2) .gt. LARGE_VALUE)) then !! undefined value
              vs_new(i,j,k,ispec) = vs(i,j,k,ispec)
            else
              ratio_depth = (depth - depth_grid(idepth)) / ddepth
              ratio_lat = (lat - lat_grid(ilat)) / dlat
              ratio_lon = (lon - lon_grid(ilon)) / dlon
              vs_new(i,j,k,ispec) = vs_square(1,1,1) * (ONE - ratio_depth) * (ONE - ratio_lat) * (ONE - ratio_lon) + &
                      vs_square(1,1,2) * ratio_depth * (ONE - ratio_lat) * (ONE - ratio_lon) + &
                      vs_square(1,2,1) * (ONE - ratio_depth) * ratio_lat * (ONE - ratio_lon) + &
                      vs_square(1,2,2) * ratio_depth * ratio_lat * (ONE - ratio_lon) + &
                      vs_square(2,1,1) * (ONE - ratio_depth) * (ONE - ratio_lat) * ratio_lon + &
                      vs_square(2,1,2) * ratio_depth * (ONE - ratio_lat) * ratio_lon + &
                      vs_square(2,2,1) * (ONE - ratio_depth) * ratio_lat * ratio_lon + &
                      vs_square(2,2,2) * ratio_depth * ratio_lat * ratio_lon
            end if
          end if
          ! check for bad vs values
          !if (vs_new(i,j,k,ispec) > vp(i,j,k,ispec)) then
          !  print *, 'bad vs value.'
          !  print *, 'new vs: ', vs_new(i,j,k,ispec), ', old vs: ',
          !  vs(i,j,k,ispec), ', vp: ', vp(i,j,k,ispec)
          !  print *, 'lat: ', lat, ', lon: ', lon, ', depth: ', depth
          !endif
          ! set bulk velocity constant
          !vp_new(i,j,k,ispec) = sqrt(vp(i,j,k,ispec) * vp(i,j,k,ispec) + 4.d0 * &
          !(vs_new(i,j,k,ispec) * vs_new(i,j,k,ispec) - vs(i,j,k,ispec) * &
          ! vs(i,j,k,ispec)) / 3.d0)
          ! set the Poisson ratio constant
          vp_new(i,j,k,ispec) = vp(i,j,k,ispec) * (vs_new(i,j,k,ispec) &
                                / vs(i,j,k,ispec))
        enddo
      enddo
    enddo
  enddo
  rho_new(:,:,:,:) = rho(:,:,:,:)
  !vp_new(:,:,:,:) = vp(:,:,:,:)
  !vs_new = vs
end subroutine update_parameters_from_netcdf

