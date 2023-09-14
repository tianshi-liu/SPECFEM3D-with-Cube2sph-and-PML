subroutine setup_gll_model_cartesian(xstore,ystore,zstore,rho_new,&
                  vp_new,vs_new,nspec)
  use meshfem3D_models_par, only: myrank
  use constants, only: NGLLX,NGLLY,NGLLZ,IIN,IOUT,IMAIN
  implicit none
  integer :: ier,nspec
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
  real, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: vp_read,vs_read,rho_read,&
          vp_new,vs_new,rho_new
  integer, parameter :: MAX_STRING_LEN = 512  ! constants.h
  character(len=MAX_STRING_LEN), parameter :: LOCAL_PATH_REF = &
          'DATABASES_MPI_REF'
  character(len=MAX_STRING_LEN) :: prname_lp,filename
  !integer, parameter :: NGLLX = 5
  !integer, parameter :: NGLLY = NGLLX
  !integer, parameter :: NGLLZ = NGLLX
  !integer, parameter :: IIN = 40,IOUT = 41
  integer :: i, j, k, ispec, nspec_irregular


  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH_REF)// '/' //'proc',myrank,'_'
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     reading mesh from: ',trim(LOCAL_PATH_REF)
  endif
  ! density
  !allocate(rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  !if (ier /= 0) stop 'error allocating array rho_read'

    ! user output
  !if (myrank == 0) write(*,*) '     reading in: rho.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'rho.bin'
  if (myrank == 0) write(IMAIN,*) 'reading ', filename
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading rho.bin file'
  endif

  read(IIN) rho_read
  close(IIN)

  ! vp
  !allocate(vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  !if (ier /= 0) stop 'error allocating array vp_read'

  ! user output
  !if (myrank == 0) write(*,*) '     reading in: vp.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vp.bin'
  if (myrank == 0) write(IMAIN,*) 'reading ', filename
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vp.bin file'
  endif

  read(IIN) vp_read
  close(IIN)

  filename = prname_lp(1:len_trim(prname_lp))//'vs.bin'
  if (myrank == 0) write(IMAIN,*) 'reading ', filename
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vs.bin file'
  endif

  read(IIN) vs_read
  close(IIN)

  if (ier/=0) then
    stop 'error allocating array for new structural parameters'
  endif
  call update_structural_parameters(xstore,ystore,zstore,&
                         rho_read,vp_read,vs_read, &
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
  !deallocate(rho_read,vp_read,vs_read)
end subroutine setup_gll_model_cartesian

subroutine update_structural_parameters(x, y, z, rho_old, vp_old, vs_old, &
                     rho_new, vp_new, vs_new, nx, ny, nz, nspec)
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
  real, dimension(nx,ny,nz,nspec) :: rho_old,vs_old,vp_old, &
                                        rho_new,vs_new,vp_new
  double precision, dimension(nx,ny,nz,nspec) :: x,y,z,xstore,ystore,zstore

  double precision xigll(nx),yigll(ny),zigll(nz)
  double precision wxgll(nx),wygll(ny),wzgll(nz)

  double precision shape3D(NGNOD,nx,ny,nz),dershape3D(NDIM,NGNOD,nx,ny,nz)
  double precision :: vs_min
  character(len=256) :: dummy_string, model_string, emc_path
  logical :: use_emc_model = .false.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !set parameters
  !myrank = 0
  scale_vel = R_EARTH * dsqrt(PI*GRAV*RHOAV) ! m/s
  scale_rho = RHOAV ! kg/m^3
  imid = 3
  jmid = 3
  kmid = 3
  LOCAL_PATH = 'DATABASES_MPI'

  open(unit=95, file='Cube2sph_model_par', action='read', form='formatted')
  call find_parameter_value('ELLIPTICITY', dummy_string, 95)
  read(dummy_string, *) ELLIPTICITY
  call find_parameter_value('TOPOGRAPHY', dummy_string, 95)
  read(dummy_string, *) TOPOGRAPHY
  call find_parameter_value('CASE_3D', dummy_string, 95)
  read(dummy_string, *) CASE_3D
  call find_parameter_value('CRUSTAL', dummy_string, 95)
  read(dummy_string, *) CRUSTAL
  call find_parameter_value('ISOTROPIC_3D_MANTLE', dummy_string, 95)
  read(dummy_string, *) ISOTROPIC_3D_MANTLE
  call find_parameter_value('ONE_CRUST', dummy_string, 95)
  read(dummy_string, *) ONE_CRUST
  call find_parameter_value('TRANSVERSE_ISOTROPY', dummy_string, 95)
  read(dummy_string, *) TRANSVERSE_ISOTROPY
  call find_parameter_value('MODEL', dummy_string, 95)
  model_string = trim(dummy_string)

  !ELLIPTICITY = .false.
  !TOPOGRAPHY = .false.
  !CASE_3D = .false.
  !CRUSTAL = .false.
  !ISOTROPIC_3D_MANTLE = .false.
  !ONE_CRUST = .true.
  !REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
  !THREE_D_MODEL = THREE_D_MODEL_GLL
  !TRANSVERSE_ISOTROPY = .false.
  !CRUSTAL = .true.
  !!!!!!! change mantle models here
  !REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
  !THREE_D_MODEL = THREE_D_MODEL_S362ANI
  !REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
  !THREE_D_MODEL = THREE_D_MODEL_S40RTS
  
  select case (trim(model_string))
    case ('default')
      rho_new(:,:,:,:) = rho_old(:,:,:,:)
      vs_new(:,:,:,:) = vs_old(:,:,:,:)
      vp_new(:,:,:,:) = vp_old(:,:,:,:)
      write(IMAIN,*) 'use default model'
      return
    case ('PREM')
      CASE_3D = .false.
      ISOTROPIC_3D_MANTLE = .false.
      REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    case ('IASP91')
      CASE_3D = .false.
      ISOTROPIC_3D_MANTLE = .false.
      REFERENCE_1D_MODEL = REFERENCE_MODEL_IASP91
    case ('1DREF')
      CASE_3D = .false.
      ISOTROPIC_3D_MANTLE = .false.
      REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    case ('S40RTS')
      CASE_3D = .true.
      ISOTROPIC_3D_MANTLE = .true.
      REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
      THREE_D_MODEL = THREE_D_MODEL_S40RTS
    case default
      if (index(trim(model_string), 'EMC:') /= 1) stop 'incorrect model'
      CASE_3D = .true.
      ISOTROPIC_3D_MANTLE = .true.
      REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
      THREE_D_MODEL = THREE_D_MODEL_S40RTS
      use_emc_model = .true.
      emc_path = model_string(5:len_trim(model_string))
  end select
  close(95)
  if (myrank==0) then
    write(IMAIN,*) 'implementing Cube2sph structrual model'
    write(IMAIN,*) 'ELLIPTICITY = ', ELLIPTICITY
    write(IMAIN,*) 'TOPOGRAPHY = ', TOPOGRAPHY
    write(IMAIN,*) 'CASE_3D = ', CASE_3D
    write(IMAIN,*) 'CRUSTAL = ', CRUSTAL
    write(IMAIN,*) 'ISOTROPIC_3D_MANTLE = ', ISOTROPIC_3D_MANTLE
    write(IMAIN,*) 'ONE_CRUST = ', ONE_CRUST
    write(IMAIN,*) 'TRANSVERSE_ISOTROPY = ', TRANSVERSE_ISOTROPY
    write(IMAIN,*) 'MODEL = ', trim(model_string)
    if (use_emc_model) then
      write(IMAIN,*) 'use IRIS EMC model at ', emc_path
    endif
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! default: PREM
  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0 ! at 15km depth
  RMOHO = 6346600.d0  ! at 24.4km depth
  R80  = 6291000.d0
  R120 = -1.d0   ! by default there is no d120 discontinuity, except in IASP91, therefore set to fictitious value
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0

  ! differing 1-D model radii
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
    ! IASP91
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6351000.d0
    RMOHO = 6336000.d0 ! ! at 35km depth
    R80  = 6291000.d0
    R120 = 6251000.d0
    R220 = 6161000.d0
    R400 = 5961000.d0
    ! there is no d600 discontinuity in IASP91 therefore this value is useless
    ! but it needs to be there for compatibility with other subroutines
    R600 = R_EARTH - 600000.d0
    R670 = 5711000.d0
    R771 = 5611000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB = 3482000.d0
    RICB = 1217000.d0
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
    ! REF
    ROCEAN = 6368000.d0
    RMIDDLE_CRUST = 6356000.d0
    RMOHO = 6346600.d0 ! at 24.4km depth
    R80  = 6291000.d0
    R220 = 6151000.d0
    R400 = 5961000.d0 ! 410km discontinuity
    R600 = 5771000.d0
    R670 = 5721000.d0 ! 650km discontinuity
    R771 = 5600000.d0
    RTOPDDOUBLEPRIME = 3630000.d0
    RCMB = 3479958.d0
    RICB = 1221491.d0
  endif
  
  ! 3D models do not honor PREM moho but a fictitious moho at 40km depth:
  ! either to make simulation cheaper or to have a 3D crustal structure
  RMOHO_FICTITIOUS_IN_MESHER = RMOHO
  if (CASE_3D) RMOHO_FICTITIOUS_IN_MESHER = (R80 + R_EARTH) / 2.0d0
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
  if (CASE_3D) then
    select case (THREE_D_MODEL)
      case (THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
            THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
        call model_s362ani_broadcast(THREE_D_MODEL)
      case (THREE_D_MODEL_S40RTS)
        call model_s40rts_broadcast()
    end select
  endif
  if (CRUSTAL .and. CASE_3D) call meshfem3D_crust_broadcast()
  call hex_nodes(iaddx, iaddy, iaddz)
  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,nx,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,ny,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,nz,GAUSSALPHA,GAUSSBETA)

  ! get the 3-D shape functions
  call get_shape3D(shape3D,dershape3D,xigll,yigll,zigll)
  ! loops over all GLL points for this spectral element
  !print *, shape(xelm), shape(x)
  do ispec = 1,nspec
    elem_in_crust = .false.
    elem_in_mantle = .false.
    elem_is_tiso = .false.
    do inode = 1, NGNOD
      i = iaddx(inode) * (imid - 1) + 1
      j = iaddy(inode) * (jmid - 1) + 1
      k = iaddz(inode) * (kmid - 1) + 1
      !if ((i.lt.1).or.(i.gt.nx).or.(j.lt.1).or.(j.gt.ny).or.&
      !      (k.lt.1).or.(k.gt.nz)) print *, nx,ny,nz
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
          call meshfem3D_models_get1D_val(iregion_code,idoubling(ispec), &
                          r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                          Qkappa,Qmu,RICB,RCMB, &
                          RTOPDDOUBLEPRIME,R80,R120,R220,R400,R600,R670,R771, &
                          RMOHO,RMIDDLE_CRUST,ROCEAN)
          ! gets the 3-D model parameters for the mantle
          if (CASE_3D) &
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
  if (use_emc_model) then
    rho_old(:,:,:,:) = rho_new(:,:,:,:)
    vs_old(:,:,:,:) = vs_new(:,:,:,:)
    vp_old(:,:,:,:) = vp_new(:,:,:,:)
    !!!!!!!!!!!!!!
    call update_parameters_from_netcdf(emc_path, &
                                     xstore, ystore, zstore, &
                                     rho_old, vp_old, vs_old, &
                                     rho_new, vp_new, vs_new, &
                                     NGLLX, NGLLY, NGLLZ, nspec)
  endif
  if (CRUSTAL .and. CASE_3D) call meshfem3D_crust_deallocate()
  deallocate(idoubling, ispec_is_tiso)


end subroutine update_structural_parameters


