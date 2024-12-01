subroutine node_stretching(nodes_coords,nodes_coords_new)
  use constants
  use meshfem3D_par, only: iregion_code, idoubling, MIN_ATTENUATION_PERIOD, &
         MAX_ATTENUATION_PERIOD,ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R120,R220,R400, &
                      R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
                      R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
                      RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER, &
                      ISOTROPIC_3D_MANTLE, USE_FULL_TISO_MANTLE, &
                      LOCAL_PATH
  use meshfem3D_models_par, only: myrank,CASE_3D, &
           REFERENCE_1D_MODEL, THREE_D_MODEL, ONE_CRUST, TRANSVERSE_ISOTROPY,&
           TOPOGRAPHY,ELLIPTICITY,CRUSTAL,ibathy_topo,nspl,rspl,espl,espl2
  use generate_databases_par, only: npts, SAVE_MOHO_MESH, nelmnts_ext_mesh, &
          elmnts_ext_mesh
  implicit none
  !! local variables
  integer :: ipts, ispec, ipass, ipts_dummy, inode, ier
  double precision, dimension(NGNOD) :: xelm, yelm, zelm
  double precision, dimension(NDIM,npts) :: nodes_coords,nodes_coords_new
  logical :: elem_in_crust, elem_in_mantle, elem_is_tiso
  logical, dimension(:), allocatable :: ispec_is_tiso
  double precision :: rmid
  character(len=256) :: dummy_string, model_string, emc_path
  logical :: use_emc_model = .false.
  SAVE_MOHO_MESH = .false.

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
  !THREE_D_MODEL = THREE_D_MODEL_S362ANI
  !TRANSVERSE_ISOTROPY = .false.
  !CRUSTAL = .true.
  !REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
  !THREE_D_MODEL = THREE_D_MODEL_S40RTS
  select case (trim(model_string))
    case ('default')
      CASE_3D = .false.
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
    print *, 'implementing Cube2sph node stretching'
    print *, 'ELLIPTICITY = ', ELLIPTICITY
    print *, 'TOPOGRAPHY = ', TOPOGRAPHY
    print *, 'CASE_3D = ', CASE_3D
    print *, 'CRUSTAL = ', CRUSTAL
    print *, 'ISOTROPIC_3D_MANTLE = ', ISOTROPIC_3D_MANTLE
    print *, 'ONE_CRUST = ', ONE_CRUST
    print *, 'TRANSVERSE_ISOTROPY = ', TRANSVERSE_ISOTROPY
    print *, 'MODEL = ', trim(model_string)
    if (use_emc_model) then
      print *, 'use IRIS EMC model at ', emc_path
    endif
  endif

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
  RMOHO_FICTITIOUS_IN_MESHER = (R80 + R_EARTH) / 2.0d0
  iregion_code = IREGION_CRUST_MANTLE
  R80_FICTITIOUS_IN_MESHER = R80
  if (CRUSTAL .and. CASE_3D) then
     !> Hejun
     ! mesh will honor 3D crustal moho topography
     ! moves MOHO up 5km to honor moho topography deeper than 35 km
     ! moves R80 down to 120km depth in order to have less squeezing for elements below moho
    RMOHO_FICTITIOUS_IN_MESHER = RMOHO_FICTITIOUS_IN_MESHER + RMOHO_STRETCH_ADJUSTMENT
    R80_FICTITIOUS_IN_MESHER = R80_FICTITIOUS_IN_MESHER + R80_STRETCH_ADJUSTMENT
  endif
  allocate(idoubling(nelmnts_ext_mesh), ispec_is_tiso(nelmnts_ext_mesh))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! sets up spline coefficients for ellipticity
  if (ELLIPTICITY) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

  if (TOPOGRAPHY) then
    ! arrays for elevations
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    if (ier /= 0) stop 'Error allocating ibathy_topo array'

    ! initializes
    ibathy_topo(:,:) = 0

    ! sets up topo/bathy
    call model_topo_bathy_broadcast(ibathy_topo,LOCAL_PATH)
  endif
  !call model_1dref_broadcast(CRUSTAL)
  !call model_s362ani_broadcast(THREE_D_MODEL)
  ! crustal model
  if (CRUSTAL) &
    call meshfem3D_crust_broadcast()
  do ispec = 1, nelmnts_ext_mesh
    do inode = 1, NGNOD
      xelm(inode) = nodes_coords(1,elmnts_ext_mesh(inode, ispec)) / R_EARTH
      yelm(inode) = nodes_coords(2,elmnts_ext_mesh(inode, ispec)) / R_EARTH
      zelm(inode) = nodes_coords(3,elmnts_ext_mesh(inode, ispec)) / R_EARTH
    enddo
    rmid = sqrt(xelm(NGNOD) * xelm(NGNOD) + yelm(NGNOD) * yelm(NGNOD) &
                 + zelm(NGNOD) * zelm(NGNOD))
    if (rmid > RMOHO / R_EARTH) then
      idoubling(ispec) = IFLAG_CRUST
    else if (rmid > R80 / R_EARTH .and. rmid <= RMOHO / R_EARTH) then
      idoubling(ispec) = IFLAG_80_MOHO
    else if (rmid > R220 / R_EARTH .and. rmid <= R80 / R_EARTH) then
      idoubling(ispec) = IFLAG_220_80
    else if (rmid > R670 / R_EARTH .and. rmid <= R220 / R_EARTH) then
      idoubling(ispec) = IFLAG_670_220
    else
      idoubling(ispec) = IFLAG_MANTLE_NORMAL
    end if
    elem_in_crust = .false.
    elem_in_mantle = .false.
    elem_is_tiso = .false.
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

    endif ! IREGION_CRUST_MANTLE

    ! sets element tiso flag
    ispec_is_tiso(ispec) = elem_is_tiso
    ! adds surface topography
    if (TOPOGRAPHY) then
      if (idoubling(ispec) == IFLAG_CRUST .or. &
          idoubling(ispec) == IFLAG_220_80 .or. &
          idoubling(ispec) == IFLAG_80_MOHO) then
        ! stretches mesh between surface and R220 accordingly
        ! stretches anchor points only, interpolates GLL points later on
        call add_topography(xelm,yelm,zelm,ibathy_topo)
      endif
    endif
    ! make the Earth elliptical
    if (ELLIPTICITY) then
      ! note: after adding ellipticity, the mesh becomes elliptical and geocentric and geodetic/geographic colatitudes differ.
      ! make the Earth's ellipticity, use element anchor points
      call get_ellipticity(xelm,yelm,zelm,nspl,rspl,espl,espl2)
    endif

    do inode = 1, NGNOD
      nodes_coords_new(1,elmnts_ext_mesh(inode, ispec)) &
              = xelm(inode) * R_EARTH
      nodes_coords_new(2,elmnts_ext_mesh(inode, ispec)) &
              = yelm(inode) * R_EARTH
      nodes_coords_new(3,elmnts_ext_mesh(inode, ispec)) &
              = zelm(inode) * R_EARTH
    enddo
  enddo
end subroutine node_stretching
