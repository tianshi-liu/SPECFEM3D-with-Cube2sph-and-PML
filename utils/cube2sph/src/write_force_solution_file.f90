program write_force_solution_file
  !! writing FORCESOLUTIONS file in Cartesian coordinate
  !! Tianshi Liu, 2022.2.21
  !! Usage: write_force_solution_file OLD_FILE NEW_FILE TOPOGRAPHY ELLIPTICITY
  !! E.g., write_force_solution_file DATA/FORCESOLUTIONS_geo DATA/FORCESOLUTIONS .true. .true.
  !! takes the DATA/FORCESOLUTIONS_geo in geographic coordinates, writes
  !! DATA/FORCESOLUTIONS in Cartesian coordinate, with the Etopo topography and the ellipticity
  use mpi
  use constants_solver
  use specfem_par, only: &
      ONE_CRUST, &
      tshift_src,theta_source,phi_source, &
      DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz,Mw,M0, &
      rspl,espl,espl2,nspl,ibathy_topo, &
      LOCAL_TMP_PATH,SIMULATION_TYPE,TOPOGRAPHY, &
      nu_source, &
      USE_FORCE_POINT_SOURCE,force_stf,factor_force_source, &
      comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS
  implicit none
  integer, parameter :: NSOURCES = 1
  logical :: ELLIPTICITY
  double precision :: ell
  double precision :: elevation
  double precision :: r0,dcost,p20
  double precision :: theta,phi
  double precision :: x_target_source,y_target_source,z_target_source
  double precision :: r_target_source
  double precision, dimension(NSOURCES) :: lat,long,depth
  double precision :: st,ct,sp,cp
  integer :: iorientation, isource
  integer :: ier
  double precision :: stazi,stdip,thetan,phin,n(3)
  character(len=MAX_STRING_LEN) :: cart_force_fn, sph_force_fn, arg_str
  double precision :: scaleF, fx, fy, fz
  double precision :: min_tshift_src_original
  double precision :: radius
  call MPI_Init(ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  USE_FORCE_POINT_SOURCE = .true.
  !NSOURCES = 1
  ONE_CRUST = .true.
  DT = 0.0
  NUMBER_OF_SIMULTANEOUS_RUNS = 1
  isource = 1
  call get_command_argument(1, sph_force_fn)
  call get_command_argument(2, cart_force_fn)
  call get_command_argument(3, arg_str)
  read(arg_str, *) TOPOGRAPHY
  call get_command_argument(4, arg_str)
  read(arg_str, *) ELLIPTICITY
  !ELLIPTICITY = .true.
  !TOPOGRAPHY = .true.
  allocate(tshift_src(NSOURCES), &
           hdur(NSOURCES), &
           theta_source(NSOURCES), &
           phi_source(NSOURCES), &
           nu_source(NDIM,NDIM,NSOURCES), &
           force_stf(NSOURCES),factor_force_source(NSOURCES), &
           comp_dir_vect_source_E(NSOURCES), &
           comp_dir_vect_source_N(NSOURCES), &
           comp_dir_vect_source_Z_UP(NSOURCES),stat=ier)
  if (TOPOGRAPHY) then
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
    ibathy_topo(:,:) = 0
    call read_topo_bathy_file(ibathy_topo)
  endif
  call get_force(tshift_src,hdur,lat,long,depth,DT,NSOURCES, &
                 min_tshift_src_original,force_stf,factor_force_source, &
                 comp_dir_vect_source_E,comp_dir_vect_source_N, &
                 comp_dir_vect_source_Z_UP,sph_force_fn)

  call lat_2_geocentric_colat_dble(lat(isource),theta)

  phi = long(isource)*DEGREES_TO_RADIANS
  call reduce(theta,phi)
  st=dsin(theta)
  ct=dcos(theta)
  sp=dsin(phi)
  cp=dcos(phi)
  ! record three components for each station
  do iorientation = 1,3

    !   North
    if (iorientation == 1) then
      stazi = 0.d0
      stdip = 0.d0
    !   East
    else if (iorientation == 2) then
      stazi = 90.d0
      stdip = 0.d0
    !   Vertical
    else if (iorientation == 3) then
      stazi = 0.d0
      stdip = - 90.d0
    else
      call exit_MPI(myrank,'incorrect orientation')
    endif

    !   get the orientation of the seismometer
    thetan=(90.0d0+stdip)*DEGREES_TO_RADIANS
    phin=stazi*DEGREES_TO_RADIANS

    ! we use the same convention as in Harvard normal modes for the orientation

    !   vertical component
    n(1) = dcos(thetan)
    !   N-S component
    n(2) = - dsin(thetan)*dcos(phin)
    !   E-W component
    n(3) = dsin(thetan)*dsin(phin)

    !   get the Cartesian components of n in the model: nu
    nu_source(iorientation,1,isource) = n(1)*st*cp+n(2)*ct*cp-n(3)*sp
    nu_source(iorientation,2,isource) = n(1)*st*sp+n(2)*ct*sp+n(3)*cp
    nu_source(iorientation,3,isource) = n(1)*ct-n(2)*st

  enddo
  ! normalized source radius
  r0 = R_UNIT_SPHERE

  ! finds elevation of position
  if (TOPOGRAPHY) then
    call get_topo_bathy(lat(isource),long(isource),elevation,ibathy_topo)
    r0 = r0 + elevation/R_EARTH
  endif
  if (ELLIPTICITY) then
    dcost = dcos(theta)
! this is the Legendre polynomial of degree two, P2(cos(theta)), see the discussion above eq (14.4) in Dahlen and Tromp (1998)
    p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
    radius = r0 - depth(isource)*1000.0d0/R_EARTH
! get ellipticity using spline evaluation
    call spline_evaluation(rspl,espl,espl2,nspl,radius,ell)
! this is eq (14.4) in Dahlen and Tromp (1998)
    r0 = r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
  endif

  ! subtracts source depth (given in km)
  r_target_source = r0 - depth(isource)*1000.0d0/R_EARTH

  ! compute the Cartesian position of the source
  x_target_source = r_target_source*dsin(theta)*dcos(phi)
  y_target_source = r_target_source*dsin(theta)*dsin(phi)
  z_target_source = r_target_source*dcos(theta)

  ! rotate force vector
  fx = comp_dir_vect_source_N(isource) * nu_source(1,1,isource) + &
       comp_dir_vect_source_E(isource) * nu_source(2,1,isource) + &
       comp_dir_vect_source_Z_UP(isource) * nu_source(3,1,isource)
  fy = comp_dir_vect_source_N(isource) * nu_source(1,2,isource) + &
       comp_dir_vect_source_E(isource) * nu_source(2,2,isource) + &
       comp_dir_vect_source_Z_UP(isource) * nu_source(3,2,isource)
  fz = comp_dir_vect_source_N(isource) * nu_source(1,3,isource) + &
       comp_dir_vect_source_E(isource) * nu_source(2,3,isource) + &
       comp_dir_vect_source_Z_UP(isource) * nu_source(3,3,isource)

  ! dimensionalize force
  scaleF = RHOAV * (R_EARTH**4) * PI*GRAV*RHOAV
  factor_force_source(:) = factor_force_source(:) * scaleF
  
  ! write to specfem cartesian FORCESOLUTION format
  !cart_force_fn = 'DATA/FORCESOLUTION_cartesian'
  open(unit=IOUT,file=trim(cart_force_fn), &
          status='unknown', form='formatted', action='write', iostat=ier)
  write(IOUT, "(a5,i4)") 'FORCE', 1
  write(IOUT, "(a11,f12.4)") 'time shift:', tshift_src(isource)
  write(IOUT, "(a3,f7.1)") 'f0:', hdur(isource)
  write(IOUT, "(a9,f25.9)") 'latorUTM:', y_target_source * R_EARTH
  write(IOUT, "(a10,f25.9)") 'longorUTM:', x_target_source * R_EARTH
  write(IOUT, "(a6,f25.9)") 'depth:', z_target_source * R_EARTH
  write(IOUT, "(a20,es15.1e3)") 'factor force source:', &
                                      factor_force_source(isource)
  write(IOUT, "(a28,f23.15)") 'component dir vect source E:', fx
  write(IOUT, "(a28,f23.15)") 'component dir vect source N:', fy
  write(IOUT, "(a31,f23.15)") 'component dir vect source Z_UP:', fz
  close(IOUT)

  deallocate(tshift_src, &
             hdur, &
             theta_source, &
             phi_source, &
             nu_source, &
             force_stf,factor_force_source, &
             comp_dir_vect_source_E, &
             comp_dir_vect_source_N, &
             comp_dir_vect_source_Z_UP)
  if (TOPOGRAPHY) deallocate(ibathy_topo)
  call MPI_Finalize(ier)
end program write_force_solution_file
