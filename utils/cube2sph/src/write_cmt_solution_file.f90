program write_cmt_file
  !! writing CMT_SOLUTIONS file in Cartesian coordinate
  !! Tianshi Liu, 2022.2.21
  !! Usage: write_cmt_solution_file OLD_FILE NEW_FILE TOPOGRAPHY ELLIPTICITY
  !! E.g., write_cmt_solution_file DATA/CMTSOLUTIONS_geo DATA/CMTSOLUTIONS .true. .true.
  !! takes the DATA/CMTSOLUTIONS_geo in geographic coordinates, writes
  !! DATA/CMTSOLUTIONS in Cartesian coordinate, with the Etopo topography and
  !the ellipticity
  use mpi
  use constants_solver
  use specfem_par, only: &
      ONE_CRUST, &
      tshift_src,theta_source,phi_source, &
      DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz,Mw,M0, &
      rspl,espl,espl2,nspl,ibathy_topo, &
      LOCAL_TMP_PATH,SIMULATION_TYPE,TOPOGRAPHY, &
      nu_source, &
      USE_FORCE_POINT_SOURCE
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS, NOISE_TOMOGRAPHY
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
  character(len=MAX_STRING_LEN) :: cart_cmt_fn, sph_cmt_fn, arg_str
  double precision :: scaleM
  double precision :: min_tshift_src_original
  double precision :: radius
  integer :: yr, jda, mo, da, ho, mi
  double precision :: sec
  double precision, dimension(6,NSOURCES) :: moment_tensor
  double precision :: Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
  call MPI_Init(ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  USE_FORCE_POINT_SOURCE = .false.
  NOISE_TOMOGRAPHY = 0
  !NSOURCES = 1
  ONE_CRUST = .true.
  DT = 0.0
  NUMBER_OF_SIMULTANEOUS_RUNS = 1
  isource = 1
  call get_command_argument(1, sph_cmt_fn)
  call get_command_argument(2, cart_cmt_fn)
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
           Mxx(NSOURCES), &
           Myy(NSOURCES), &
           Mzz(NSOURCES), &
           Mxy(NSOURCES), &
           Mxz(NSOURCES), &
           Myz(NSOURCES), stat=ier)
  allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
  call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  ibathy_topo(:,:) = 0
  call read_topo_bathy_file(ibathy_topo)
  call get_cmt(yr,jda,mo,da,ho,mi,sec,tshift_src,hdur,lat,long,depth,&
               moment_tensor,DT,NSOURCES,min_tshift_src_original,sph_cmt_fn)
  ! dimensionalize
  scaleM = 1.d7 * RHOAV * (R_EARTH**5) * PI*GRAV*RHOAV
  moment_tensor(:,:) = moment_tensor(:,:) * scaleM
  call lat_2_geocentric_colat_dble(lat(isource),theta)

  phi = long(isource)*DEGREES_TO_RADIANS
  call reduce(theta,phi)
  st=dsin(theta)
  ct=dcos(theta)
  sp=dsin(phi)
  cp=dcos(phi)
  ! get the moment tensor
  Mrr = moment_tensor(1,isource)
  Mtt = moment_tensor(2,isource)
  Mpp = moment_tensor(3,isource)
  Mrt = moment_tensor(4,isource)
  Mrp = moment_tensor(5,isource)
  Mtp = moment_tensor(6,isource)

  Mxx(isource)=st*st*cp*cp*Mrr+ct*ct*cp*cp*Mtt+sp*sp*Mpp &
          +2.0d0*st*ct*cp*cp*Mrt-2.0d0*st*sp*cp*Mrp-2.0d0*ct*sp*cp*Mtp
  Myy(isource)=st*st*sp*sp*Mrr+ct*ct*sp*sp*Mtt+cp*cp*Mpp &
          +2.0d0*st*ct*sp*sp*Mrt+2.0d0*st*sp*cp*Mrp+2.0d0*ct*sp*cp*Mtp
  Mzz(isource)=ct*ct*Mrr+st*st*Mtt-2.0d0*st*ct*Mrt
  Mxy(isource)=st*st*sp*cp*Mrr+ct*ct*sp*cp*Mtt-sp*cp*Mpp &
          +2.0d0*st*ct*sp*cp*Mrt+st*(cp*cp-sp*sp)*Mrp+ct*(cp*cp-sp*sp)*Mtp
  Mxz(isource)=st*ct*cp*Mrr-st*ct*cp*Mtt &
          +(ct*ct-st*st)*cp*Mrt-ct*sp*Mrp+st*sp*Mtp
  Myz(isource)=st*ct*sp*Mrr-st*ct*sp*Mtt &
          +(ct*ct-st*st)*sp*Mrt+ct*cp*Mrp-st*cp*Mtp
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
  !fx = comp_dir_vect_source_N(isource) * nu_source(1,1,isource) + &
  !     comp_dir_vect_source_E(isource) * nu_source(2,1,isource) + &
  !     comp_dir_vect_source_Z_UP(isource) * nu_source(3,1,isource)
  !fy = comp_dir_vect_source_N(isource) * nu_source(1,2,isource) + &
  !     comp_dir_vect_source_E(isource) * nu_source(2,2,isource) + &
  !     comp_dir_vect_source_Z_UP(isource) * nu_source(3,2,isource)
  !fz = comp_dir_vect_source_N(isource) * nu_source(1,3,isource) + &
  !     comp_dir_vect_source_E(isource) * nu_source(2,3,isource) + &
  !     comp_dir_vect_source_Z_UP(isource) * nu_source(3,3,isource)

  
  ! write to specfem cartesian FORCESOLUTION format
  !cart_cmt_fn = 'DATA/CMTSOLUTION_cartesian'
  open(unit=IOUT,file=trim(cart_cmt_fn), &
          status='unknown', form='formatted', action='write', iostat=ier)
  write(IOUT, "(a3,i5,i3,i3,i3,i3,f6.2,f9.4,f10.4,f6.1)") 'PDE', &
              yr,mo,da,ho,mi,sec,lat(isource),long(isource),depth(isource)
  !write(IOUT, *)
  write(IOUT, "(a11,a12)") "event name:","     000000A"
  write(IOUT, "(a11,f12.4)") "time shift:", tshift_src(isource)
  write(IOUT, "(a14,f9.4)") "half duration:", hdur(isource)
  write(IOUT, "(a9,f14.4)") "latorUTM:", y_target_source*R_EARTH
  write(IOUT, "(a10,f13.4)") "longorUTM:", x_target_source*R_EARTH
  write(IOUT, "(a6,f17.4)") "depth:", z_target_source*R_EARTH
  write(IOUT, "(a4, es19.6E3)") "Mrr:", Mzz(isource)
  write(IOUT, "(a4, es19.6E3)") "Mtt:", Myy(isource)
  write(IOUT, "(a4, es19.6E3)") "Mpp:", Mxx(isource)
  write(IOUT, "(a4, es19.6E3)") "Mrt:", -Myz(isource)
  write(IOUT, "(a4, es19.6E3)") "Mrp:", Mxz(isource)
  write(IOUT, "(a4, es19.6E3)") "Mtp:", -Mxy(isource)
  close(IOUT)

  deallocate(tshift_src, &
             hdur, &
             theta_source, &
             phi_source, &
             nu_source, &
             Mxx, Myy, Mzz, Mxy, Mxz, Myz)
  deallocate(ibathy_topo)
  call MPI_Finalize(ier)
end program write_cmt_file
