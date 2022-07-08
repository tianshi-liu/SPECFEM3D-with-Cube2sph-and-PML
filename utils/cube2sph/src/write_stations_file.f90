program write_stations_file
  use mpi
  use constants_solver, only: &
    !ELLIPTICITY_VAL, &
    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,NDIM, &
    HUGEVAL,IMAIN,IIN,IOUT,IOUT_VTK,MIDX,MIDY,MIDZ, &
    DEGREES_TO_RADIANS,RADIANS_TO_DEGREES,TWO_PI,R_UNIT_SPHERE, &
    R_EARTH,R_EARTH_KM,NX_BATHY,NY_BATHY,MAX_STRING_LEN
  use shared_input_parameters, only: OUTPUT_FILES
  use specfem_par, only: &
    ONE_CRUST, &
    myrank,DT,NSTEP, &
    STATIONS_FILE,nrec, &
    station_name,network_name, &
    stlat,stlon,stele,stbur,nu, &
    rspl,espl,espl2,nspl,ibathy_topo, &
    TOPOGRAPHY,RECEIVERS_CAN_BE_BURIED
  implicit none
  logical :: ELLIPTICITY
  integer :: iorientation
  double precision :: stazi,stdip
  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  integer :: irec
  integer :: ier
  double precision :: ell
  double precision :: elevation
  double precision :: n(3)
  double precision :: thetan,phin
  double precision :: sint,cost,sinp,cosp
  double precision :: r0,p20
  double precision :: theta,phi
  double precision :: x_target_rec,y_target_rec,z_target_rec
  character(len=256) :: string
  character(len=MAX_STRING_LEN) :: dummystring
  double precision :: lat,lon,r
  call MPI_Init(ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  STATIONS_FILE = 'DATA/STATIONS'
  ONE_CRUST = .true.
  TOPOGRAPHY = .true.
  ELLIPTICITY = .true.
  open(unit=IIN,file=trim(STATIONS_FILE),status='old',action='read',iostat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Stations file '//trim(STATIONS_FILE)//' could not be found, please check your setup')
  ! counts records
  nrec = 0
  do while(ier == 0)
    read(IIN,"(a)",iostat=ier) dummystring
    if (ier == 0) then
      ! excludes empty lines
      if (len_trim(dummystring) > 0 ) nrec = nrec + 1
    endif
  enddo
  close(IIN)
  allocate(station_name(nrec), &
           network_name(nrec), &
           stlat(nrec), &
           stlon(nrec), &
           stele(nrec), &
           stbur(nrec),stat=ier)
  allocate(nu(NDIM,NDIM,nrec),stat=ier)
  allocate(x_target(nrec), &
           y_target(nrec), &
           z_target(nrec),stat=ier)
  allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
  call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  ibathy_topo(:,:) = 0
  call read_topo_bathy_file(ibathy_topo)
  open(unit=IIN,file=trim(STATIONS_FILE),status='old',action='read',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening STATIONS file')

  ! loop on all the stations to read station information
  do irec = 1,nrec

    ! old line:
    !read(IIN,*,iostat=ier) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)

    ! reads in line as string
    read(IIN,"(a256)",iostat=ier) string
    if (ier /= 0) then
      write(IMAIN,*) 'Error reading in station ',irec
      call exit_MPI(myrank,'Error reading in station in STATIONS file')
    endif

    ! skips empty lines
    do while( len_trim(string) == 0 )
      read(IIN,"(a256)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading in station ',irec
        call exit_MPI(myrank,'Error reading in station in STATIONS file')
      endif
    enddo

    ! reads in station information
    read(string(1:len_trim(string)),*,iostat=ier) station_name(irec),network_name(irec), &
                                                  stlat(irec),stlon(irec),stele(irec),stbur(irec)
    if (ier /= 0) then
      write(IMAIN,*) 'Error reading in station ',irec
      call exit_MPI(myrank,'Error reading in station in STATIONS file')
    endif

    ! checks latitude
    if (stlat(irec) < -90.d0 .or. stlat(irec) > 90.d0) then
      write(IMAIN,*) 'Error station ',trim(station_name(irec)),': latitude ',stlat(irec), &
                     ' is invalid, please check STATIONS record'
      call exit_MPI(myrank,'Error station latitude invalid')
    endif

  enddo
  ! close receiver file
  close(IIN)
  STATIONS_FILE = 'DATA/STATIONS_cartesian'
  open(unit=IOUT,file=trim(STATIONS_FILE), &
          status='unknown', form='formatted', action='write', iostat=ier)
  do irec = 1,nrec


    ! station lat/lon in degrees
    lat = stlat(irec)
    lon = stlon(irec)

    ! limits longitude to [0.0,360.0]
    if (lon < 0.d0 ) lon = lon + 360.d0
    if (lon > 360.d0 ) lon = lon - 360.d0

    ! converts geographic latitude stlat (degrees) to geocentric colatitude theta (radians)
    call lat_2_geocentric_colat_dble(lat,theta)

    phi = lon*DEGREES_TO_RADIANS
    call reduce(theta,phi)

    ! record three components for each station
    do iorientation = 1,3
      !     North
      if (iorientation == 1) then
        stazi = 0.d0
        stdip = 0.d0
      !     East
      else if (iorientation == 2) then
        stazi = 90.d0
        stdip = 0.d0
      !     Vertical
      else if (iorientation == 3) then
        stazi = 0.d0
        stdip = - 90.d0
      else
        call exit_MPI(myrank,'incorrect orientation')
      endif

      !     get the orientation of the seismometer
      thetan=(90.0d0+stdip)*DEGREES_TO_RADIANS
      phin=stazi*DEGREES_TO_RADIANS

      ! we use the same convention as in Harvard normal modes for the orientation

      !     vertical component
      n(1) = cos(thetan)
      !     N-S component
      n(2) = - sin(thetan)*cos(phin)
      !     E-W component
      n(3) = sin(thetan)*sin(phin)

      !     get the Cartesian components of n in the model: nu
      sint = sin(theta)
      cost = cos(theta)
      sinp = sin(phi)
      cosp = cos(phi)
      nu(iorientation,1,irec) = n(1)*sint*cosp+n(2)*cost*cosp-n(3)*sinp
      nu(iorientation,2,irec) = n(1)*sint*sinp+n(2)*cost*sinp+n(3)*cosp
      nu(iorientation,3,irec) = n(1)*cost-n(2)*sint
    enddo

    ! normalized receiver radius
    r0 = R_UNIT_SPHERE

    ! finds elevation of receiver
    if (TOPOGRAPHY) then
       call get_topo_bathy(lat,lon,elevation,ibathy_topo)
       r0 = r0 + elevation/R_EARTH
    endif
    ! ellipticity
    if (ELLIPTICITY) then
      cost=cos(theta)
      ! this is the Legendre polynomial of degree two, P2(cos(theta)),
      ! see the discussion above eq (14.4) in Dahlen and Tromp (1998)
      p20=0.5d0*(3.0d0*cost*cost-1.0d0)
      ! get ellipticity using spline evaluation
      call spline_evaluation(rspl,espl,espl2,nspl,r0,ell)
      ! this is eq (14.4) in Dahlen and Tromp (1998)
      r0=r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
    endif

    ! subtract station burial depth (in meters)
    r0 = r0 - stbur(irec)/R_EARTH

    ! compute the Cartesian position of the receiver
    x_target_rec = r0*sin(theta)*cos(phi)
    y_target_rec = r0*sin(theta)*sin(phi)
    z_target_rec = r0*cos(theta)

    x_target(irec) = x_target_rec * R_EARTH
    y_target(irec) = y_target_rec * R_EARTH
    z_target(irec) = z_target_rec * R_EARTH
    write(IOUT, "(a4,a7,f23.4,f23.4,f5.1,f23.4)") network_name(irec), &
            station_name(irec), y_target(irec), x_target(irec), &
            0.0, z_target(irec) 
  enddo
  close(IOUT)
  STATIONS_FILE = 'DATA/rotation_nu'
  open(unit=IOUT,file=trim(STATIONS_FILE), &
          status='unknown', form='formatted', action='write', iostat=ier)
  do irec = 1,nrec
    write(IOUT, "(a4,a7)") network_name(irec), station_name(irec)
    write(IOUT, "(f23.9,f23.9,f23.9)") nu(1,1,irec),nu(1,2,irec),nu(1,3,irec)
    write(IOUT, "(f23.9,f23.9,f23.9)") nu(2,1,irec),nu(2,2,irec),nu(2,3,irec)
    write(IOUT, "(f23.9,f23.9,f23.9)") nu(3,1,irec),nu(3,2,irec),nu(3,3,irec)
  enddo
  close(IOUT)
  
  deallocate(station_name, &
           network_name, &
           stlat, &
           stlon, &
           stele, &
           stbur)
  deallocate(nu)
  deallocate(x_target, &
           y_target, &
           z_target)
  deallocate(ibathy_topo)
  call MPI_Finalize(ier)
end program write_stations_file

