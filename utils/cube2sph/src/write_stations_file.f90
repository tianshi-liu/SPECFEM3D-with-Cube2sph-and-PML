program write_stations_file
  !! writing STATIONS file in Cartesian coordinate
  !! Tianshi Liu, 2022.2.21
  !! Usage: write_stations_file OLD_FILE NEW_FILE ROTATION_NU_FILE TOPOGRAPHY ELLIPTICITY
  !! E.g., write_stations_file DATA/STATIONS_geo DATA/STATIONS DATA/rotation_nu .true. .true.
  !! takes the DATA/STATIONS_geo in geographic coordinates, writes
  !! DATA/STATIONS in Cartesian coordinate,  and
  !! DATA/rotation_nu the rotation matrix file
  !! with the Etopo topography and the ellipticity
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
  !double precision, allocatable, dimension(:) :: x_target,y_target,z_target
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
  character(len=MAX_STRING_LEN) :: dummystring, arg_str
  double precision :: lat,lon,r
  call MPI_Init(ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  call get_command_argument(1, STATIONS_FILE)
  !STATIONS_FILE = 'DATA/STATIONS'
  !ONE_CRUST = .true.
  call get_command_argument(4, arg_str)
  read(arg_str, *) TOPOGRAPHY
  call get_command_argument(5, arg_str)
  read(arg_str, *) ELLIPTICITY
  !TOPOGRAPHY = .true.
  !ELLIPTICITY = .true.
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
  call prepare_topography_ellipticity(TOPOGRAPHY)
  !if (TOPOGRAPHY) then
  !  allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
  !  call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  !  ibathy_topo(:,:) = 0
  !  call read_topo_bathy_file(ibathy_topo)
  !  call prepare_topography_ellipticity(TOPOGRAPHY)
  !endif
  
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
  call get_command_argument(2, STATIONS_FILE)
  !STATIONS_FILE = 'DATA/STATIONS_cartesian'
  open(unit=IOUT,file=trim(STATIONS_FILE), &
          status='unknown', form='formatted', action='write', iostat=ier)
  do irec = 1,nrec

    call geographic_to_cartesian(stlat(irec), stlon(irec), stele(irec), &
                                 stbur(irec), TOPOGRAPHY, ELLIPTICITY, &
                                 x_target_rec, y_target_rec, z_target_rec, &
                                 nu(:,:,irec))
    write(IOUT, "(a7,a7,f23.4,f23.4,f5.1,f23.4)") station_name(irec), &
            network_name(irec), y_target_rec, x_target_rec, &
            0.0, z_target_rec
  enddo
  close(IOUT)
  call get_command_argument(3, STATIONS_FILE)
  !STATIONS_FILE = 'DATA/rotation_nu'
  open(unit=IOUT,file=trim(STATIONS_FILE), &
          status='unknown', form='formatted', action='write', iostat=ier)
  do irec = 1,nrec
    write(IOUT, "(a7,a7)") network_name(irec), station_name(irec)
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
  if (TOPOGRAPHY) call deallocate_topography()
  call MPI_Finalize(ier)
end program write_stations_file

