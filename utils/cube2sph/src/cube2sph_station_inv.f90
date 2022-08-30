program cube2sph_station_inv
  character(300) :: station_name,network_name,infn, outfn,string
  integer :: IIN, IOUT, ier
  double precision, dimension(3,1) :: nodes_coords, nodes_coords_new
  double precision, parameter :: r_earth = 6371000.0
  double precision, parameter :: DEGREES_TO_RADIANS = 3.1415926535 / 180.0
  double precision :: center_lat,center_lon,rotation_azi,rs,dummy
  double precision ::  lat, lon, r, elev, depth, theta, phi
  call get_command_argument(1, infn)
  call get_command_argument(2, outfn)
  call get_command_argument(3, string)
  read(string, *) center_lat
  call get_command_argument(4, string)
  read(string, *) center_lon
  call get_command_argument(5, string)
  read(string, *) rotation_azi
  IIN = 20
  IOUT = 21
  !infn = 'STATIONS_sph'
  !outfn = 'STATIONS_cart'
  open(unit=IIN,file=trim(infn),status='old',action='read',iostat=ier)
  open(unit=IOUT,file=trim(outfn), &
          status='unknown', form='formatted', action='write', iostat=ier)
  ier = 0
  do while (ier == 0)
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) exit

    if (len_trim(string) > 0) then
      string = trim(string)
      !read(string, *) station_name, network_name, nodes_coords(2,1), &
      !        nodes_coords(1,1), dummy, nodes_coords(3,1)
      read(string, *) station_name, network_name, lat, lon, elev, depth
      ! limits longitude to [0.0,360.0]
      if (lon < 0.d0 ) lon = lon + 360.d0
      if (lon > 360.d0 ) lon = lon - 360.d0
      
      ! converts geographic latitude stlat (degrees) to geocentric colatitude
      ! theta (radians)
      call lat_2_geocentric_colat_dble(lat,theta)
      phi = lon * DEGREES_TO_RADIANS
      call reduce(theta,phi)
      nodes_coords_new(1,1) = r_earth * sin(theta)*cos(phi)
      nodes_coords_new(2,1) = r_earth * sin(theta)*sin(phi)
      nodes_coords_new(3,1) = r_earth * cos(theta)

      call cube2sph_trans_inv(nodes_coords,nodes_coords_new,1,&
            r_earth,center_lat,center_lon,rotation_azi)
      
      !call xyz_2_rlatlon_dble(nodes_coords_new(1,1),nodes_coords_new(2,1),&
      !                    nodes_coords_new(3,1),r,lat,lon)
 
      !write(IOUT,'(a10,1x,a10,4e20.7)') trim(station_name),trim(network_name), &
      !         sngl(nodes_coords_new(2,1)),sngl(nodes_coords_new(1,1)), &
      !         0.0,sngl(nodes_coords_new(3,1))
      !write(IOUT,'(a10,1x,a10,4e20.7)') trim(station_name),trim(network_name), &
      !         lat,lon, 0.0, 0.0
      write(IOUT,'(a10,1x,a10,4e20.7)') trim(station_name),trim(network_name),&
               nodes_coords(2,1), nodes_coords(1,1), 0.0, nodes_coords(3,1)
    endif
  enddo
  close(IIN)
  close(IOUT)
end program cube2sph_station_inv
