program cube2latlon
  implicit none
  
  !! ./cube2sph_station infn outfn center_lat center_lon rotation_azi
  character(300) :: infn, outfn,string
  integer :: IIN, IOUT, ier
  double precision, dimension(3,1) :: nodes_coords, nodes_coords_new
  double precision :: r_earth=6371000.0,center_lat,center_lon,rotation_azi
  double precision ::  lat, lon, r

  if (command_argument_count() /= 5) then
    print*, 'Usage: ./this coordinates_cube(eta,xi,zeta) coordinates_sph(lat,lon,r) lat0 lon0 rot'
    stop
  endif
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
  !infn = 'STATIONS_cart'
  !outfn = 'STATIONS_sph.out'
  open(unit=IIN,file=trim(infn),status='old',action='read',iostat=ier)
  open(unit=IOUT,file=trim(outfn), &
          status='unknown', form='formatted', action='write', iostat=ier)
  ier = 0
  do while (ier == 0)
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) exit

    if (len_trim(string) > 0) then
      string = trim(string)
      read(string, *) nodes_coords(2,1), &
              nodes_coords(1,1), nodes_coords(3,1)

      call cube2sph_trans(nodes_coords,nodes_coords_new,1,&
            r_earth,center_lat,center_lon,rotation_azi)
      
      call xyz_2_rlatlon_dble(nodes_coords_new(1,1),nodes_coords_new(2,1),&
                          nodes_coords_new(3,1),r,lat,lon)
 
      !write(IOUT,'(a10,1x,a10,4e20.7)') trim(station_name),trim(network_name), &
      !         sngl(nodes_coords_new(2,1)),sngl(nodes_coords_new(1,1)), &
      !         0.0,sngl(nodes_coords_new(3,1))
      write(IOUT,'(g0,1x,g0,1x,g0)') lat,lon,r
    endif
  enddo
  close(IIN)
  close(IOUT)
end program cube2latlon
