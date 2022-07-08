program cube2sph_station
  character(300) :: station_name,network_name,infn, outfn,string
  integer :: IIN, IOUT, ier
  double precision, dimension(3,1) :: nodes_coords, nodes_coords_new
  double precision :: r_earth=6371000.0,center_lat=62.5,&
          center_lon=-151.0,rotation_azi=20.0,rs,dummy
  double precision ::  lat, lon, r
  IIN = 20
  IOUT = 21
  infn = 'DATA/STATIONS_cart'
  outfn = 'DATA/STATIONS_sph'
  open(unit=IIN,file=trim(infn),status='old',action='read',iostat=ier)
  open(unit=IOUT,file=trim(outfn), &
          status='unknown', form='formatted', action='write', iostat=ier)
  ier = 0
  do while (ier == 0)
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) exit

    if (len_trim(string) > 0) then
      string = trim(string)
      read(string, *) station_name, network_name, nodes_coords(2,1), &
              nodes_coords(1,1), dummy, nodes_coords(3,1)

      call cube2sph_trans(nodes_coords,nodes_coords_new,1,&
            r_earth,center_lat,center_lon,rotation_azi)
      
      call xyz_2_rlatlon_dble(nodes_coords_new(1,1),nodes_coords_new(2,1),&
                          nodes_coords_new(3,1),r,lat,lon)
 
      !write(IOUT,'(a10,1x,a10,4e20.7)') trim(station_name),trim(network_name), &
      !         sngl(nodes_coords_new(2,1)),sngl(nodes_coords_new(1,1)), &
      !         0.0,sngl(nodes_coords_new(3,1))
      write(IOUT,'(a10,1x,a10,4e20.7)') trim(station_name),trim(network_name), &
               lat,lon, 0.0, 0.0
    endif
  enddo
  close(IIN)
  close(IOUT)
end program cube2sph_station
