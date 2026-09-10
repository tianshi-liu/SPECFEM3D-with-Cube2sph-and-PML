program cube2latlon
  use specfem_par, only: &
  ONE_CRUST, &
  myrank,DT,NSTEP, &
  STATIONS_FILE,nrec, &
  station_name,network_name, &
  stlat,stlon,stele,stbur,nu, &
  rspl,espl,espl2,nspl,ibathy_topo, &
  TOPOGRAPHY,RECEIVERS_CAN_BE_BURIED
  implicit none
  
  !! ./cube2sph_station infn outfn center_lat center_lon rotation_azi
  character(300) :: infn, outfn,string
  integer :: IIN, IOUT, ier
  double precision, dimension(3,1) :: nodes_coords, nodes_coords_new
  double precision :: r_earth=6371000.0,center_lat,center_lon,rotation_azi
  logical :: ENABLE_ELLIPTICITY
  double precision ::  lat, lon, r,dep 
  double precision :: theta,phi,cost,p20,ell,factor 

  if (command_argument_count() /= 5 .and. command_argument_count() /= 6) then
    print*, 'Usage: ./this coordinates_cube(eta,xi,zeta) coordinates_sph(lat,lon,r) lat0 lon0 rot (ellipticity optional)'
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
  if(command_argument_count() == 6) then
    call get_command_argument(6, string)
    read(string, *) ENABLE_ELLIPTICITY
  else
    ENABLE_ELLIPTICITY = .false.
  endif
  IIN = 20
  IOUT = 21

  if(ENABLE_ELLIPTICITY) then 
    call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  endif

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

      ! get x/y/z for sphere 
      call cube2sph_trans(nodes_coords,nodes_coords_new,1,&
            r_earth,center_lat,center_lon,rotation_azi)

      ! get r/theta/phi from x/y/z
      call xyz_2_rthetaphi_dble(nodes_coords_new(1,1), &
            nodes_coords_new(2,1),nodes_coords_new(3,1), &
            r,theta,phi)

      ! stretch x/y/z for ellipticity
      if(ENABLE_ELLIPTICITY) then
        cost = dcos(theta)
        p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)
        call spline_evaluation(rspl,espl,espl2,nspl,r/r_earth,ell)
        
        factor = 1.0d0-(2.0d0/3.0d0)*ell*p20
        print*,r,ell,factor,p20
        ! this is eq (14.4) in Dahlen and Tromp (1998)
        nodes_coords_new(:,1) = nodes_coords_new(:,1)*factor

        call get_geodetic_lat_lon_depth( &
            nodes_coords_new(1,1), &
            nodes_coords_new(2,1), &
            nodes_coords_new(3,1), &
            lat,lon,dep)
      else 

        ! nqdu added, convert to geodetic latitude, longitude, and depth for the node
        call geocentric_2_geographic_dble(theta,lat)
        lat = 90.0d0 - lat*180.0d0/dacos(-1.0d0)
        lon = phi*180.0d0/dacos(-1.0d0)
        dep = r_earth - r
      endif
 
      !write(IOUT,'(a10,1x,a10,4e20.7)') trim(station_name),trim(network_name), &
      !         sngl(nodes_coords_new(2,1)),sngl(nodes_coords_new(1,1)), &
      !         0.0,sngl(nodes_coords_new(3,1))
      write(IOUT,'(g0,1x,g0,1x,g0)') lat,lon,dep 
    endif
  enddo
  close(IIN)
  close(IOUT)
end program cube2latlon

subroutine get_geodetic_lat_lon_depth(x,y,z,lat,lon,dep)
  implicit none
  double precision,intent(in) :: x,y,z
  double precision,intent(out) :: lat,lon,dep
  double precision, parameter :: f = 1.0d0/298.257223563 ! WGS84 flattening
  double precision,parameter :: e2 = 2.0d0 * f - f*f ! eccentricity squared
  double precision,parameter :: ep2 = e2/(1.0d0 - e2) ! second eccentricity squared
  double precision,parameter :: r_earth=6371000.0
  double precision :: r,theta,phi,p,x_den,y_num,nn  

  double precision :: a,b  

  ! axis 
  a = r_earth * 3. / (3. - f)  ! semi-major axis
  b = a * (1. - f)              ! semi-minor axis
  p = sqrt(x*x + y*y)

  ! bowring's formula 
  theta = atan2(z*a,p*b)

  ! geodetic phi 
  y_num = z + ep2 * b * sin(theta)**3
  x_den = p - e2 * a * cos(theta)**3
  phi = atan2(y_num, x_den)

  ! calculate height 
  nn = a / sqrt(1.0d0 - e2 * sin(phi)**2)
  dep = -p / cos(phi) + nn

  ! lat/lon in degrees
  lat = phi * 180.0d0 / dacos(-1.0d0)
  lon = atan2(y, x) * 180.0d0 / dacos(-1.0d0)

end subroutine get_geodetic_lat_lon_depth