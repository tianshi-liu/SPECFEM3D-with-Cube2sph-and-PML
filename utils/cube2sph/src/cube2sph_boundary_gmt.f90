program cube2sph_boundary_gmt
  implicit none
  character(300) :: outfn,line
  integer :: IOUT, i, ier
  integer :: N = 20
  double precision, dimension(3,1) :: nodes_coords, nodes_coords_new
  double precision, parameter :: r_earth=6371000.0,pi = atan(1.0) * 4.0
  double precision :: center_lat=43.5,&
          center_lon=121.8,rotation_azi=0.0,rs,dummy
  double precision ::  lat, lon, r, margin=2.0, &
                       x1=-10.0, x2=10.0, y1=-5.0, y2=5.0
  double precision, parameter :: meter_per_deg = r_earth * pi / 180.0

  if (command_argument_count() /= 8) then 
    print*, 'Usage: ./cube2sph_boundary_gmt cen_lat cen_lon rotangl xmin xmax ymin ymax margin'
    stop
  endif

  call get_command_argument(1,line)
  read(line,*) center_lat
  call get_command_argument(2,line)
  read(line,*) center_lon 
  call get_command_argument(3,line)
  read(line,*) rotation_azi
  call get_command_argument(4,line)
  read(line,*) x1 
  call get_command_argument(5,line)
  read(line,*) x2
  call get_command_argument(6,line)
  read(line,*) y1
  call get_command_argument(7,line)
  read(line,*) y2
  call get_command_argument(8,line)
  read(line,*),margin
  print*, 'cen_lat cen_lon rotangl = ',center_lat,center_lon,rotation_azi
  print*,'study region in deg = ', x1,x2,y1,y2
  print*,'margin in deg = ', margin

  IOUT = 21
  outfn = 'boundary_2d.gmt'
  open(unit=IOUT,file=trim(outfn), &
          status='unknown', form='formatted', action='write', iostat=ier)
  ! left side
  nodes_coords(1,1) = (x1 + margin) * meter_per_deg
  nodes_coords(3,1) = 0.0
  write(IOUT, '(a1)') '>'
  do i = 1, N
    nodes_coords(2,1) = ((y1 + margin) + (i-1) * (y2-y1-2.0*margin) / (N-1)) * &
                        meter_per_deg
    call cube2sph_trans(nodes_coords,nodes_coords_new,1,&
            r_earth,center_lat,center_lon,rotation_azi)
      
    call xyz_2_rlatlon_dble(nodes_coords_new(1,1),nodes_coords_new(2,1),&
                          nodes_coords_new(3,1),r,lat,lon)
    write(IOUT, '(2e20.7)') lon, lat
  enddo

  ! right side
  nodes_coords(1,1) = (x2 - margin) * meter_per_deg
  nodes_coords(3,1) = 0.0
  write(IOUT, '(a1)') '>'
  do i = 1, N
    nodes_coords(2,1) = ((y1 + margin) + (i-1) * (y2-y1-2.0*margin) / (N-1)) * &
                        meter_per_deg
    call cube2sph_trans(nodes_coords,nodes_coords_new,1,&
            r_earth,center_lat,center_lon,rotation_azi)
      
    call xyz_2_rlatlon_dble(nodes_coords_new(1,1),nodes_coords_new(2,1),&
                          nodes_coords_new(3,1),r,lat,lon)
    write(IOUT, '(2e20.7)') lon, lat
  enddo

  ! bottom side
  nodes_coords(2,1) = (y1 + margin) * meter_per_deg
  nodes_coords(3,1) = 0.0
  write(IOUT, '(a1)') '>'
  do i = 1, N
    nodes_coords(1,1) = ((x1 + margin) + (i-1) * (x2-x1-2.0*margin) / (N-1)) * &
                        meter_per_deg
    call cube2sph_trans(nodes_coords,nodes_coords_new,1,&
            r_earth,center_lat,center_lon,rotation_azi)
      
    call xyz_2_rlatlon_dble(nodes_coords_new(1,1),nodes_coords_new(2,1),&
                          nodes_coords_new(3,1),r,lat,lon)
    write(IOUT, '(2e20.7)') lon, lat
  enddo

  ! top side
  nodes_coords(2,1) = (y2 - margin) * meter_per_deg
  nodes_coords(3,1) = 0.0
  write(IOUT, '(a1)') '>'
  do i = 1, N
    nodes_coords(1,1) = ((x1 + margin) + (i-1) * (x2-x1-2.0*margin) / (N-1)) * &
                        meter_per_deg
    call cube2sph_trans(nodes_coords,nodes_coords_new,1,&
            r_earth,center_lat,center_lon,rotation_azi)
      
    call xyz_2_rlatlon_dble(nodes_coords_new(1,1),nodes_coords_new(2,1),&
                          nodes_coords_new(3,1),r,lat,lon)
    write(IOUT, '(2e20.7)') lon, lat
  enddo
  close(IOUT)
end program cube2sph_boundary_gmt
