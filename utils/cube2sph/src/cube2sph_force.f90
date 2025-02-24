program cube2sph_force
  !! ./cube2sph_force infn outfn center_lat center_lon rotation_azi
  implicit none
  character(len=300) :: string, infn, outfn, dummy
  integer :: dummyval, IFILE,ier
  double precision :: t_shift, hdur, factor_force_source,&
       comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP
  double precision, dimension(3,1) :: nodes_coords, nodes_coords_new, &
          nu, nu_new
  double precision :: r_earth=6371000.0,center_lat,center_lon,rotation_azi,rs
  double precision :: r, lat, lon
  call get_command_argument(1, infn)
  call get_command_argument(2, outfn)
  call get_command_argument(3, string)
  read(string, *) center_lat
  call get_command_argument(4, string)
  read(string, *) center_lon
  call get_command_argument(5, string)
  read(string, *) rotation_azi
  IFILE = 20
  !infn = 'DATA/FORCESOLUTION_cart'
  open(unit=IFILE,file=trim(infn),status='old',action='read',iostat=ier)
  read(IFILE,"(a)") string
  ! skips empty lines
  do while (len_trim(string) == 0)
  read(IFILE,"(a)") string
  enddo

  ! read header with event information
  read(string,"(a6,i4)") dummy,dummyval

  ! read time shift
  read(IFILE,"(a)") string
  read(string(12:len_trim(string)),*) t_shift

  ! read f0 (stored in hdur() array for convenience, to use the same array as for CMTSOLUTION)
  ! Please be careful, if you meet an error in reading the file FORCESOLUTION,
  ! such as you still write "hdur:" instead of "f0:"
  ! Please change your file or do following change in the code, such as changing
  ! read(string(4:len_trim(string)),*) hdur(isource)
  ! to
  ! read(string(6:len_trim(string)),*) hdur(isource)
  read(IFILE,"(a)") string
  read(string(4:len_trim(string)),*) hdur

  ! read latitude
  read(IFILE,"(a)") string
  read(string(10:len_trim(string)),*) nodes_coords(2,1)

  ! read longitude
  read(IFILE,"(a)") string
  read(string(11:len_trim(string)),*) nodes_coords(1,1)

  ! read depth
  read(IFILE,"(a)") string
  read(string(7:len_trim(string)),*) nodes_coords(3,1)

  ! read magnitude
  read(IFILE,"(a)") string
  read(string(21:len_trim(string)),*) factor_force_source

  ! read direction vector's East component
  read(IFILE,"(a)") string
  read(string(29:len_trim(string)),*) comp_dir_vect_source_E

  ! read direction vector's North component
  read(IFILE,"(a)") string
  read(string(29:len_trim(string)),*) comp_dir_vect_source_N

  ! read direction vector's vertical component
  read(IFILE,"(a)") string
  read(string(32:len_trim(string)),*) comp_dir_vect_source_Z_UP

  close(IFILE)

  call cube2sph_trans(nodes_coords,nodes_coords_new,1,&
          r_earth,center_lat,center_lon,rotation_azi)

  rs = sqrt(sum(nodes_coords_new(:,1)**2))
  comp_dir_vect_source_E = nodes_coords_new(1,1) / rs
  comp_dir_vect_source_N = nodes_coords_new(2,1) / rs
  comp_dir_vect_source_Z_UP = nodes_coords_new(3,1) / rs

  call xyz_2_rlatlon_dble(nodes_coords_new(1,1),nodes_coords_new(2,1),&
                          nodes_coords_new(3,1),r,lat,lon)

  !nu(1,1) = comp_dir_vect_source_E
  !nu(2,1) = comp_dir_vect_source_N
  !nu(3,1) = comp_dir_vect_source_Z_UP
  !call cube2sph_trans(nu,nu_new,1,&
  !        r_earth,center_lat,center_lon,rotation_azi)
  !comp_dir_vect_source_E = nu_new(1,1)
  !comp_dir_vect_source_N = nu_new(2,1)
  !comp_dir_vect_source_Z_UP = nu_new(3,1)


  !outfn = 'DATA/FORCESOLUTION_sph'
  open(unit=IFILE,file=trim(outfn), &
          status='unknown', form='formatted', action='write', iostat=ier)
  write(IFILE, "(a5,i4)") 'FORCE', 1
  write(IFILE, "(a11,f12.4)") 'time shift:', t_shift
  write(IFILE, "(a3,f7.1)") 'f0:', hdur
  !write(IFILE, "(a9,f25.9)") 'latorUTM:', nodes_coords_new(2,1)
  !write(IFILE, "(a10,f25.9)") 'longorUTM:', nodes_coords_new(1,1)
  !write(IFILE, "(a6,f25.9)") 'depth:', nodes_coords_new(3,1)
  write(IFILE, "(a9,f25.9)") 'latorUTM:', lat
  write(IFILE, "(a10,f25.9)") 'longorUTM:', lon
  write(IFILE, "(a6,f25.9)") 'depth:', 0.0
  write(IFILE, "(a20,es15.1e3)") 'factor force source:', &
                                      factor_force_source
  !write(IFILE, "(a28,f23.15)") 'component dir vect source E:', comp_dir_vect_source_E
  !write(IFILE, "(a28,f23.15)") 'component dir vect source N:', comp_dir_vect_source_N
  !write(IFILE, "(a31,f23.15)") 'component dir vect source Z_UP:', comp_dir_vect_source_Z_UP
  write(IFILE, "(a28,f23.15)") 'component dir vect source E:', 0.0
  write(IFILE, "(a28,f23.15)") 'component dir vect source N:', 0.0
  write(IFILE, "(a31,f23.15)") 'component dir vect source Z_UP:', 1.0
  close(IFILE)
end program cube2sph_force

