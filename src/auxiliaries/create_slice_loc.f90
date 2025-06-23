program create_slice_loc
  use mpi
  use constants
  use shared_parameters
  use specfem_par, only: xigll, yigll, zigll, wxgll, wygll, wzgll
  use projection_on_FD_grid
  implicit none
  logical, parameter :: NEAREST_GLL_POINT = .false.
  integer, parameter :: NARGS = 3, MAX_POINTS_IN_SLICE = 20000000
  character(len=MAX_STRING_LEN) slice_fn, mesh_dir, out_fn, prname
  integer :: myrank, sizeprocs, NSPEC_IRREGULAR, ios, err, ier, ourtag=1, irank

  real(kind=CUSTOM_REAL), dimension(MAX_POINTS_IN_SLICE, 3) :: slice_coords
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            x_min_glob, x_max_glob, &
                            y_min_glob, y_max_glob, &
                            z_min_glob, z_max_glob
  double precision :: xi_loc, eta_loc, gamma_loc, x, y, z,&
                            x_found, y_found, z_found
  real(kind=CUSTOM_REAL) :: xmin, xmax, ymin, ymax, zmin, zmax, dist_found,dist
  real(kind=CUSTOM_REAL) :: xmin_proc, xmax_proc, ymin_proc, ymax_proc, &
                            zmin_proc, zmax_proc
  integer :: npts, ipt, ispec, iglob, i, j, k, ispec_selected, n_pts_proc
  integer, dimension(:), allocatable :: iaddx, iaddy, iaddz
  integer, dimension(MAX_POINTS_IN_SLICE) :: ispec_rec, ipt_rec
  real(kind=CUSTOM_REAL), dimension(MAX_POINTS_IN_SLICE) :: xi_rec, eta_rec, &
                                                  gamma_rec, dist_rec
  integer :: n_pts0, i_found, j_found, k_found
  integer, dimension(MAX_POINTS_IN_SLICE) :: ispec0, ipt0
  real(kind=CUSTOM_REAL), dimension(MAX_POINTS_IN_SLICE) :: xi0, eta0, &
                                                  gamma0, dist0
  integer :: n_pts_this_rank
  integer, dimension(MAX_POINTS_IN_SLICE) :: ipt_this_rank
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)
  call get_command_argument(1, slice_fn)
  !call get_command_argument(2, mesh_dir)
  call get_command_argument(2, out_fn)
  call read_parameter_file(myrank, .true.)
  if (myrank==0) print *, 'NGNOD=', NGNOD
  allocate(iaddx(NGNOD), iaddy(NGNOD), iaddz(NGNOD))
  mesh_dir = 'DATABASES_MPI'
  n_pts_proc = 0
  open(11, file=slice_fn, iostat=ios)
  ipt = 0
  do while (1 == 1)
    ipt=ipt+1
    read(11,*,iostat=ios) slice_coords(ipt, 1), slice_coords(ipt, 2), &
                          slice_coords(ipt, 3)
    if (ios /= 0) exit
  enddo
  close(11)
  npts=ipt-1
  if (myrank == 0) then
    if (npts > MAX_POINTS_IN_SLICE .or. npts < 0) &
            call exit_mpi(myrank,'Npts error ...')
    write(*,*) 'Total number of points = ', npts
  endif
  write(prname,'(a,i6.6,a)') trim(mesh_dir)//'/proc',myrank,'_'
  open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  read(27) NSPEC_AB
  read(27) NGLOB_AB
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1102')
  if (ier /= 0) stop 'error allocating array ibool'
  allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1103')
  if (ier /= 0) stop 'error allocating array xstore etc.'

  read(27) NSPEC_IRREGULAR
  read(27) ibool
  read(27) xstore
  read(27) ystore
  read(27) zstore
  close(27)
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
         x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
         elemsize_min_glob,elemsize_max_glob, &
         distance_min_glob,distance_max_glob)
  n_pts_this_rank = 0
  xmin_proc = minval(xstore)
  xmax_proc = maxval(xstore)
  ymin_proc = minval(ystore)
  ymax_proc = maxval(ystore)
  zmin_proc = minval(zstore)
  zmax_proc = maxval(zstore)
  if (myrank == 0) then
    print *, 'xmin=', xmin_proc, ', xmax=', xmax_proc
    print *, 'ymin=', ymin_proc, ', ymax=', ymax_proc
    print *, 'zmin=', zmin_proc, ', zmax=', zmax_proc
  endif


  do ipt = 1, npts
    x = slice_coords(ipt, 1)
    y = slice_coords(ipt, 2)
    z = slice_coords(ipt, 3)
    if ((x < xmin_proc) .or. (x > xmax_proc)) cycle
    if ((y < ymin_proc) .or. (y > ymax_proc)) cycle
    if ((z < zmin_proc) .or. (z > zmax_proc)) cycle
    n_pts_this_rank = n_pts_this_rank + 1
    ipt_this_rank(n_pts_this_rank) =  ipt
  enddo
  call synchronize_all()
  print *, 'found ', n_pts_this_rank, ' in rank ', myrank
  do ispec = 1, NSPEC_AB
    xmin = HUGEVAL
    xmax = -HUGEVAL
    ymin =  HUGEVAL
    ymax = -HUGEVAL
    zmin =  HUGEVAL
    zmax = -HUGEVAL
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          xmin  = min( xmin, xstore(iglob))
          xmax  = max( xmax, xstore(iglob))
          ymin  = min( ymin, ystore(iglob))
          ymax  = max( ymax, ystore(iglob))
          zmin  = min( zmin, zstore(iglob))
          zmax  = max( zmax, zstore(iglob))
        enddo
      enddo
    enddo
    do ipt = 1, n_pts_this_rank
      x = slice_coords(ipt_this_rank(ipt), 1)
      y = slice_coords(ipt_this_rank(ipt), 2)
      z = slice_coords(ipt_this_rank(ipt), 3)
      if ((x < xmin) .or. (x > xmax)) cycle
      if ((y < ymin) .or. (y > ymax)) cycle
      if ((z < zmin) .or. (z > zmax)) cycle
      if (NEAREST_GLL_POINT) then
        dist_found = HUGEVAL
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool(i,j,k,ispec)
              dist = sqrt((x-xstore(iglob))**2+(y-ystore(iglob))**2 &
                          +(z-zstore(iglob))**2)
              if (dist_found > dist) then
                i_found = i
                j_found = j
                k_found = k
                dist_found = dist
              endif
            enddo
          enddo
        enddo
        ispec_selected = ispec
        xi_loc = xigll(i_found)
        eta_loc = yigll(j_found)
        gamma_loc = zigll(k_found)
        iglob = ibool(i_found, j_found, k_found, ispec_selected)
        x_found = xstore(iglob)
        y_found = ystore(iglob)
        z_found = zstore(iglob)
      else
        call locate_point_in_element(x, y, z, iaddx, iaddy, iaddz, &
           elemsize_max_glob, ispec_selected, xi_loc, eta_loc, gamma_loc, &
                                   x_found, y_found, z_found, myrank, ispec)
      endif
      if (abs(xi_loc) < 1.05d0 .and. abs(eta_loc) < 1.05d0 &
                          .and. abs(gamma_loc) < 1.05d0) then
        n_pts_proc = n_pts_proc + 1
        ispec_rec(n_pts_proc) = ispec_selected
        xi_rec(n_pts_proc) = xi_loc
        eta_rec(n_pts_proc) = eta_loc
        gamma_rec(n_pts_proc) = gamma_loc
        ipt_rec(n_pts_proc) = ipt_this_rank(ipt)
        dist_rec(n_pts_proc) = sqrt((x-x_found)**2+(y-y_found)**2&
                                               +(z-z_found)**2)
      endif
    enddo
  enddo
  call synchronize_all()
  !if (myrank == 0) print *, 'in total ', n_pts_proc, ' times in rank 0'
  if (myrank == 0) print *, 'write to ', out_fn
  if (myrank == 0) open(11, file=out_fn, iostat=ios, action='write')
  if (myrank == 0) then
    do ipt = 1, n_pts_proc
      write(11, '(3i12, 4e18.6)')ipt_rec(ipt), 0, ispec_rec(ipt), &
               xi_rec(ipt), eta_rec(ipt), gamma_rec(ipt), dist_rec(ipt)
    enddo
  else
    call MPI_Ssend(n_pts_proc, 1, MPI_INT, 0, ourtag, &
                   MPI_COMM_WORLD, err)
    call MPI_Ssend(ispec_rec, n_pts_proc, MPI_INT, 0, ourtag, &
                   MPI_COMM_WORLD, err)
    call MPI_Ssend(ipt_rec, n_pts_proc, MPI_INT, 0, ourtag, &
                   MPI_COMM_WORLD, err)
    call MPI_Ssend(xi_rec, n_pts_proc, MPI_REAL, 0, ourtag, &
                   MPI_COMM_WORLD, err)
    call MPI_Ssend(eta_rec, n_pts_proc, MPI_REAL, 0, ourtag, &
                   MPI_COMM_WORLD, err)
    call MPI_Ssend(gamma_rec, n_pts_proc, MPI_REAL, 0, ourtag, &
                   MPI_COMM_WORLD, err)
    call MPI_Ssend(dist_rec, n_pts_proc, MPI_REAL, 0, ourtag, &
                   MPI_COMM_WORLD, err)
  endif
  print *, 'in total ', n_pts_proc, ' times in rank ', myrank
  if (myrank == 0) then
    do irank = 1, sizeprocs - 1
      call MPI_Recv(n_pts0, 1, MPI_INT, irank, ourtag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      call MPI_Recv(ispec0, n_pts0, MPI_INT, irank, ourtag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      call MPI_Recv(ipt0, n_pts0, MPI_INT, irank, ourtag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      call MPI_Recv(xi0, n_pts0, MPI_REAL, irank, ourtag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      call MPI_Recv(eta0, n_pts0, MPI_REAL, irank, ourtag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      call MPI_Recv(gamma0, n_pts0, MPI_REAL, irank, ourtag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      call MPI_Recv(dist0, n_pts0, MPI_REAL, irank, ourtag, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
      do ipt = 1, n_pts0
        write(11, '(3i12, 4e18.6)')ipt0(ipt), irank, ispec0(ipt), &
                 xi0(ipt), eta0(ipt), gamma0(ipt), dist0(ipt)
      enddo
    enddo
    close(11)
  endif
!  do irank = 0, sizeprocs - 1
!    if (irank == 0) then
!      if (myrank == 0) then
!      do ipt = 1, n_pts_proc
!        write(11, '(3i10, 3e18.6)')ipt_rec(ipt), irank, ispec_rec(ipt), &
!                 xi_rec(ipt), eta_rec(ipt), gamma_rec(ipt)
!      enddo
!      endif
!    else
!      call MPI_Sendrecv(n_pts_proc, 1, MPI_INT, 0, ourtag, &
!                        n_pts0, 1, MPI_INT, irank, ourtag, &
!                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
!      call MPI_Sendrecv(ispec_rec, MAX_POINTS_IN_SLICE*8, MPI_INT, 0, ourtag, &
!                        ispec0, MAX_POINTS_IN_SLICE*8, MPI_INT, irank, ourtag, &
!                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
!      call MPI_Sendrecv(ipt_rec, MAX_POINTS_IN_SLICE*8, MPI_INT, 0, ourtag, &
!                        ipt0, MAX_POINTS_IN_SLICE*8, MPI_INT, irank, ourtag, &
!                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
!      call MPI_Sendrecv(xi_rec, MAX_POINTS_IN_SLICE*8, MPI_REAL, 0, ourtag, &
!                        xi0, MAX_POINTS_IN_SLICE*8, MPI_REAL, irank, ourtag, &
!                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
!      call MPI_Sendrecv(eta_rec, MAX_POINTS_IN_SLICE*8, MPI_REAL, 0, ourtag, &
!                        eta0, MAX_POINTS_IN_SLICE*8, MPI_REAL, irank, ourtag, &
!                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
!      call MPI_Sendrecv(gamma_rec, MAX_POINTS_IN_SLICE*8, MPI_REAL, 0, ourtag, &
!                        gamma0, MAX_POINTS_IN_SLICE*8, MPI_REAL, irank, ourtag, &
!                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
!      if (myrank == 0) then
!      do ipt = 1, n_pts0
!        write(11, '(3i10, 3e18.6)')ipt0(ipt), irank, ispec0(ipt), &
!                 xi0(ipt), eta0(ipt), gamma0(ipt)
!      enddo
!      endif
!    endif
!  enddo
!  if (myrank == 0) close(11)
  deallocate(iaddx, iaddy, iaddz)
  deallocate(ibool, xstore, ystore, zstore)
  call finalize_mpi()
end program create_slice_loc
