program create_slice
  use,intrinsic :: ieee_arithmetic
  use mpi
  use constants
  use shared_parameters, only: NGNOD
  use specfem_par, only: xigll, yigll, zigll, wxgll, wygll, wzgll, &
            NSPEC_IRREGULAR, NSPEC_AB, NGLOB_AB, ibool, xstore, ystore, zstore
  implicit none
  logical, parameter :: USE_GLL_POINTS = .true.
  logical, parameter :: LINEAR_INTERP = .false.
  !logical, parameter :: normalize = .true.
  integer, parameter :: MAX_POINTS_IN_SLICE = 10000000
  integer, dimension(:), allocatable :: iaddx, iaddy, iaddz
  double precision, dimension(:), allocatable :: shape3D
  double precision, dimension(NGLLX) :: hxis, hpxis
  double precision, dimension(NGLLY) :: hetas, hpetas
  double precision, dimension(NGLLZ) :: hgammas, hpgammas
  double precision :: xi, eta, gamma, dist, v
  real(kind=CUSTOM_REAL), dimension(MAX_POINTS_IN_SLICE) :: val, dist_rec
  real(kind=CUSTOM_REAL), dimension(3) :: real_dummy 
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable ::  database
  real(kind=CUSTOM_REAL) :: FILL_VAL
  integer :: irank, current_rank, ios, ispec, iglob, ier
  character(len=MAX_STRING_LEN) slice_fn, loc_fn, out_fn, prname, &
                                model_name, model_dir, norm_string,&
                                database_dir
  logical :: normalize
  integer :: ipt, npts, i, j, k
  real(kind=CUSTOM_REAL) :: norm
  integer :: igll, jgll, kgll
  double precision :: ratio_x, ratio_y, ratio_z
  call MPI_INIT(ier)
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  current_rank = -1

  if ( command_argument_count() /= 7) then 
    print*,'usage: ./xcreate_slice model_name model_dir database_dir coordfile glllocfile  outfile do_norm'
    stop
  endif

  ! set value of fill_val
  FILL_VAL = ieee_value(FILL_VAL,ieee_quiet_nan)

  call get_command_argument(1, model_name)
  call get_command_argument(2, model_dir)
  call get_command_argument(3, database_dir)
  call get_command_argument(4, slice_fn)
  call get_command_argument(5, loc_fn)
  call get_command_argument(6, out_fn)
  call get_command_argument(7, norm_string)
  read(norm_string, *) normalize
  !call read_parameter_file(0, .true.)
  NGNOD = 27
  print *, 'NGNOD=', NGNOD
  allocate(iaddx(NGNOD), iaddy(NGNOD), iaddz(NGNOD), shape3D(NGNOD))
  npts = 0
  dist_rec(:) = HUGEVAL
  val(:) = FILL_VAL
  open(11, file=slice_fn, iostat=ios)
  do while (1 == 1)
    npts = npts + 1
    read(11, *, iostat=ios) real_dummy(1), real_dummy(2), real_dummy(3)
    if(ios /= 0) exit
  enddo
  close(11)
  npts = npts - 1
  open(11, file=loc_fn, iostat=ios)
  do while(1 == 1)
    read(11, *, iostat=ios) ipt, irank, ispec, xi, eta, gamma, dist
    if (ios /= 0) exit
    if (irank /= current_rank) then
      !if (current_rank >= 0) deallocate(ibool, xstore, ystore, zstore, database)
      if (allocated(ibool)) &
        deallocate(ibool, xstore, ystore, zstore, database)
      write(prname,'(a,i6.6,a)') TRIM(database_dir) // '/proc',irank,'_'
      open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
      read(27) NSPEC_AB
      read(27) NGLOB_AB
      allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
      read(27) NSPEC_IRREGULAR
      read(27) ibool
      read(27) xstore
      read(27) ystore
      read(27) zstore
      close(27)
      allocate(database(NGLLX, NGLLY, NGLLZ, NSPEC_AB))
      write(prname,'(a,i6.6,a)') trim(model_dir)//'/proc',irank,'_'
      open(unit = 27,file = trim(prname)//trim(model_name)//'.bin', &
           status='old', action='read', iostat=ios,form ='unformatted')
      if (ios /= 0) print *, 'error opening database file ', &
                                    trim(prname)//trim(model_name)//'.bin'
      print *, 'read file ', trim(prname)//trim(model_name)//'.bin'
      read(27) database
      close(27) 
      current_rank = irank
    endif
    if (dist > dist_rec(ipt)) cycle
    v = 0.d0
    if (USE_GLL_POINTS) then
      call lagrange_any(xi, NGLLX, xigll, hxis, hpxis)
      call lagrange_any(eta, NGLLY, yigll, hetas, hpetas)
      call lagrange_any(gamma, NGLLZ, zigll, hgammas, hpgammas)
      if (LINEAR_INTERP) then
        igll = 1
        do i = 2, NGLLX-1
          if (xi > xigll(i)) igll = i
        enddo
        ratio_x = (xi - xigll(igll)) / (xigll(igll+1)-xigll(igll))
        jgll = 1
        do j = 2, NGLLY-1
          if (eta > yigll(j)) jgll = j
        enddo
        ratio_y = (eta - yigll(jgll)) / (yigll(jgll+1)-yigll(jgll))
        kgll = 1
        do k = 2, NGLLZ-1
          if (gamma > zigll(k)) kgll = k
        enddo
        ratio_z = (gamma - zigll(kgll)) / (zigll(kgll+1)-zigll(kgll))
        v=dble(database(igll,jgll,kgll,ispec))*(1.0-ratio_x)*(1.0-ratio_y)*(1.0-ratio_z) + &
          dble(database(igll+1,jgll,kgll,ispec))*ratio_x*(1.0-ratio_y)*(1.0-ratio_z)+ &
          dble(database(igll,jgll+1,kgll,ispec))*(1.0-ratio_x)*ratio_y*(1.0-ratio_z)+ &
          dble(database(igll+1,jgll+1,kgll,ispec))*ratio_x*ratio_y*(1.0-ratio_z) + &
          dble(database(igll,jgll,kgll+1,ispec))*(1.0-ratio_x)*(1.0-ratio_y)*ratio_z+ &
          dble(database(igll+1,jgll,kgll+1,ispec))*ratio_x*(1.0-ratio_y)*ratio_z + &
          dble(database(igll,jgll+1,kgll+1,ispec))*(1.0-ratio_x)*ratio_y*ratio_z + &
          dble(database(igll+1,jgll+1,kgll+1,ispec))*ratio_x*ratio_y*ratio_z
      else
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool(i, j, k, ispec)
              v = v + hxis(i) * hetas(j) * hgammas(k) * dble(database(i,j,k,ispec))
            enddo
          enddo
        enddo
      endif
    else
      call usual_hex_nodes(NGNOD, iaddx, iaddy, iaddz)
      call eval_shape3D_single(0, shape3D, xi, eta, gamma, NGNOD)
      iaddx = iaddx * 2 + 1
      iaddy = iaddy * 2 + 1
      iaddz = iaddz * 2 + 1
      do i = 1, NGNOD
        iglob = ibool(iaddx(i), iaddy(i), iaddz(i), ispec)
        v = v + shape3D(i) * dble(database(iaddx(i),iaddy(i),iaddz(i),ispec))
      enddo
    endif
    val(ipt) = v
    dist_rec(ipt) = dist
    if (abs(v) > norm) norm = abs(v)
  enddo
  close(11)
  open(12, file=out_fn, iostat=ios, action='write')
  print *, 'max abs:', norm
  do ipt = 1, npts
    if (normalize) then
      write(12, '(2(G0,1x))') val(ipt) / norm, dist_rec(ipt)
    else
      write(12,'(2(G0,1x))') val(ipt), dist_rec(ipt)
    endif
  enddo
  close(12)
  !deallocate(ibool, xstore, ystore, zstore, database)
  if (allocated(ibool)) &
        deallocate(ibool, xstore, ystore, zstore, database)
  deallocate(iaddx, iaddy, iaddz, shape3D)
  call MPI_FINALIZE(ier)
end program create_slice
