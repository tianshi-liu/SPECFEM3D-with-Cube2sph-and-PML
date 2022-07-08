program div_kernel_by_nmeas
  use postprocess_par, only : MAX_KERNEL_NAMES, MAX_KERNEL_PATHS
  !use shared_parameters
  use specfem_par

  implicit none
  integer, parameter :: NARGS = 3
  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES), &
                                   meas_paths(MAX_KERNEL_PATHS), &
                                   kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: sline,kernel_dir,input_file,&
                                   kernel_name,filename
  character(len=MAX_STRING_LEN) :: arg(NARGS)
  integer :: npath, nker
  integer :: i, ier,iker,ipath
  integer :: nmeas,nmeas_all
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
                                   array
  double precision :: chi_dummy
  logical :: BROADCAST_AFTER_READ
  real(kind=CUSTOM_REAL) :: max_val,min_val

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  call synchronize_all()

  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo
  
  read(arg(1),'(a)') kernel_names_comma_delimited
  read(arg(2),'(a)') input_file
  read(arg(3),'(a)') kernel_dir

  ! parse names from KERNEL_NAMES
  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)

  ! parse paths from INPUT_FILE
  npath=0
  open(unit = IIN, file = trim(input_file), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(input_file), myrank
     stop 1
  endif
  do while (1 == 1)
     read(IIN,'(a)',iostat=ier) sline
     if (ier /= 0) exit
     npath = npath+1
     if (npath > MAX_KERNEL_PATHS) stop 'Error number of paths exceeds MAX_KERNEL_PATHS'
     meas_paths(npath) = sline
  enddo
  close(IIN)
  if (myrank == 0) then
    write(*,*) '  ',npath,' events'
    write(*,*)
  endif

  ! read simulation parameters
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! checks number of MPI processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *
      print *,'Expected number of MPI processes: ', NPROC
      print *,'Actual number of MPI processes: ', sizeprocs
      print *
    endif
    call synchronize_all()
    stop 'Error wrong number of MPI processes'
  endif
  call synchronize_all()

  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some
  ! array sizes below
  call read_mesh_for_init()

  allocate(array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)

  array = 0._CUSTOM_REAL

  call read_kernel_from_path(kernel_dir,kernel_names,nker,array)
  if (myrank == 0) print *, 'Before smoothing'
  do iker = 1, nker
    call min_all_cr(minval(array(:,:,:,:,iker)),min_val)
    call max_all_cr(maxval(array(:,:,:,:,iker)),max_val)
    if (myrank == 0) then
      print *, 'kernel ', trim(kernel_names(iker))
      print *, 'min:', min_val, 'max:', max_val
    endif
  enddo 

  nmeas_all = 0
  do ipath = 1, npath
    filename = trim(meas_paths(ipath)) // '/sum_chi'
    open(unit=IIN, file=trim(filename), form='formatted', status='old', &
         action='read', iostat=ier)
    read(IIN, *) chi_dummy
    read(IIN, *) nmeas
    close(IIN)
    if (myrank == 0) then
      print *, 'measurement path:', trim(meas_paths(ipath))
      print *, '  nmeas = ', nmeas
    endif
    nmeas_all = nmeas_all + nmeas
  enddo
  if (myrank == 0) then
    print *, 'total number of measurements:', nmeas_all
  endif

  array(:,:,:,:,:) = array(:,:,:,:,:) / (nmeas_all * 1.0)
  if (myrank == 0) print *, 'After smoothing'
  do iker = 1, nker
    call min_all_cr(minval(array(:,:,:,:,iker)),min_val)
    call max_all_cr(maxval(array(:,:,:,:,iker)),max_val)
    if (myrank == 0) then
      print *, 'kernel ', trim(kernel_names(iker))
      print *, 'min:', min_val, 'max:', max_val
    endif
  enddo 

  do iker = 1, nker
    kernel_name = kernel_names(iker)
    write(filename,'(a,i6.6,a)') trim(kernel_dir) // &
                '/proc',myrank,'_'// trim(kernel_name)//'_norm.bin'
    open(IOUT,file=trim(filename),form='unformatted',action='write',&
              status='unknown', iostat=ier)
    write(IOUT) array(:,:,:,:,iker)
    close(IOUT)
  enddo

  deallocate(array)
  call finalize_mpi()
end program div_kernel_by_nmeas

subroutine read_kernel_from_path(kernel_path,kernel_names,nker,array)
  use postprocess_par, only: MAX_KERNEL_NAMES, MAX_KERNEL_PATHS
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_AB, &
                         MAX_STRING_LEN, IIN, myrank
  implicit none
  integer, intent(in) :: nker
  character(len=MAX_STRING_LEN), intent(in) :: kernel_path
  character(len=MAX_STRING_LEN), intent(in) :: kernel_names(MAX_KERNEL_NAMES)
  real(kind=CUSTOM_REAL), intent(out) :: array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker)
  ! local variables
  integer :: iker, ier
  character(len=MAX_STRING_LEN) :: kernel_name, filename

  do iker = 1, nker
    kernel_name = kernel_names(iker)
    write(filename,'(a,i6.6,a)') trim(kernel_path) &
                           //'/proc',myrank,'_'//trim(kernel_name)//'.bin'
    open(IIN,file=trim(filename),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) '  array not found: ',trim(filename)
      stop 'Error array file not found'
    endif
    read(IIN) array(:,:,:,:,iker)
    close(IIN)
  enddo

end subroutine read_kernel_from_path
