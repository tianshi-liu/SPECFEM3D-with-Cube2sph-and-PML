program write_zero_array
  use postprocess_par, only: MAX_KERNEL_NAMES
  use specfem_par
  implicit none
  integer, parameter :: NARGS = 2
  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: dir_out, kernel_name,&
                                   kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: arg(NARGS)

  integer :: nker, i, ier, iker
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
                                   ker_array
  logical :: BROADCAST_AFTER_READ

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  call synchronize_all()

  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo

  read(arg(1),'(a)') dir_out
  read(arg(2),'(a)') kernel_names_comma_delimited

  ! parse names from KERNEL_NAMES
  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)

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

  allocate(ker_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)

  ker_array = 0._CUSTOM_REAL


  do iker = 1, nker
    kernel_name = kernel_names(iker)
    write(filename,'(a,i6.6,a)') trim(dir_out) // &
                '/proc',myrank,'_'// trim(kernel_name)//'.bin'
    open(IOUT,file=trim(filename),form='unformatted',action='write',&
              status='unknown', iostat=ier)
    write(IOUT) ker_array(:,:,:,:,iker)
    close(IOUT)
  enddo

  deallocate(ker_array)

  call finalize_mpi()

end program write_zero_array

