program get_isotropic_vbulk
  use specfem_par
  implicit none
  integer, parameter :: NARGS = 5
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                               vp_array, vs_array, vbulk_array
  character(len=MAX_STRING_LEN) :: filename, dir_in, vp_fn, vs_fn, &
                                   dir_out,vbulk_fn
  character(len=MAX_STRING_LEN) :: arg(NARGS)
  logical :: BROADCAST_AFTER_READ
  integer :: i, ier  

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  call synchronize_all()
  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo
  read(arg(1),'(a)') dir_in
  read(arg(2),'(a)') vp_fn
  read(arg(3),'(a)') vs_fn
  read(arg(4),'(a)') dir_out
  read(arg(5),'(a)') vbulk_fn

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

  allocate(vp_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  allocate(vs_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  allocate(vbulk_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)

  vp_array = 0._CUSTOM_REAL
  vs_array = 0._CUSTOM_REAL
  vbulk_array = 0._CUSTOM_REAL

  ! read vp
  write(filename, '(a,i6.6,a)') trim(dir_in) &
               //'/proc',myrank,'_'//trim(vp_fn)//'.bin'
  open(IIN,file=trim(filename),status='old',form='unformatted',&
       action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) '  array not found: ',trim(filename)
    stop 'Error array file not found'
  endif
  read(IIN) vp_array(:,:,:,:)
  close(IIN)

  ! read vs
  write(filename, '(a,i6.6,a)') trim(dir_in) &
               //'/proc',myrank,'_'//trim(vs_fn)//'.bin'
  open(IIN,file=trim(filename),status='old',form='unformatted',&
       action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) '  array not found: ',trim(filename)
    stop 'Error array file not found'
  endif
  read(IIN) vs_array(:,:,:,:)
  close(IIN)

  vbulk_array(:,:,:,:) = sqrt(vp_array(:,:,:,:)*vp_array(:,:,:,:) -&
                       4.0*vs_array(:,:,:,:) * vs_array(:,:,:,:) / 3.0)

  ! write vbulk
  write(filename, '(a,i6.6,a)') trim(dir_out) &
               //'/proc',myrank,'_'//trim(vbulk_fn)//'.bin'
  open(IOUT,file=trim(filename),status='unknown',form='unformatted',&
       action='write',iostat=ier)
  write(IOUT) vbulk_array
  close(IOUT)

  deallocate(vp_array,vs_array,vbulk_array)
  call finalize_mpi()
  
end program get_isotropic_vbulk
