program scale_rho_with_vs
  use specfem_par
  implicit none
  integer, parameter :: NARGS = 7
  real(kind=CUSTOM_REAL), parameter :: SCALE_RHO = 0.33
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    vsv_array, vsh_array, vs_voigt_sq_array, dvsv_array, dvsh_array, drho_array
  character(len=MAX_STRING_LEN) :: filename, kernel_dir, model_dir,&
                                   vsv_fn, vsh_fn, &
                                   dvsv_fn, dvsh_fn, drho_fn
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
  read(arg(1),'(a)') model_dir
  read(arg(2),'(a)') vsv_fn
  read(arg(3),'(a)') vsh_fn
  read(arg(4),'(a)') kernel_dir
  read(arg(5),'(a)') dvsv_fn
  read(arg(6),'(a)') dvsh_fn
  read(arg(7),'(a)') drho_fn

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

  allocate(vsv_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  allocate(vsh_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  allocate(vs_voigt_sq_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  allocate(dvsv_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  allocate(dvsh_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  allocate(drho_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)

  vsv_array = 0._CUSTOM_REAL
  vsh_array = 0._CUSTOM_REAL
  vs_voigt_sq_array = 0._CUSTOM_REAL
  dvsv_array = 0._CUSTOM_REAL
  dvsh_array = 0._CUSTOM_REAL
  drho_array = 0._CUSTOM_REAL

  ! read vsv
  write(filename, '(a,i6.6,a)') trim(model_dir) &
               //'/proc',myrank,'_'//trim(vsv_fn)//'.bin'
  open(IIN,file=trim(filename),status='old',form='unformatted',&
       action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) '  array not found: ',trim(filename)
    stop 'Error array file not found'
  endif
  read(IIN) vsv_array(:,:,:,:)
  close(IIN)

  ! read vsh
  write(filename, '(a,i6.6,a)') trim(model_dir) &
               //'/proc',myrank,'_'//trim(vsh_fn)//'.bin'
  open(IIN,file=trim(filename),status='old',form='unformatted',&
       action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) '  array not found: ',trim(filename)
    stop 'Error array file not found'
  endif
  read(IIN) vsh_array(:,:,:,:)
  close(IIN)

  ! read dvsv
  write(filename, '(a,i6.6,a)') trim(kernel_dir) &
               //'/proc',myrank,'_'//trim(dvsv_fn)//'.bin'
  open(IIN,file=trim(filename),status='old',form='unformatted',&
       action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) '  array not found: ',trim(filename)
    stop 'Error array file not found'
  endif
  read(IIN) dvsv_array(:,:,:,:)
  close(IIN)

  ! read dvsh
  write(filename, '(a,i6.6,a)') trim(kernel_dir) &
               //'/proc',myrank,'_'//trim(dvsh_fn)//'.bin'
  open(IIN,file=trim(filename),status='old',form='unformatted',&
       action='read',iostat=ier)
  if (ier /= 0) then
    write(*,*) '  array not found: ',trim(filename)
    stop 'Error array file not found'
  endif
  read(IIN) dvsh_array(:,:,:,:)
  close(IIN)

  vs_voigt_sq_array(:,:,:,:) = 2.0*vsv_array(:,:,:,:)*vsv_array(:,:,:,:)/3.0 +&
                       vsh_array(:,:,:,:) * vsh_array(:,:,:,:) / 3.0
  drho_array(:,:,:,:) = 2.0*vsv_array(:,:,:,:)*vsv_array(:,:,:,:)/3.0/&
                        vs_voigt_sq_array(:,:,:,:)*dvsv_array(:,:,:,:) +&
                        vsh_array(:,:,:,:)*vsh_array(:,:,:,:)/3.0/&
                        vs_voigt_sq_array(:,:,:,:)*dvsh_array(:,:,:,:)
  drho_array(:,:,:,:) = SCALE_RHO * drho_array(:,:,:,:)

  ! write drho
  write(filename, '(a,i6.6,a)') trim(kernel_dir) &
               //'/proc',myrank,'_'//trim(drho_fn)//'.bin'
  open(IOUT,file=trim(filename),status='unknown',form='unformatted',&
       action='write',iostat=ier)
  write(IOUT) drho_array
  close(IOUT)

  deallocate(vsv_array,vsh_array,vs_voigt_sq_array,dvsv_array,dvsh_array,&
             drho_array)
  call finalize_mpi()
  
end program scale_rho_with_vs
