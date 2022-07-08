program write_lbfgs_direction_iso
  use tomography_par
  use tomography_kernels_iso
  !use specfem_par, only: jacobian,jacobian_regular,irregular_element_number, &
  !    wxgll, wygll, wzgll
  integer :: ier
  !!!
  ! allocate arrays for storing gradient
  call initialize()
  allocate(model_dbulk(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_dbeta(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_drho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)

  if (ier /= 0) stop 'error allocating gradient arrays'

  ! initializes arrays
  model_dbulk = 0.0_CUSTOM_REAL
  model_dbeta = 0.0_CUSTOM_REAL
  model_drho = 0.0_CUSTOM_REAL
  call read_parameters_lbfgs()
  call get_lbfgs_direction_iso()
  call write_gradient_iso()

  call synchronize_all()
  call finalize_mpi()
end program write_lbfgs_direction_iso

subroutine initialize()

! initializes arrays
! TL: modify it to accomodate for external mesher

  use tomography_par

  use specfem_par, only: NPROC,ADIOS_ENABLED,LOCAL_PATH

  implicit none

  logical :: BROADCAST_AFTER_READ
  character(len=MAX_STRING_LEN) :: database_name, prname
  integer :: ier

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! reads the parameter file
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  if (ADIOS_ENABLED) stop 'Flag ADIOS_ENABLED not supported yet for xadd_model, please rerun program...'

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *, 'Error number of processors supposed to run on: ',NPROC
      print *, 'Error number of MPI processors actually run on: ',sizeprocs
      print *
      print *, 'please rerun with: mpirun -np ',NPROC,' bin/xadd_model .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif
  call create_name_database(prname,myrank,LOCAL_PATH)
  database_name = prname(1:len_trim(prname))//'external_mesh.bin'
  open(unit=IIN,file=trim(database_name),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error could not open database file: ',trim(database_name)
    call exit_mpi(myrank,'Error opening database file')
  endif
  read(IIN) NSPEC
  read(IIN) NGLOB
  close(IIN)

  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some
  ! array sizes below
  !call read_mesh_for_init()

  ! sets tomography array dimensions
  !NSPEC = NSPEC_AB
  !NGLOB = NGLOB_AB

end subroutine initialize
