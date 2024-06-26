program get_max_absolute_value
  use postprocess_par, only : MAX_KERNEL_NAMES
  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic
      !nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic
  integer, parameter :: NARGS = 3
  character(len=MAX_STRING_LEN) :: arg(NARGS)
  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: in_kernel_path, out_fn, &
                                   kernel_names_comma_delimited
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: kernel_array
  integer :: i,ier,nker
  real(kind=CUSTOM_REAL) :: max_val
  logical :: BROADCAST_AFTER_READ
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
    write(*,*) 'getting maximum absolute value'
    write(*,*)
  endif
  call synchronize_all()

  do i = 1, NARGS
    call get_command_argument(i,arg(i),status=ier)
  enddo

  read(arg(1),'(a)') in_kernel_path
  read(arg(2),'(a)') kernel_names_comma_delimited
  read(arg(3),'(a)') out_fn
  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  
! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  call initialize_simulation()

  ! reads in external mesh
  call read_mesh_databases()

  allocate(kernel_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)

  call read_kernel_from_path(in_kernel_path,kernel_names,nker,kernel_array)
  if (myrank == 0) then
    print *, in_kernel_path
  endif

  call max_all_cr(maxval(abs(kernel_array)), max_val)

  if (myrank == 0) then
    open(unit=IOUT, file=trim(out_fn), form='formatted', &
         action='write', status='unknown', iostat=ier)
    write(IOUT, *) max_val
    close(IOUT)
  endif
  deallocate(ibool,irregular_element_number)
  deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian)
  deallocate(xstore,ystore,zstore,kappastore,mustore)
  deallocate(ispec_is_acoustic, ispec_is_elastic, ispec_is_poroelastic)
  deallocate(kernel_array)
  call finalize_mpi()
  
end program get_max_absolute_value


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
