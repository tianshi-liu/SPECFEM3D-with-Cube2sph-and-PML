program subtract_kernel_normalize
  use postprocess_par, only : MAX_KERNEL_NAMES
  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic
      !nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic
  integer, parameter :: NARGS = 4
  character(len=MAX_STRING_LEN) :: arg(NARGS)
  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: in_kernel_path1, in_kernel_path2, dir_out, &
                                   kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: filename
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: kernel_array1, &
                      kernel_array2, kernel_array
  integer :: i,ier,nker
  real(kind=CUSTOM_REAL) :: max_val
  logical :: BROADCAST_AFTER_READ
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
    write(*,*) 'subtracting kernel and normalize'
    write(*,*)
  endif
  call synchronize_all()

  do i = 1, NARGS
    call get_command_argument(i,arg(i),status=ier)
  enddo

  read(arg(1),'(a)') in_kernel_path1
  read(arg(2),'(a)') in_kernel_path2
  read(arg(3),'(a)') kernel_names_comma_delimited
  read(arg(4),'(a)') dir_out

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  
! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  call initialize_simulation()

  ! reads in external mesh
  call read_mesh_databases()

  allocate(kernel_array1(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), &
           kernel_array2(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), &
           kernel_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)

  call read_kernel_from_path(in_kernel_path1,kernel_names,nker,kernel_array1)
  if (myrank == 0) then
    print *, in_kernel_path1
  endif

  call read_kernel_from_path(in_kernel_path2,kernel_names,nker,kernel_array2)
  if (myrank == 0) then
    print *, in_kernel_path2
  endif

  kernel_array(:,:,:,:,:) = kernel_array1(:,:,:,:,:) - kernel_array2(:,:,:,:,:)

  call max_all_all_cr(maxval(abs(kernel_array)), max_val)

  if (myrank == 0) then
    print *, 'normalize with max absolute value = ', max_val, ' at rank ', myrank
  endif

  if (myrank == (sizeprocs - 1)) then
    print *, 'normalize with max absolute value = ', max_val, ' at rank ', myrank
  endif

  kernel_array(:,:,:,:,:) = kernel_array(:,:,:,:,:) / max_val

  do i = 1, nker
    write(filename,'(a,i6.6,a)') trim(dir_out) // &
                '/proc',myrank,'_'// trim(kernel_names(i))//'.bin'
    open(IOUT,file=trim(filename),form='unformatted',action='write',&
              status='unknown', iostat=ier)
    write(IOUT) kernel_array(:,:,:,:,i)
    close(IOUT)
  enddo
  deallocate(ibool,irregular_element_number)
  deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian)
  deallocate(xstore,ystore,zstore,kappastore,mustore)
  deallocate(ispec_is_acoustic, ispec_is_elastic, ispec_is_poroelastic)
  deallocate(kernel_array1, kernel_array2, kernel_array)
  call finalize_mpi()
  
end program subtract_kernel_normalize


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
