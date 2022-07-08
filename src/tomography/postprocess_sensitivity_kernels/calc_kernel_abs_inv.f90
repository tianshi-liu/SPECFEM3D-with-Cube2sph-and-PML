program calc_kernel_abs_inv
!use tomography_par, only: MAX_STRING_LEN, MAX_KERNEL_PATHS, IIN, myrank, &
!    sizeprocs, NGLOB, NSPEC, USE_ALPHA_BETA_RHO, USE_ISO_KERNELS
use postprocess_par, only: MAX_STRING_LEN,MAX_KERNEL_PATHS,IIN, IOUT,&
                     myrank,sizeprocs,NGLOB,NSPEC,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
use shared_parameters

implicit none
integer, parameter :: NARGS = 3
character(len=MAX_STRING_LEN) :: kernel_path, &
                 kernel_names_comma_delimited, kernel_names(MAX_KERNEL_PATHS)
character(len=MAX_STRING_LEN) :: prname_lp,output_dir,input_dir
character(len=MAX_STRING_LEN) :: arg(NARGS)
integer :: nker
integer :: i,ier,iker

logical :: BROADCAST_AFTER_READ
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kernel

call init_mpi()
call world_size(sizeprocs)
call world_rank(myrank)
call synchronize_all()
! check command line arguments
if (command_argument_count() /= NARGS) then
  if (myrank == 0) then
    print *, 'USAGE: mpirun -np NPROC bin/calc_kernel_abs KERNEL_NAMES INPUT_DIR OUTPUT_DIR'
    stop ' Please check command line arguments'
  endif
endif
call synchronize_all()
do i = 1, NARGS
  call get_command_argument(i,arg(i), status=ier)
enddo

read(arg(1),'(a)') kernel_names_comma_delimited
read(arg(2),'(a)') input_dir
read(arg(3),'(a)') output_dir

! parse names from KERNEL_NAMES
call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)

BROADCAST_AFTER_READ = .true.
call read_parameter_file(myrank,BROADCAST_AFTER_READ)
write(prname_lp,'(a,i6.6,a)') &
    trim(LOCAL_PATH)//'/proc',myrank,'_'//'external_mesh.bin'
open(unit=IIN, file=trim(prname_lp), &
    status='old',action='read',form='unformatted',iostat=ier)
read(IIN) NSPEC
read(IIN) NGLOB
close(IIN)
allocate(kernel(NGLLX, NGLLY, NGLLZ, NSPEC), stat=ier)
!allocate(hess(NGLLX, NGLLY, NGLLZ, NSPEC), stat=ier)
!kernel = 0.0
!hess = 0.0
!write(hess_path, '(a,i6.6,a)') trim(input_dir)//'/proc', myrank, &
!    '_'//trim(kernel_names(nker))//'.bin'
!open(IIN, file=trim(hess_path), status='old',form='unformatted',&
!             action='read',iostat=ier)
!read(IIN) hess
!close(IIN)
!hess = abs(hess)
!call invert_hess(hess)
do iker = 1, nker
  kernel = 0.0
  write(kernel_path, '(a,i6.6,a)') trim(input_dir)//'/proc', myrank, &
    '_'//trim(kernel_names(iker))//'.bin'
  open(IIN, file=trim(kernel_path), status='old',form='unformatted',&
             action='read',iostat=ier)
  read(IIN) kernel
  close(IIN)
  kernel = abs(kernel)
  call invert_hess(kernel)
  !kernel = kernel * hess
  write(kernel_path, '(a,i6.6,a)') trim(output_dir)//'/proc', myrank, &
    '_'//trim(kernel_names(iker))//'_abs_inv.bin'
  open(IOUT, file=trim(kernel_path), status='unknown',form='unformatted',&
             action='write',iostat=ier)
  write(IOUT) abs(kernel)
  close(IOUT)
enddo
deallocate(kernel)
call finalize_mpi()
end program calc_kernel_abs_inv


subroutine invert_hess( hess_matrix )

! inverts the Hessian matrix
! the approximate Hessian is only defined for diagonal elements: like
! H_nn = \frac{ \partial^2 \chi }{ \partial \rho_n \partial \rho_n }
! on all GLL points, which are indexed (i,j,k,ispec)

  use tomography_par, only: THRESHOLD_HESS
  use postprocess_par, only: myrank,NSPEC,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: hess_matrix

  ! local parameters
  real(kind=CUSTOM_REAL) :: maxh,maxh_all

  ! maximum value of Hessian
  maxh = maxval( abs(hess_matrix) )

  ! determines maximum from all slices on master
  call max_all_all_cr(maxh, maxh_all)

  ! user output
  if (myrank == 0) then
    !print *
    print *,'Hessian maximum: ',maxh_all
    !print *
  endif

  ! normalizes Hessian
  if (maxh_all < 1.e-18) then
    ! Hessian is zero, re-initializes
    hess_matrix = 1.0_CUSTOM_REAL
    !stop 'Error Hessian too small'
  else
    ! since Hessian has absolute values, this scales between [0,1]
    hess_matrix = hess_matrix / maxh_all
  endif

  !!!!!!TL: change to adding a small value
  !! inverts Hessian values
  !where( abs(hess_matrix(:,:,:,:)) > THRESHOLD_HESS )
  !  hess_matrix = 1.0_CUSTOM_REAL / hess_matrix
  !elsewhere
  !  hess_matrix = 1.0_CUSTOM_REAL / THRESHOLD_HESS
  !endwhere
  hess_matrix = 1.0_CUSTOM_REAL / (hess_matrix + THRESHOLD_HESS)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! rescales Hessian
  !hess_matrix = hess_matrix * maxh_all

end subroutine invert_hess
