program calc_kernel_mul_inv_hess
!use tomography_par, only: MAX_STRING_LEN, MAX_KERNEL_PATHS, IIN, myrank, &
!    sizeprocs, NGLOB, NSPEC, USE_ALPHA_BETA_RHO, USE_ISO_KERNELS
use postprocess_par, only: MAX_STRING_LEN,MAX_KERNEL_PATHS,IIN, IOUT,&
                     myrank,sizeprocs,NGLOB,NSPEC,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
use shared_parameters

implicit none
integer, parameter :: NARGS = 3
character(len=MAX_STRING_LEN) :: kernel_path, hess_path, &
                 kernel_names_comma_delimited, kernel_names(MAX_KERNEL_PATHS)
character(len=MAX_STRING_LEN) :: prname_lp,output_dir,input_dir
character(len=MAX_STRING_LEN) :: arg(NARGS)
integer :: nker
integer :: i,ier,iker

logical :: BROADCAST_AFTER_READ
real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kernel, hess

call init_mpi()
call world_size(sizeprocs)
call world_rank(myrank)
call synchronize_all()
! check command line arguments
if (command_argument_count() /= NARGS) then
  if (myrank == 0) then
    print *, 'USAGE: mpirun -np NPROC bin/calc_kernel_div_hess KERNEL_NAMES INPUT_DIR OUTPUT_DIR'
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
allocate(hess(NGLLX, NGLLY, NGLLZ, NSPEC), stat=ier)
!kernel = 0.0
hess = 0.0
write(hess_path, '(a,i6.6,a)') trim(input_dir)//'/proc', myrank, &
    '_'//trim(kernel_names(nker))//'.bin'
open(IIN, file=trim(hess_path), status='old',form='unformatted',&
             action='read',iostat=ier)
read(IIN) hess
close(IIN)
!hess = abs(hess)
!call invert_hess(hess)
do iker = 1, nker - 1
  kernel = 0.0
  write(kernel_path, '(a,i6.6,a)') trim(input_dir)//'/proc', myrank, &
    '_'//trim(kernel_names(iker))//'.bin'
  open(IIN, file=trim(kernel_path), status='old',form='unformatted',&
             action='read',iostat=ier)
  read(IIN) kernel
  close(IIN)
  kernel = kernel * hess
  write(kernel_path, '(a,i6.6,a)') trim(output_dir)//'/proc', myrank, &
    '_'//trim(kernel_names(iker))//'_precon.bin'
  open(IOUT, file=trim(kernel_path), status='unknown',form='unformatted',&
             action='write',iostat=ier)
  write(IOUT) kernel
  close(IOUT)
enddo
deallocate(hess, kernel)
call finalize_mpi()
end program calc_kernel_mul_inv_hess
