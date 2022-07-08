! This program is used to compute lbfgs update direction
! Author: Kai Wang, wangkai@physics.utoronto.ca
! University of Toronto, Ontario, CA
! Last modified: Dec 23, 2018
! 
! 2018/12/23, Kai, revised from compute_direction_lbfgs.f90
!     written by Hejun Zhu
! 2020/11/09, modified by Tianshi Liu
! -- proper definition of inner product
! -- more flexible file names and kernel names
! -- compute kernels as local arrays
!=====================================================================

  program write_lbfgs_direction_smooth

! calculates TI gradient based on a conjugate gradient method
!
! based on: Tarantola, Inverse problem theory, 2005.
!                  section 6.22.7 conjugate directions, page 217.
!                  formula for alpha_n based on Polak & Ribiere (1969)
!
! note: we use a preconditioner F_0 = 1, thus lambda_n = gamma_n in (6.322)
!          and use gamma_n as the smoothed kernel (for bulk_c, bulk_betav,..).
!
!          however, one could see smoothing as preconditioner F_0, thus
!          gamma_n would be un-smoothed kernel and lambda_n would be smoothed one...
!          i'm not sure if this makes a difference.

  use postprocess_par, only : MAX_KERNEL_NAMES, MAX_KERNEL_PATHS
  !use shared_parameters
  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic
      !nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic
  !!!
  implicit none
  !!!!!!!! originally from tomography_par.f90 !!!!!!!!!!
  integer, parameter :: max_store = 15
  integer, parameter :: NARGS = 6
  real(CUSTOM_REAL), parameter :: max_angle = 85.0
  !integer:: iker
  !character(len=MAX_STRING_LEN) :: filename,dirname

  character(len=MAX_STRING_LEN) :: grad_names(MAX_KERNEL_NAMES), &
                                   kernel_names(MAX_KERNEL_NAMES), &
                                   kernel0_paths(MAX_KERNEL_PATHS), &
                                   kernel1_paths(MAX_KERNEL_PATHS), &
                                   model_update_names(MAX_KERNEL_NAMES), &
                                   model_update_paths(MAX_KERNEL_PATHS), &
                                   out_kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: grad_names_comma_delimited, &
                                   kernel_names_comma_delimited, &
                                   model_update_names_comma_delimited, &
                                   out_kernel_names_comma_delimited, &
                                   out_kernel_path,&
                                   out_kernel_name,filename,&
                                   current_kernel_path, fn_store
  !logical :: is_perturbation(MAX_KERNEL_NAMES)

  character(len=MAX_STRING_LEN) :: arg(NARGS)

  integer :: istore_used(max_store)
  integer :: nstore_used
  real(kind=CUSTOM_REAL), dimension(MAX_KERNEL_PATHS) :: step_len
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!! Kai added !!!
  ! L-BFGS
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: &
          q_vector,r_vector,gradient1,gradient0,&
          gradient_diff,model_diff
  real(kind=CUSTOM_REAL),dimension(max_store) :: p,a
  real(kind=CUSTOM_REAL) :: p_tmp,p_sum,a_tmp,a_sum,b_tmp,b_sum
  real(kind=CUSTOM_REAL) :: g_tmp,g_sum,m_tmp,m_sum,min_val,max_val
  real(kind=CUSTOM_REAL) :: b,p_k_up,p_k_down,p_k_up_sum,p_k_down_sum,p_k
  !integer :: iglob
  integer :: i,ier,nker,nmod,iker,itemp
  integer :: istore,nstore_provided
  
  logical :: BROADCAST_AFTER_READ

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
    write(*,*) 'computing L-BFGS action'
    write(*,*)  
  endif 
  call synchronize_all()

  do i = 1, NARGS
    call get_command_argument(i,arg(i),status=ier)
  enddo
  read(arg(1),'(a)') grad_names_comma_delimited
  read(arg(2),'(a)') kernel_names_comma_delimited 
  read(arg(3),'(a)') model_update_names_comma_delimited
  read(arg(4),'(a)') out_kernel_names_comma_delimited
  !read(arg(5),'(a)') is_perturbation_comma_delimited
  read(arg(5),'(a)') out_kernel_path
  read(arg(6),'(a)') fn_store ! file that contain paths of current gradient and
                              ! all the s-y pairs

  call parse_kernel_names(grad_names_comma_delimited,grad_names,nker)
  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  call parse_kernel_names(model_update_names_comma_delimited,&
                          model_update_names,nmod)
  call parse_kernel_names(out_kernel_names_comma_delimited,out_kernel_names,nker)
  !call parse_kernel_names(is_perturbation_comma_delimited,is_perturbation_str,nker)
  !do iker = 1, nker
  !  read(is_perturbation_str(iker), *) is_perturbation(iker)
  !enddo
  
  if (nmod /= nker) stop 'number of kernel fields must match number of model fields'
  
  open(unit=IIN, file=trim(fn_store), action='read', form='formatted', &
       status='old', iostat=ier)
  read(IIN,'(a)') current_kernel_path
  read(IIN, *) nstore_provided
  do istore = 1, nstore_provided
    read(IIN,'(a)') kernel0_paths(istore)
    read(IIN,'(a)') kernel1_paths(istore)
    read(IIN,'(a)') model_update_paths(istore)
    read(IIN,*) step_len(istore)
  enddo
  
! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

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

  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),irregular_element_number(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 980')

  if (NSPEC_IRREGULAR > 0) then
    ! allocate arrays for storing the databases
    allocate(xix(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 981')
    allocate(xiy(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 982')
    allocate(xiz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 983')
    allocate(etax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 984')
    allocate(etay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 985')
    allocate(etaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 986')
    allocate(gammax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 987')
    allocate(gammay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 988')
    allocate(gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 989')
    allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 990')
   else
       ! allocate arrays for storing the databases
    allocate(xix(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 991')
    allocate(xiy(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 992')
    allocate(xiz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 993')
    allocate(etax(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 994')
    allocate(etay(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 995')
    allocate(etaz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 996')
    allocate(gammax(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 997')
    allocate(gammay(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 998')
    allocate(gammaz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 999')
    allocate(jacobian(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1000')
  endif
  ! mesh node locations
  allocate(xstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1001')
  allocate(ystore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1002')
  allocate(zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1003')
  if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1004')
  allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1005')
  if (ier /= 0) stop 'Error allocating arrays for material properties'

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1006')
  allocate(ispec_is_elastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1007')
  allocate(ispec_is_poroelastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1008')
  if (ier /= 0) stop 'Error allocating arrays for material flags'
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.


  ! reads in external mesh
  call read_mesh_databases()

  allocate(gradient1(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  allocate(gradient0(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  !allocate(model1(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  !allocate(model0(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  allocate(gradient_diff(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  allocate(model_diff(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  allocate(q_vector(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  allocate(r_vector(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  ! initialize arrays
  a(:)=0.0
  p(:)=0.0
  gradient1(:,:,:,:,:)=0.0
  gradient0(:,:,:,:,:)=0.0
  !model1(:,:,:,:,:)=0.0
  !model0(:,:,:,:,:)=0.0
  gradient_diff(:,:,:,:,:)=0.0
  model_diff(:,:,:,:,:)=0.0
  q_vector(:,:,:,:,:)=0.0
  r_vector(:,:,:,:,:)=0.0
  nstore_used = 0
  istore_used(:) = 0

  call read_kernel_from_path(current_kernel_path,grad_names,nker,q_vector)
  if (myrank == 0) then
    print *, current_kernel_path
  endif
  
  if (myrank == 0) then
    print *,'************************************************'
    print *,'*******starting backward store *****************'
    print *,'************************************************'
  endif

  do istore=nstore_provided,1,-1
    if (nstore_used >= max_store) then
      if (myrank == 0) print *, 'max number of s-y pairs is achieved'
      exit
    endif
    call read_kernel_from_path(kernel0_paths(istore),kernel_names,&
                               nker,gradient0)
    call read_kernel_from_path(kernel1_paths(istore),kernel_names,&
                               nker,gradient1)
    call read_kernel_from_path(model_update_paths(istore),model_update_names,&
                               nmod,model_diff)
    model_diff(:,:,:,:,:) = model_diff(:,:,:,:,:) * step_len(istore)
    !call read_kernel_from_path(model1_paths(istore),model_names,&
    !                           nmod,model1)
    !do iker = 1, nker
    !  if (is_perturbation(iker)) &
    !    model0(:,:,:,:,iker) = log(model0(:,:,:,:,iker))
    !    model1(:,:,:,:,iker) = log(model1(:,:,:,:,iker))
    !enddo
    gradient_diff=gradient1-gradient0
    !model_diff=model1-model0
    if (myrank == 0) then
      print *, kernel0_paths(istore)
      print *, kernel1_paths(istore)
      print *, model_update_paths(istore)
      print *, step_len(istore)
    endif

    call calc_inner_product(gradient_diff,model_diff,nker,p_tmp)
    call sum_all_all_cr(p_tmp,p_sum)

    call calc_inner_product(gradient_diff,gradient_diff,nker,g_tmp)
    call sum_all_all_cr(g_tmp,g_sum)

    call calc_inner_product(model_diff,model_diff,nker,m_tmp)
    call sum_all_all_cr(m_tmp,m_sum)

    if ((abs(p_sum)/(sqrt(g_sum)*sqrt(m_sum)) < cos(max_angle/180.0*PI)) & 
        .or. (abs(p_sum) < 1.e-22)) then
      ! the gradient and model diff are almost orthogonal,
      ! this s-y pair is unused, go to the next pair
      if (myrank == 0) &
        print *, 'angle between gradient and model diff is smaller than ', &
                  max_angle, 'degree, exclude pair', istore
      cycle 
    else
      ! use this s-y pair
      if (myrank == 0) print *, 'use pair', istore
      nstore_used = nstore_used + 1
      istore_used(nstore_used) = istore
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    p(nstore_used) = 1.0_CUSTOM_REAL / p_sum !1.0/p_sum

    call calc_inner_product(model_diff,q_vector,nker,a_tmp)
    call sum_all_all_cr(a_tmp,a_sum)

    a(nstore_used)=p(nstore_used)*a_sum    !a_i <--- p_i*s_i^T*q

    q_vector=q_vector-a(nstore_used)*gradient_diff  !q<--- q-a_i*y_i
  enddo
  if (nstore_used > 0) then
    istore = istore_used(1)
    call read_kernel_from_path(kernel0_paths(istore),kernel_names,&
                               nker,gradient0)
    call read_kernel_from_path(kernel1_paths(istore),kernel_names,&
                               nker,gradient1)
    call read_kernel_from_path(model_update_paths(istore),model_update_names,&
                               nmod,model_diff)
    model_diff(:,:,:,:,:) = model_diff(:,:,:,:,:) * step_len(istore)
    !model0(:,:,:,:,:) = log(model0(:,:,:,:,:))
    !call read_kernel_from_path(model1_paths(istore),model_names,&
    !                           nmod,model1)
    !model1(:,:,:,:,:) = log(model1(:,:,:,:,:))
    !do iker = 1, nker
    !  if (is_perturbation(iker)) &
    !    model0(:,:,:,:,iker) = log(model0(:,:,:,:,iker))
    !    model1(:,:,:,:,iker) = log(model1(:,:,:,:,iker))
    !enddo
    gradient_diff=gradient1-gradient0
    !model_diff=model1-miodel0

    if (myrank == 0) then
      print *, 'estimate scaling using:'
      print *, kernel0_paths(istore)
      print *, kernel1_paths(istore)
      print *, model_update_paths(istore)
      print *, step_len(istore)
    endif

! this implements Algorithm equation (9.6) on page 226 of the book of
! Jorge Nocedal and Stephen Wright, "Numerical Optimization", Springer (2006)
    call calc_inner_product(gradient_diff, model_diff, nker, p_k_up)
    call sum_all_all_cr(p_k_up, p_k_up_sum)

    call calc_inner_product(gradient_diff, gradient_diff, nker, p_k_down)
    call sum_all_all_cr(p_k_down, p_k_down_sum)
    
    p_k=p_k_up_sum/p_k_down_sum
    r_vector=p_k*q_vector
  else
    r_vector=1.0*q_vector
  endif  !nstore_used > 0
  if (myrank == 0) then
     print *,'******************************************'
     print *,'********starting forward store ***********'
     print *,'******************************************'
  endif

  do itemp = nstore_used,1,-1
    istore = istore_used(itemp)
    call read_kernel_from_path(kernel0_paths(istore),kernel_names,&
                               nker,gradient0)
    call read_kernel_from_path(kernel1_paths(istore),kernel_names,&
                               nker,gradient1)
    call read_kernel_from_path(model_update_paths(istore),model_update_names,&
                               nmod,model_diff)
    model_diff(:,:,:,:,:) = model_diff(:,:,:,:,:) * step_len(istore)
    !model0(:,:,:,:,:) = log(model0(:,:,:,:,:))
    !call read_kernel_from_path(model1_paths(istore),model_names,&
    !                           nmod,model1)
    !model1(:,:,:,:,:) = log(model1(:,:,:,:,:))
    !do iker = 1, nker
    !  if (is_perturbation(iker)) &
    !    model0(:,:,:,:,iker) = log(model0(:,:,:,:,iker))
    !    model1(:,:,:,:,iker) = log(model1(:,:,:,:,iker))
    !enddo
     
    gradient_diff=gradient1-gradient0
    !model_diff=model1-model0

    if (myrank == 0) then
      print *, kernel0_paths(istore)
      print *, kernel1_paths(istore)
      print *, model_update_paths(istore)
      print *, step_len(istore)
    endif

    call calc_inner_product(gradient_diff, r_vector, nker, b_tmp)
    call sum_all_all_cr(b_tmp, b_sum) 
    b=p(itemp)*b_sum

    r_vector=r_vector+model_diff*(a(itemp)-b)

  enddo

  r_vector=-1.0*r_vector
  
  do iker = 1, nker
    out_kernel_name = out_kernel_names(iker)
    write(filename,'(a,i6.6,a)') trim(out_kernel_path) //'/proc',myrank,'_'//&
                        trim(out_kernel_name)//'.bin'
    open(IOUT,file=trim(filename),form='unformatted',action='write',iostat=ier)
    write(IOUT) r_vector(:,:,:,:,iker)
    close(IOUT)
  enddo

  do iker = 1, nker
    call min_all_cr(minval(r_vector(:,:,:,:,iker)),min_val)
    call max_all_cr(maxval(r_vector(:,:,:,:,iker)),max_val)
    if (myrank == 0) then
      print *, trim(out_kernel_names(iker))
      print *, '    min:', min_val
      print *, '    max:', max_val
    endif
  enddo
  
  !call max_all_cr(maxval(abs(r_vector)), max_val)
  call read_kernel_from_path(current_kernel_path,kernel_names,nker,q_vector)
  call calc_inner_product(r_vector,q_vector,nker,p_tmp)
  call sum_all_all_cr(p_tmp,p_sum)

  if (myrank == 0) then
    open(unit=IOUT, file='line_search_derivative', form='formatted', &
         action='write', status='unknown', iostat=ier)
    !write(IOUT, *) max_val
    write(IOUT, *) p_sum
    close(IOUT)
  endif
  deallocate(ibool,irregular_element_number)
  deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian)
  deallocate(xstore,ystore,zstore,kappastore,mustore)
  deallocate(ispec_is_acoustic, ispec_is_elastic, ispec_is_poroelastic)
  deallocate(q_vector,r_vector,gradient1,gradient0,&
             gradient_diff,model_diff)
  call finalize_mpi()
  

!********************************************************

  !call synchronize_all()

  end program write_lbfgs_direction_smooth


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

subroutine calc_inner_product(vect1, vect2, nv, q)
  use specfem_par, only: jacobian, irregular_element_number, &
                         wxgll, wygll, wzgll, NGLLX, NGLLY, NGLLZ, NSPEC_AB, &
                         CUSTOM_REAL
  integer, intent(in) :: nv
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nv), &
                          intent(in):: vect1, vect2
  real(kind=CUSTOM_REAL), intent(out) :: q
  ! local variables
  integer :: iv, i,j,k,ispec,ispec_irreg
  real(kind=CUSTOM_REAL) :: jacobianl, weight, coeff_n1, coeff_n2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nv) :: &
                          vect1n, vect2n

  ! nornalize
  coeff_n1 = maxval(abs(vect1(:,:,:,:,:)))
  if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL
  vect1n(:,:,:,:,:) = vect1(:,:,:,:,:) / coeff_n1

  coeff_n2 = maxval(abs(vect2(:,:,:,:,:)))
  if (coeff_n2 == 0._CUSTOM_REAL) coeff_n2=1._CUSTOM_REAL
  vect2n(:,:,:,:,:) = vect2(:,:,:,:,:) / coeff_n2

  q = 0._CUSTOM_REAL

  do iv = 1, nv
    do ispec = 1, NSPEC_AB
      ispec_irreg = irregular_element_number(ispec)
      do k=1, NGLLZ; do j=1,NGLLY; do i=1,NGLLX
        weight = wxgll(i)*wygll(j)*wzgll(k)
        if (ispec_irreg /= 0) jacobianl = jacobian(i,j,k,ispec_irreg)
        q = q + jacobianl * weight * vect1n(i,j,k,ispec,iv) * &
                vect2n(i,j,k,ispec,iv)
      enddo;enddo;enddo
    enddo
  enddo
  q = q * coeff_n1 * coeff_n2

end subroutine calc_inner_product
