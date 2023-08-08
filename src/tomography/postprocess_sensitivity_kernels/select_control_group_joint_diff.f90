program select_control_group
  use postprocess_par, only : MAX_KERNEL_NAMES, MAX_KERNEL_PATHS
  !use shared_parameters
  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic
      !nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic

  implicit none
  integer, parameter :: NARGS = 6
  integer, parameter :: MAX_DATASETS = 2

  character(len=MAX_STRING_LEN) :: &
          kernel_paths(MAX_KERNEL_PATHS, MAX_DATASETS), &
          meas_paths(MAX_KERNEL_PATHS, MAX_DATASETS), &
          kernel_names(MAX_KERNEL_NAMES), &
          kernel_names_comma_delimited, &
          dir_lists_comma_delimited, & 
          meas_lists_comma_delimited
  logical :: kept(MAX_KERNEL_PATHS, MAX_DATASETS)
  character(len=MAX_STRING_LEN) :: sline,output_dir, &
                                   fns_dir_list(MAX_KERNEL_NAMES), &
                                   fns_meas_list(MAX_KERNEL_NAMES),&
                                   kernel_name,filename
  character(len=MAX_STRING_LEN) :: arg(NARGS)
  integer :: npath(MAX_DATASETS),nker,min_num_path,npath_left(MAX_DATASETS),&
             nset, nmeas(MAX_KERNEL_PATHS, MAX_DATASETS), &
             nmeas_all(MAX_DATASETS), nmeas_left(MAX_DATASETS), &
             ndummy
  integer :: i,ier,iker,ipath,iset, ipath_min,iset_min
  !integer :: list_path(MAX_KERNEL_PATHS)
  real(kind=CUSTOM_REAL) :: threshold_diff, diff, norm_sum_all,&
                           norm_sum_sub, norm_sum_all_all, &
                           norm_sum_sub_all, min_diff
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: &
                                   sum_all_array, array, sum_ctrlgrp_array
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:),allocatable :: &
                                   sum_sub_array
  !real(kind=CUSTOM_REAL), parameter :: SMALL_VAL = 1.0e-20
  double precision :: chi_dummy
  logical :: BROADCAST_AFTER_READ

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
    write(*,*) 'selecting control group'
    write(*,*)
  endif
  call synchronize_all()

  ! check command line arguments
  !if (command_argument_count() /= NARGS) then
  !  if (myrank == 0) then
  !    print *, 'USAGE: mpirun -np NPROC bin/xselect_control_group KERNEL_NAMES INPUT_FILE OUTPUT_DIR'
  !    stop ' Please check command line arguments'
  !  endif
  !endif
  !call synchronize_all()

  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo

  read(arg(1),'(a)') kernel_names_comma_delimited
  read(arg(2),'(a)') dir_lists_comma_delimited
  read(arg(3),'(a)') meas_lists_comma_delimited
  read(arg(4),'(a)') output_dir
  read(arg(5),*) min_num_path
  read(arg(6),*) threshold_diff
  !threshold_cos = cos(threshold_cos*PI/180.0) ! get cos value

  ! parse names from KERNEL_NAMES
  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  call parse_kernel_names(dir_lists_comma_delimited,fns_dir_list,nset)
  call parse_kernel_names(meas_lists_comma_delimited,fns_meas_list,ndummy)
  if (myrank == 0) then
    if (.not. (nset == ndummy)) stop 'inconsist numbers of dir_list and meas_list'
  endif

  ! parse paths from INPUT_FILE
  do iset = 1, nset
  npath(iset)=0
  open(unit = IIN, file = trim(fns_dir_list(iset)), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(fns_dir_list(iset)), myrank
     stop 1
  endif
  do while (1 == 1)
     read(IIN,'(a)',iostat=ier) sline
     if (ier /= 0) exit
     npath(iset) = npath(iset)+1
     if (npath(iset) > MAX_KERNEL_PATHS) &
       stop 'Error number of paths exceeds MAX_KERNEL_PATHS'
     kernel_paths(npath(iset), iset) = sline
  enddo
  close(IIN)
  if (myrank == 0) then
    write(*,*) '  set ', iset, ' ' ,npath(iset),' events'
    write(*,*)
  endif
  ndummy = 0
  open(unit = IIN, file = trim(fns_meas_list(iset)), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(fns_meas_list(iset)), myrank
     stop 1
  endif
  do while (1 == 1)
     read(IIN,'(a)',iostat=ier) sline
     if (ier /= 0) exit
     ndummy = ndummy+1
     if (ndummy > MAX_KERNEL_PATHS) &
       stop 'Error number of paths exceeds MAX_KERNEL_PATHS'
     meas_paths(ndummy, iset) = sline
  enddo
  close(IIN)
  if (myrank == 0) then
    if (.not. (npath(iset) == ndummy)) &
      stop 'inconsist numbers of kernel paths and measurement paths'
  endif
  nmeas_all(iset) = 0
  do ipath = 1, npath(iset)
    filename = trim(meas_paths(ipath, iset)) // '/sum_chi'
    open(unit=IIN, file=trim(filename), form='formatted', status='old', &
         action='read', iostat=ier)
    read(IIN, *) chi_dummy
    read(IIN, *) nmeas(ipath, iset)
    close(IIN)
    if (myrank == 0) then
      print *, 'measurement path:', trim(meas_paths(ipath, iset))
      print *, '  nmeas = ', nmeas(ipath, iset)
    endif
    nmeas_all(iset) = nmeas_all(iset) + nmeas(ipath, iset)
  enddo
  if (myrank == 0) then
    print *, 'total number of measurements in dataset ', &
             iset, ': ', nmeas_all(iset)
  endif
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

  !allocate(sum_set_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker,nset),stat=ier)
  allocate(sum_all_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker),stat=ier)
  allocate(sum_ctrlgrp_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker),stat=ier)
  allocate(array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  allocate(sum_sub_array(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker,nset),stat=ier)

  sum_sub_array = 0._CUSTOM_REAL
  sum_all_array = 0._CUSTOM_REAL
  kept(:,:) = .true.
  do iset = 1, nset
  do ipath = 1, npath(iset)
    array = 0._CUSTOM_REAL
    call read_kernel_from_path(kernel_paths(ipath, iset),kernel_names,&
                               nker,array)
    sum_sub_array(:,:,:,:,:,iset) = sum_sub_array(:,:,:,:,:,iset) + &
                                    array(:,:,:,:,:)
  enddo
  sum_all_array(:,:,:,:,:) = sum_all_array(:,:,:,:,:) + &
                      sum_sub_array(:,:,:,:,:,iset) / (nmeas_all(iset) * 1.0)
  enddo
  
  nmeas_left(1:nset) = nmeas_all(1:nset)
  npath_left(1:nset) = npath(1:nset)
  !sum_sub_array(:,:,:,:,:,:) = sum_set_array(:,:,:,:,:,:)
  call calc_inner_product(sum_all_array,sum_all_array,nker,norm_sum_all)
  if (myrank == 0) norm_sum_all_all = 0.0
  call sum_all_cr(norm_sum_all,norm_sum_all_all)
  if (myrank == 0) norm_sum_all_all = sqrt(norm_sum_all_all)
  if (myrank == 0) print *, 'norm of total gradiant = ', norm_sum_all_all
  do while (any(npath_left(1:nset) > min_num_path))
    if (myrank == 0) then
      do iset = 1, nset
      print *, npath_left(iset), 'event kernels left in set ', iset
      enddo
      print *, 'attempting to remove one ...'
      min_diff = 100.0
      ipath_min = 0
      iset_min = 0
    endif
    do iset = 1, nset
    if (npath_left(iset) == min_num_path) cycle
    do ipath = 1, npath(iset)
      if (.not. kept(ipath, iset)) cycle
      array = 0._CUSTOM_REAL
      if (myrank == 0) print *, '    trying to remove event kernel', &
                       trim(kernel_paths(ipath, iset))
      call read_kernel_from_path(kernel_paths(ipath, iset),kernel_names,&
                                 nker,array)
      array(:,:,:,:,:) = (sum_sub_array(:,:,:,:,:,iset) - array(:,:,:,:,:))&
                       / ((nmeas_left(iset) - nmeas(ipath, iset)) * 1.0)
      do i = 1, nset
        if (i == iset) cycle
        array(:,:,:,:,:) = array(:,:,:,:,:) + &
                           sum_sub_array(:,:,:,:,:,i) / (nmeas_left(i) * 1.0)
      enddo
      array(:,:,:,:,:) = sum_all_array(:,:,:,:,:) - array(:,:,:,:,:)
      !call calc_inner_product(sum_all_array,array,nker,inner_prod)
      !if (myrank == 0) inner_prod_all = 0.0
      !call sum_all_cr(inner_prod,inner_prod_all)
      call calc_inner_product(array,array,nker,norm_sum_sub)
      if (myrank == 0) norm_sum_sub_all = 0.0
      call sum_all_cr(norm_sum_sub,norm_sum_sub_all)
      if (myrank == 0) then
        norm_sum_sub_all = sqrt(norm_sum_sub_all)
        diff = norm_sum_sub_all / norm_sum_all_all
        !cos_angle = inner_prod_all / (norm_sum_sub_all * norm_sum_all_all)
        !if (cos_angle > 1.0-SMALL_VAL) then
        !  print *, 'WARNING: cos out of range:', cos_angle
        !  cos_angle = 1.0-SMALL_VAL 
        !else if (cos_angle < -1.0+SMALL_VAL) then
        !  print *, 'WARNING: cos out of range:', cos_angle
        !  cos_angle = -1.0+SMALL_VAL
        !endif
        print *, '    diff = ', diff
        if (diff < min_diff) then
          min_diff = diff
          ipath_min = ipath
          iset_min = iset
        endif
      endif
    enddo
    enddo
    if (myrank == 0) then
      print *, 'minimum diff = ', diff
      if (min_diff > threshold_diff) then ! reach threshold
        print *, 'threshold has reached, no event kernel is removed'
        ipath_min = 0 
        iset_min = 0
      endif
    endif
    call bcast_all_singlei(ipath_min)
    call bcast_all_singlei(iset_min)
    if (ipath_min == 0) exit ! reach angle threshold, break directly
    if (myrank == 0) print *, 'remove event kernel ', &
                     trim(kernel_paths(ipath_min, iset_min))
    array = 0._CUSTOM_REAL
    call read_kernel_from_path(kernel_paths(ipath_min, iset_min),kernel_names,&
                               nker,array)
    sum_sub_array(:,:,:,:,:,iset_min) = sum_sub_array(:,:,:,:,:,iset_min) - &
                                  array(:,:,:,:,:)    
    nmeas_left(iset_min) = nmeas_left(iset_min) - nmeas(ipath_min, iset_min)
    npath_left(iset_min) = npath_left(iset_min) - 1
    kept(ipath_min, iset_min) = .false.
  enddo
  
  sum_ctrlgrp_array = 0._CUSTOM_REAL
  do iset = 1, nset
    sum_ctrlgrp_array(:,:,:,:,:) = sum_ctrlgrp_array(:,:,:,:,:) + &
                   sum_sub_array(:,:,:,:,:,iset) / (nmeas_left(iset) * 1.0)
  enddo
  call calc_inner_product(sum_ctrlgrp_array,sum_ctrlgrp_array,nker,norm_sum_sub)
  if (myrank == 0) norm_sum_sub_all = 0.0
  call sum_all_cr(norm_sum_sub,norm_sum_sub_all)
  if (myrank == 0) norm_sum_sub_all = sqrt(norm_sum_sub_all)
  if (myrank == 0) print *, 'norm of ctrlgrp gradiant = ', norm_sum_sub_all
  ! output kernels
  do iker = 1, nker
    kernel_name = kernel_names(iker)
    write(filename,'(a,i6.6,a)') trim(output_dir) //'/proc',myrank,'_'//&
                        trim(kernel_name)//'.bin'
    open(IOUT,file=trim(filename),form='unformatted',action='write',iostat=ier)
    write(IOUT) sum_all_array(:,:,:,:,iker)
    close(IOUT)
    write(filename,'(a,i6.6,a)') trim(output_dir) //'/proc',myrank,'_'//&
                        trim(kernel_name)//'_ctrlgrp.bin'
    open(IOUT,file=trim(filename),form='unformatted',action='write',&
              status='unknown', iostat=ier)
    write(IOUT) sum_ctrlgrp_array(:,:,:,:,iker)
    close(IOUT)
  enddo
  ! write components of control group
  if (myrank == 0) then
    do iset = 1, nset
    open(unit=IOUT, file=trim(fns_dir_list(iset))//'.ctrlgrp',form='formatted',&
                 action='write',status='unknown',iostat=ier)
    do ipath = 1, npath(iset)
      if (kept(ipath, iset)) &
        write(IOUT,'(a)') trim(kernel_paths(ipath, iset))
    enddo
    close(IOUT)
    enddo
  endif

  deallocate(ibool,irregular_element_number)
  deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian)
  deallocate(sum_all_array, array, sum_sub_array, &
             sum_ctrlgrp_array)
  deallocate(xstore,ystore,zstore,kappastore,mustore)
  deallocate(ispec_is_acoustic, ispec_is_elastic, ispec_is_poroelastic)

  call finalize_mpi()
  
end program select_control_group

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
