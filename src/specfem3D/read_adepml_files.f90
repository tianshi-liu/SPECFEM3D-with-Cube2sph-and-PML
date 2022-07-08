subroutine read_adepml_files()
  use pml_par
  use specfem_par
  use specfem_par_elastic

  implicit none
  integer :: ier, iout_vtk=998
  character(len=MAX_STRING_LEN) :: database_name, vtk_name
  integer :: ispec_CPML,i,j,k,ispec,iglob,iglob_CPML,iface, iside
  !integer, dimension(:), allocatable :: iglob_CPML_temp
  real, dimension(:,:,:,:), allocatable :: rvolume_loc
  logical :: POINT_EXIST

  database_name = prname(1:len_trim(prname))//'adepml_damping_indexing.bin'

  open(unit=27,file=trim(database_name),status='old', &
       action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error could not open database file: ',trim(database_name)
    call exit_mpi(myrank,'Error opening database file')
  endif

  read(27) NSPEC_CPML
  if (NSPEC_CPML > 0) then
    allocate(CPML_regions(NSPEC_CPML),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1505')
    if (ier /= 0) stop 'Error allocating array CPML_regions'
    allocate(CPML_to_spec(NSPEC_CPML),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1506')
    if (ier /= 0) stop 'Error allocating array CPML_to_spec'

    allocate(pml_d(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1507')
    if (ier /= 0) stop 'Error allocating array pml_d'
    allocate(pml_kappa(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1510')
    if (ier /= 0) stop 'Error allocating array pml_kappa'
    allocate(pml_beta(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1513')
    if (ier /= 0) stop 'Error allocating array pml_beta'
    allocate(ibool_CPML(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array ibool_CPML')
    if (ier /= 0) stop 'Error allocating array ibool_CPML'
    read(27) CPML_regions
    read(27) CPML_to_spec
    read(27) is_CPML
    read(27) pml_d
    read(27) pml_beta
    read(27) pml_kappa
    !read(27) nglob_CPML
    !read(27) ibool_CPML
    !allocate(CPML_to_glob(nglob_CPML))
    !read(27) CPML_to_glob
    !allocate(iglob_CPML_temp(NGLLX*NGLLY*NGLLZ*nspec_cpml))
    !nglob_CPML = 0
    !ibool_CPML(:,:,:,:) = 0
    !do ispec_CPML = 1, nspec_cpml
    !  do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
    !    ispec = CPML_to_spec(ispec_CPML)
    !    iglob = ibool(i,j,k,ispec)
    !    POINT_EXIST = .false.
    !    do iglob_CPML = 1, nglob_CPML
    !      if (iglob == iglob_CPML_temp(iglob_CPML)) then
    !        POINT_EXIST = .true.
    !        exit
    !      endif
    !    enddo
    !    if (.not. POINT_EXIST) then
    !      nglob_CPML = nglob_CPML + 1
    !      iglob_CPML_temp(nglob_CPML) = iglob
    !      ibool_CPML(i,j,k,ispec_CPML) = nglob_CPML
    !    else
    !      ibool_CPML(i,j,k,ispec_CPML) = iglob_CPML
    !    endif
    !  enddo;enddo;enddo
    !enddo
    !allocate(CPML_to_glob(nglob_CPML))
    !CPML_to_glob(1:nglob_CPML) = iglob_CPML_temp(1:nglob_CPML)
    !deallocate(iglob_CPML_temp)
  endif
  read(27) num_pml_physical
  if (num_pml_physical > 0) then
    allocate(pml_physical_ispec(num_pml_physical),stat=ier)
    allocate(pml_physical_ijk(NDIM, NGLLSQUARE, num_pml_physical),stat=ier)
    read(27) pml_physical_ispec
    read(27) pml_physical_ijk
  endif
  !read(27) num_interfaces_PML
  !if (num_interfaces_PML > 0) then
  !  read(27) max_nibool_interfaces_PML
  !  allocate(nibool_interfaces_PML(num_interfaces_PML),stat=ier)
  !  allocate(ibool_interfaces_PML(max_nibool_interfaces_PML,&
  !                                    num_interfaces_PML),stat=ier)
  !  allocate(my_neighbors_PML(num_interfaces_PML),stat=ier)
  !  read(27) nibool_interfaces_PML
  !  read(27) ibool_interfaces_PML
  !  read(27) my_neighbors_PML
  !endif
  close(27)
  
  database_name = prname(1:len_trim(prname))//'adepml_glob_indexing.bin'
  open(unit=27,file=trim(database_name),status='old', &
       action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error could not open database file: ',trim(database_name)
    call exit_mpi(myrank,'Error opening database file')
  endif
  read(27) nglob_CPML
  if (nglob_CPML > 0) then
    read(27) ibool_CPML
    allocate(CPML_to_glob(nglob_CPML))
    read(27) CPML_to_glob
  endif
  read(27) num_interfaces_PML
  if (num_interfaces_PML > 0) then
    read(27) max_nibool_interfaces_PML
    allocate(nibool_interfaces_PML(num_interfaces_PML),stat=ier)
    allocate(ibool_interfaces_PML(max_nibool_interfaces_PML,&
                                      num_interfaces_PML),stat=ier)
    allocate(my_neighbors_PML(num_interfaces_PML),stat=ier)
    read(27) nibool_interfaces_PML
    read(27) ibool_interfaces_PML
    read(27) my_neighbors_PML
  endif
  read(27) nglob_pml_in
  if (nglob_pml_in > 0) then
    allocate(pml_in_iglob(nglob_pml_in))
    read(27) pml_in_iglob
  endif
  read(27) nglob_dirichlet
  if (nglob_dirichlet > 0) then
    allocate(iglob_dirichlet(nglob_dirichlet))
    read(27) iglob_dirichlet
  endif
  close(27)  

  if (NSPEC_CPML > 0) then
    database_name = prname(1:len_trim(prname))//'adepml_param.bin'

    open(unit=27,file=trim(database_name),status='old', &
       action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error could not open database file: ',trim(database_name)
      call exit_mpi(myrank,'Error opening database file')
    endif
    if (num_pml_physical > 0) then
      allocate(pml_physical_normal(NDIM,NGLLSQUARE,num_pml_physical),stat=ier)
      allocate(pml_physical_jacobian2Dw(NGLLSQUARE,num_pml_physical),stat=ier)
      read(27) pml_physical_normal
      read(27) pml_physical_jacobian2Dw
    endif
    allocate(r_trans(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    allocate(r_trans_inv(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    allocate(rvolume(nglob_CPML),stat=ier)
    allocate(rvolume_loc(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    read(27) rvolume_loc
    rvolume(:) = 0.0
    do ispec_cpml = 1, NSPEC_CPML
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            rvolume(ibool_CPML(i,j,k,ispec_cpml)) = & 
                    rvolume(ibool_CPML(i,j,k,ispec_cpml)) + &
                    rvolume_loc(i,j,k,ispec_cpml)
          enddo
        enddo
      enddo
    enddo
    deallocate(rvolume_loc)
    read(27) r_trans
    read(27) r_trans_inv
    close(27)
    allocate(pml_spec_physical(2*NDIM,NSPEC_CPML))
    pml_spec_physical(:,:) = 0
    do iface = 1, num_pml_physical
      ispec_cpml = pml_physical_ispec(iface)
      do iside = 1, 2*NDIM
        if (pml_spec_physical(iside,ispec_cpml) .eq. 0) then
          pml_spec_physical(iside,ispec_cpml) = iface
          exit
        else
          if (iside .eq. 2*NDIM) print *, 'wrong setting pml_spec_physical'
        endif
      enddo
    enddo
  endif

  if (num_pml_physical > 0) then
    vtk_name = prname(1:len_trim(prname)) // 'pml_physical.vtk'
    open(unit=iout_vtk,file=trim(vtk_name),action='write',&
         form='formatted',iostat=ier)
    write(iout_vtk,'(a)') '# vtk DataFile Version 3.1'
    write(iout_vtk,'(a)') 'PML outer surface VTK file'
    write(iout_vtk,'(a)') 'ASCII'
    write(iout_vtk,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(iout_vtk,'(a,i15,a)') 'POINTS ', num_pml_physical*4, ' float'
    do iface = 1, num_pml_physical
      ispec_CPML = pml_physical_ispec(iface)
      if (ispec_CPML > NSPEC_CPML) &
        print *, 'ispec_CPML=',ispec_CPML,' NSPEC_CPML=',NSPEC_CPML,' myrank=',myrank
      i = pml_physical_ijk(1,1,iface)
      j = pml_physical_ijk(2,1,iface)
      k = pml_physical_ijk(3,1,iface)
      iglob_CPML = ibool_CPML(i,j,k,ispec_CPML)
      if (iglob_CPML > nglob_CPML) &
        print *, 'iglob_CPML=',iglob_CPML,' nglob_CPML=',nglob_CPML,&
                 ' myrank=',myrank,'\nispec_CPML=',ispec_CPML,' i=',i, ' j=',j, ' k=',k
      iglob = CPML_to_glob(iglob_CPML)
      if (iglob > NGLOB_AB) &
        print *, 'iglob=',iglob,' NGLOB_AB=',NGLOB_AB,' myrank=',myrank
      write(iout_vtk,'(3e18.6)') xstore(iglob), ystore(iglob), zstore(iglob)
      i = pml_physical_ijk(1,NGLLX,iface)
      j = pml_physical_ijk(2,NGLLX,iface)
      k = pml_physical_ijk(3,NGLLX,iface)
      iglob_CPML = ibool_CPML(i,j,k,ispec_CPML)
      if (iglob_CPML > nglob_CPML) &
        print *, 'iglob_CPML=',iglob_CPML,' nglob_CPML=',nglob_CPML,&
                 ' myrank=',myrank,'\nispec_CPML=',ispec_CPML,' i=',i, ' j=',j, ' k=',k
      iglob = CPML_to_glob(iglob_CPML)
      if (iglob > NGLOB_AB) &
        print *, 'iglob=',iglob,' NGLOB_AB=',NGLOB_AB,' myrank=',myrank
      write(iout_vtk,'(3e18.6)') xstore(iglob), ystore(iglob), zstore(iglob)
      i = pml_physical_ijk(1,NGLLX*NGLLY,iface)
      j = pml_physical_ijk(2,NGLLX*NGLLY,iface)
      k = pml_physical_ijk(3,NGLLX*NGLLY,iface)
      iglob_CPML = ibool_CPML(i,j,k,ispec_CPML)
      if (iglob_CPML > nglob_CPML) &
        print *, 'iglob_CPML=',iglob_CPML,' nglob_CPML=',nglob_CPML,&
                 ' myrank=',myrank,'\nispec_CPML=',ispec_CPML,' i=',i, ' j=',j, ' k=',k
      iglob = CPML_to_glob(iglob_CPML)
      if (iglob > NGLOB_AB) &
        print *, 'iglob=',iglob,' NGLOB_AB=',NGLOB_AB,' myrank=',myrank
      write(iout_vtk,'(3e18.6)') xstore(iglob), ystore(iglob), zstore(iglob)
      i = pml_physical_ijk(1,NGLLX*(NGLLY-1)+1,iface)
      j = pml_physical_ijk(2,NGLLX*(NGLLY-1)+1,iface)
      k = pml_physical_ijk(3,NGLLX*(NGLLY-1)+1,iface)
      iglob_CPML = ibool_CPML(i,j,k,ispec_CPML)
      if (iglob_CPML > nglob_CPML) &
        print *, 'iglob_CPML=',iglob_CPML,' nglob_CPML=',nglob_CPML,&
                 ' myrank=',myrank,'\nispec_CPML=',ispec_CPML,' i=',i, ' j=',j, ' k=',k
      iglob = CPML_to_glob(iglob_CPML)
      if (iglob > NGLOB_AB) &
        print *, 'iglob=',iglob,' NGLOB_AB=',NGLOB_AB,' myrank=',myrank
      write(iout_vtk,'(3e18.6)') xstore(iglob), ystore(iglob), zstore(iglob)
    enddo
    write(iout_vtk,'(a)') ''
    write(iout_vtk,'(a,2i12)') 'CELLS',num_pml_physical,num_pml_physical*5
    do iface = 1, num_pml_physical
      write(iout_vtk,'(5i12)') 4, (iface-1)*4, (iface-1)*4+1, &
                                  (iface-1)*4+2, (iface-1)*4+3
    enddo
    write(iout_vtk,'(a)') ''
    write(iout_vtk,'(a,i12)') 'CELL_TYPES', num_pml_physical
    do iface = 1, num_pml_physical
      write(iout_vtk,'(i12)') 9
    enddo
    close(iout_vtk)
  endif

  if (num_interfaces_PML > 0) then
    vtk_name = prname(1:len_trim(prname)) // 'pml_interface.vtk'
    open(unit=iout_vtk,file=trim(vtk_name),action='write',&
         form='formatted',iostat=ier)
    write(iout_vtk,'(a)') '# vtk DataFile Version 3.1'
    write(iout_vtk,'(a)') 'PML interface VTK file'
    write(iout_vtk,'(a)') 'ASCII'
    write(iout_vtk,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(iout_vtk,'(a,i15,a)') 'POINTS ', sum(nibool_interfaces_PML), ' float'
    do iface = 1, num_interfaces_PML
      do iside = 1, nibool_interfaces_PML(iface)
        iglob_CPML = ibool_interfaces_PML(iside,iface)
        iglob = CPML_to_glob(iglob_CPML)
        write(iout_vtk,'(3e18.6)') xstore(iglob), ystore(iglob), zstore(iglob)
      enddo
    enddo
    write(iout_vtk,'(a)') ''
    write(iout_vtk,'(a,2i12)') 'CELLS', sum(nibool_interfaces_PML), &
                                        sum(nibool_interfaces_PML) * 2
    do iside = 1, sum(nibool_interfaces_PML)
      write(iout_vtk,'(2i12)') 1, iside - 1
    enddo
    write(iout_vtk,'(a)') ''
    write(iout_vtk,'(a,i12)') 'CELL_TYPES', sum(nibool_interfaces_PML)
    do iside = 1, sum(nibool_interfaces_PML)
      write(iout_vtk,'(i12)') 1
    enddo
    write(iout_vtk,'(a)') ''
    write(iout_vtk,'(a,i12)') 'POINT_DATA', sum(nibool_interfaces_PML)
    write(iout_vtk,'(a)') 'SCALARS i_neighbour int'
    write(iout_vtk,'(a)') 'LOOKUP_TABLE default'
    do iface = 1, num_interfaces_PML
      do iside = 1, nibool_interfaces_PML(iface)
        write(iout_vtk,'(i12)') iface
      enddo
    enddo
    close(iout_vtk)
  endif

  if (nglob_pml_in > 0) then
    vtk_name = prname(1:len_trim(prname)) // 'pml_internal.vtk'
    open(unit=iout_vtk,file=trim(vtk_name),action='write',&
         form='formatted',iostat=ier)
    write(iout_vtk,'(a)') '# vtk DataFile Version 3.1'
    write(iout_vtk,'(a)') 'PML interface VTK file'
    write(iout_vtk,'(a)') 'ASCII'
    write(iout_vtk,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(iout_vtk,'(a,i15,a)') 'POINTS ', nglob_pml_in, ' float'
    do iside = 1, nglob_pml_in
      iglob_CPML = pml_in_iglob(iside)
      iglob = CPML_to_glob(iglob_CPML)
      write(iout_vtk,'(3e18.6)') xstore(iglob), ystore(iglob), zstore(iglob)
    enddo
    write(iout_vtk,'(a)') ''
    write(iout_vtk,'(a,2i12)') 'CELLS', nglob_pml_in, nglob_pml_in * 2
    do iside = 1, nglob_pml_in
      write(iout_vtk,'(2i12)') 1, iside - 1
    enddo
    write(iout_vtk,'(a)') ''
    write(iout_vtk,'(a,i12)') 'CELL_TYPES',nglob_pml_in
    do iside = 1, nglob_pml_in
      write(iout_vtk,'(i12)') 1
    enddo
    close(iout_vtk)
  endif

  allocate(buffer_send_matrix_PML(NDIM,NDIM,max_nibool_interfaces_PML,&
            num_interfaces_PML),stat=ier)
  allocate(buffer_recv_matrix_PML(NDIM,NDIM,max_nibool_interfaces_PML,&
            num_interfaces_PML),stat=ier)
  allocate(request_send_matrix_PML(num_interfaces_PML),stat=ier)
  allocate(request_recv_matrix_PML(num_interfaces_PML),stat=ier)
  !if (NGLOB_AB > 0) then
  !  allocate(rvolume_ext(NGLOB_AB))
  !  allocate(Qt_t_ext(NDIM,NDIM,NGLOB_AB))
  !  rvolume_ext(:) = 0.0
  !  Qt_t_ext(:,:,:) = 0.0
  !endif

end subroutine read_adepml_files
