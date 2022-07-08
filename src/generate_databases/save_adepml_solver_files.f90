subroutine save_adepml_solver_files()
  use create_regions_mesh_ext_par, only: MAX_STRING_LEN, prname
  use generate_databases_par, only: IOUT
  use generate_databases_par, only: nspec_CPML, CPML_to_spec, CPML_regions,&
         is_CPML
  use generate_databases_par, only: &
    pml_d,pml_beta,pml_kappa,num_pml_out,pml_out_ispec,pml_out_ijk
    !nglob_CPML,CPML_to_glob,ibool_CPML
    !,num_interfaces_PML,&
    ! max_nibool_interfaces_PML,nibool_interfaces_PML,ibool_interfaces_PML,&
    ! my_neighbors_PML
  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename
  filename = prname(1:len_trim(prname))//'adepml_damping_indexing.bin'
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening database proc######_adepml_damping_indexing.bin'
  write(IOUT) nspec_cpml
  if (nspec_cpml > 0) then
    write(IOUT) CPML_regions
    write(IOUT) CPML_to_spec
    write(IOUT) is_CPML
    write(IOUT) pml_d
    write(IOUT) pml_beta
    write(IOUT) pml_kappa
    !write(IOUT) nglob_CPML
    !write(IOUT) ibool_CPML
    !write(IOUT) CPML_to_glob
  endif
  write(IOUT) num_pml_out
  if (num_pml_out > 0) then
    write(IOUT) pml_out_ispec
    write(IOUT) pml_out_ijk
  endif
  !write(IOUT) num_interfaces_PML
  !if (num_interfaces_PML > 0) then
  !  write(IOUT) max_nibool_interfaces_PML
  !  write(IOUT) nibool_interfaces_PML
  !  write(IOUT) ibool_interfaces_PML
  !  write(IOUT) my_neighbors_PML
  !endif
  deallocate(pml_d)
  deallocate(pml_kappa)
  deallocate(pml_beta)
  !deallocate(ibool_CPML,CPML_to_glob)
  if (num_pml_out > 0) then
    deallocate(pml_out_ispec)
    deallocate(pml_out_ijk)
  endif
  !if (num_interfaces_PML > 0) then
  !  deallocate(nibool_interfaces_PML)
  !  deallocate(ibool_interfaces_PML)
  !  deallocate(my_neighbors_PML)
  !endif
  close(IOUT)
end subroutine save_adepml_solver_files

subroutine save_adepml_solver_files_glob()
  use create_regions_mesh_ext_par, only: MAX_STRING_LEN, prname
  use generate_databases_par, only: IOUT
  use generate_databases_par, only: &
    nglob_CPML,CPML_to_glob,ibool_CPML,num_interfaces_PML,&
    max_nibool_interfaces_PML,nibool_interfaces_PML,ibool_interfaces_PML,&
    my_neighbors_PML,nglob_pml_in,pml_in_iglob,nglob_dirichlet,iglob_dirichlet
  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename
  filename = prname(1:len_trim(prname))//'adepml_glob_indexing.bin'
  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening database proc######_adepml_glob_indexing.bin'
  write(IOUT) nglob_CPML
  if (nglob_CPML > 0) then
    write(IOUT) ibool_CPML
    write(IOUT) CPML_to_glob
  endif
  write(IOUT) num_interfaces_PML
  if (num_interfaces_PML > 0) then
    write(IOUT) max_nibool_interfaces_PML
    write(IOUT) nibool_interfaces_PML
    write(IOUT) ibool_interfaces_PML
    write(IOUT) my_neighbors_PML
  endif
  write(IOUT) nglob_pml_in
  if (nglob_pml_in > 0) then
    write(IOUT) pml_in_iglob
    deallocate(pml_in_iglob)
  endif
  write(IOUT) nglob_dirichlet
  if (nglob_dirichlet > 0) then
    write(IOUT) iglob_dirichlet
    deallocate(iglob_dirichlet)
  endif
  if (nglob_CPML > 0) deallocate(ibool_CPML,CPML_to_glob)
  if (num_interfaces_PML > 0) then
    deallocate(nibool_interfaces_PML)
    deallocate(ibool_interfaces_PML)
    deallocate(my_neighbors_PML)
  endif
  
  close(IOUT)
end subroutine save_adepml_solver_files_glob
