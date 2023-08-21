!! Solving the wavefield discontinuity problem with a non-split-node
!! scheme
!! Tianshi Liu, 2023.5

subroutine read_mesh_databases_wavefield_discontinuity()
  use specfem_par, only: IIN, CUSTOM_REAL, &
                         NSPEC_AB, NGLLX, NGLLY, NGLLZ, NDIM, NGLLSQUARE
  use wavefield_discontinuity_solver
  implicit none
  allocate(ispec_to_elem_wd(NSPEC_AB))
  read(IIN) ispec_to_elem_wd
  read(IIN) nglob_wd
  read(IIN) nspec_wd
  allocate(ibool_wd(NGLLX, NGLLY, NGLLZ, nspec_wd), &
           boundary_to_iglob_wd(nglob_wd), &
           mass_in_wd(nglob_wd))
  read(IIN) ibool_wd
  read(IIN) boundary_to_iglob_wd
  read(IIN) mass_in_wd
  read(IIN) nfaces_wd
  allocate(face_ijk_wd(NDIM, NGLLSQUARE, nfaces_wd), &
           face_ispec_wd(nfaces_wd), &
           face_normal_wd(NDIM, NGLLSQUARE, nfaces_wd), &
           face_jacobian2Dw_wd(NGLLSQUARE, nfaces_wd))
  read(IIN) face_ijk_wd
  read(IIN) face_ispec_wd
  read(IIN) face_normal_wd
  read(IIN) face_jacobian2Dw_wd
  allocate(displ_wd(NDIM, nglob_wd), &
           accel_wd(NDIM, nglob_wd), &
           traction_wd(NDIM, NGLLSQUARE, nfaces_wd))
end subroutine read_mesh_databases_wavefield_discontinuity

subroutine open_wavefield_discontinuity_file()
  use specfem_par, only: prname
  use wavefield_discontinuity_solver, only: IFILE_WAVEFIELD_DISCONTINUITY
  implicit none
  open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
       file=trim(prname)//'wavefield_discontinuity.bin', &
       status='old',action='read',form='unformatted')
end subroutine open_wavefield_discontinuity_file

subroutine read_wavefield_discontinuity_file()
  use wavefield_discontinuity_solver, only: IFILE_WAVEFIELD_DISCONTINUITY, &
                 displ_wd, accel_wd, traction_wd
  implicit none
  read(IFILE_WAVEFIELD_DISCONTINUITY) displ_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) accel_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) traction_wd
end subroutine read_wavefield_discontinuity_file

subroutine finalize_wavefield_discontinuity()
  use wavefield_discontinuity_solver
  implicit none
  close(IFILE_WAVEFIELD_DISCONTINUITY)
  deallocate(ispec_to_elem_wd, ibool_wd, boundary_to_iglob_wd, mass_in_wd, &
             face_ijk_wd, face_ispec_wd, face_normal_wd, face_jacobian2Dw_wd, &
             displ_wd, accel_wd, traction_wd)
end subroutine finalize_wavefield_discontinuity

subroutine add_displacement_discontinuity_element(ispec, dummyx_loc, &
                                                  dummyy_loc, dummyz_loc)
  use specfem_par, only: CUSTOM_REAL,NGLLX, NGLLY, NGLLZ
  use wavefield_discontinuity_solver, only: ispec_to_elem_wd, ibool_wd, &
                                            displ_wd
  implicit none
  integer, intent(in) :: ispec
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: dummyx_loc, &
                        dummyy_loc, dummyz_loc
  integer :: ispec_wd, i, j, k, iglob_wd
  ispec_wd = ispec_to_elem_wd(ispec)
  if (ispec_wd /= 0) then
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob_wd = ibool_wd(i,j,k,ispec_wd)
          if (iglob_wd /= 0) then
            dummyx_loc(i,j,k) = dummyx_loc(i,j,k) + displ_wd(1, iglob_wd)
            dummyy_loc(i,j,k) = dummyy_loc(i,j,k) + displ_wd(2, iglob_wd)
            dummyz_loc(i,j,k) = dummyz_loc(i,j,k) + displ_wd(3, iglob_wd)
          endif
        enddo
      enddo
    enddo
  endif 
end subroutine add_displacement_discontinuity_element

subroutine add_traction_discontinuity()
  use wavefield_discontinuity_solver, only: accel_wd, mass_in_wd, nglob_wd, &
                                    boundary_to_iglob_wd, traction_wd,&
                                    nfaces_wd, face_ijk_wd, face_ispec_wd, &
                                    face_jacobian2Dw_wd
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NGLLSQUARE, ibool
  use specfem_par_elastic, only: accel
  implicit none
  integer :: iglob_wd, iglob, ispec, i, j, k, iface_wd, igll
  real(kind=CUSTOM_REAL) :: jacobianw
  do iglob_wd = 1, nglob_wd
    iglob = boundary_to_iglob_wd(iglob_wd)
    accel(:,iglob) = accel(:,iglob) - &
                     accel_wd(:,iglob_wd) * mass_in_wd(iglob_wd)
  enddo
  do iface_wd = 1, nfaces_wd
    do igll = 1, NGLLSQUARE
      i = face_ijk_wd(1, igll, iface_wd)
      j = face_ijk_wd(2, igll, iface_wd)
      k = face_ijk_wd(3, igll, iface_wd)
      ispec = face_ispec_wd(iface_wd)
      iglob = ibool(i,j,k,ispec)
      jacobianw = face_jacobian2Dw_wd(igll, iface_wd)
      accel(:,iglob) = accel(:,iglob) + &
                     traction_wd(:,igll,iface_wd) * jacobianw
    enddo
  enddo
end subroutine add_traction_discontinuity
