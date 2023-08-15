!! Setting up the surface on which wavefield discontinuity condition is
!! imposed
!! Tianshi Liu, 2023.5


subroutine read_partition_files_wavefield_discontinuity()
  use generate_databases_par, only: IIN
  use wavefield_discontinuity_generate_databases, only: &
                              nb_wd, boundary_to_ispec_wd, side_wd
  implicit none
  !integer :: ier
  read(IIN) nb_wd
  allocate(boundary_to_ispec_wd(nb_wd), side_wd(nb_wd))
  read(IIN) boundary_to_ispec_wd
  read(IIN) side_wd
end subroutine read_partition_files_wavefield_discontinuity

subroutine save_arrays_solver_mesh_wavefield_discontinuity()
  use generate_databases_par, only: IOUT
  use wavefield_discontinuity_generate_databases
  implicit none
  write(IOUT) ispec_to_elem_wd
  write(IOUT) nglob_wd
  write(IOUT) nspec_wd
  write(IOUT) ibool_wd
  write(IOUT) boundary_to_iglob_wd
  write(IOUT) mass_in_wd
  write(IOUT) nfaces_wd
  write(IOUT) face_ijk_wd
  write(IOUT) face_ispec_wd
  write(IOUT) face_normal_wd
  write(IOUT) face_jacobian2dw_wd
  deallocate(boundary_to_ispec_wd, side_wd)
  deallocate(ispec_to_elem_wd, ibool_wd, boundary_to_iglob_wd, mass_in_wd, &
             face_ijk_wd, face_ispec_wd, face_normal_wd, face_jacobian2dw_wd)
end subroutine save_arrays_solver_mesh_wavefield_discontinuity

subroutine setup_boundary_wavefield_discontinuity()
  use generate_databases_par, only: NDIM, NGLLX, NGLLY, NGLLZ, CUSTOM_REAL, &
                                    NGLLSQUARE, ibool, NSPEC_AB
                                    
  use generate_databases_par, only: myrank ! debug
  use wavefield_discontinuity_generate_databases
  implicit none
  integer :: i1, i2, j1, j2, k1, k2, i, j, k
  integer :: ib, iside, ispec, iglob, ispec_wd, iglob_wd
  logical :: is_face
  integer, dimension(:), allocatable :: elem_list_temp, gllp_list_temp
  real(kind=CUSTOM_REAL) :: val_mass
  allocate(ispec_to_elem_wd(NSPEC_AB))
  allocate(elem_list_temp(nb_wd), gllp_list_temp(nb_wd*NGLLSQUARE))
  nspec_wd = 0
  nglob_wd = 0
  nfaces_wd = 0
  ispec_to_elem_wd(:) = 0
  elem_list_temp(:) = 0
  gllp_list_temp(:) = 0  
  do ib = 1, nb_wd
    ispec = boundary_to_ispec_wd(ib)
    iside = side_wd(ib)
    call get_points_boundary_wd(iside, i1, i2, j1, j2, k1, k2, is_face)
    call find_point_wd(ispec, elem_list_temp, nspec_wd, ispec_wd)
    if (ispec_wd == 0) then
      ispec_to_elem_wd(ispec) = nspec_wd + 1
      elem_list_temp(nspec_wd+1) = ispec
      nspec_wd = nspec_wd + 1
    endif
    do i=i1,i2; do j=j1,j2; do k=k1,k2
      iglob = ibool(i,j,k,ispec)
      call find_point_wd(iglob, gllp_list_temp, nglob_wd, iglob_wd)
      if (iglob_wd == 0) then
        gllp_list_temp(nglob_wd+1) = iglob
        nglob_wd = nglob_wd + 1
      endif
    enddo; enddo; enddo
    if (is_face) nfaces_wd = nfaces_wd + 1
  enddo
  if (myrank == 0) then
    ! debug
    print *, nspec_wd, ' elements, ', nfaces_wd, ' faces'
  endif
  allocate(ibool_wd(NGLLX, NGLLY, NGLLZ, nspec_wd), &
           boundary_to_iglob_wd(nglob_wd), &
           mass_in_wd(nglob_wd), &
           face_ijk_wd(NDIM, NGLLSQUARE, nfaces_wd), &
           face_ispec_wd(nfaces_wd), &
           face_normal_wd(NDIM, NGLLSQUARE, nfaces_wd), &
           face_jacobian2Dw_wd(NGLLSQUARE, nfaces_wd))
  ibool_wd(:,:,:,:) = 0
  mass_in_wd(:) = 0.0
  nfaces_wd = 0
  do ib = 1, nb_wd
    ispec = boundary_to_ispec_wd(ib)
    iside = side_wd(ib)
    call get_points_boundary_wd(iside, i1, i2, j1, j2, k1, k2, is_face)
    call find_point_wd(ispec, elem_list_temp, nspec_wd, ispec_wd)
    do i=i1,i2; do j=j1,j2; do k=k1,k2
      iglob = ibool(i,j,k,ispec)
      call find_point_wd(iglob, gllp_list_temp, nglob_wd, iglob_wd)
      ibool_wd(i,j,k,ispec_wd) = iglob_wd
      boundary_to_iglob_wd(iglob_wd) = iglob
      !call get_mass_wd(i, j, k, ispec, val_mass)
      !mass_in_wd(iglob_wd) = mass_in_wd(iglob_wd) + val_mass
    enddo; enddo; enddo
    if (is_face) then
      nfaces_wd = nfaces_wd + 1
      face_ispec_wd(nfaces_wd) = ispec
      call get_face_wd(ispec, iside, face_ijk_wd(:,:,nfaces_wd), &
                       face_normal_wd(:,:,nfaces_wd), &
                       face_jacobian2Dw_wd(:,nfaces_wd))
    endif
  enddo
  do ispec = 1, NSPEC_AB
    ispec_wd = ispec_to_elem_wd(ispec)
    if (ispec_wd > 0) then
      do k=1,NGLLZ; do j=1,NGLLY; do i=1, NGLLX
        iglob_wd= ibool_wd(i,j,k,ispec_wd)
        if (iglob_wd > 0) then
          call get_mass_wd(i, j, k, ispec, val_mass)
          mass_in_wd(iglob_wd) = mass_in_wd(iglob_wd) + val_mass
        endif
      enddo; enddo; enddo
    endif
  enddo
  call write_discontinuity_surface_file()
end subroutine setup_boundary_wavefield_discontinuity

  subroutine write_discontinuity_surface_file()
  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,NGLLSQUARE,CUSTOM_REAL
  use generate_databases_par, only: ibool, prname, myrank, LOCAL_PATH
  use create_regions_mesh_ext_par, only: xstore_dummy, ystore_dummy, &
                                         zstore_dummy
  use wavefield_discontinuity_generate_databases, only: nglob_wd, nfaces_wd, &
          boundary_to_iglob_wd, face_ijk_wd, face_normal_wd, face_ispec_wd, &
          face_jacobian2Dw_wd, mass_in_wd, &
          IFILE_WAVEFIELD_DISCONTINUITY
  implicit none
  integer :: i,j,k,ispec,iglob,igll,iglob_wd,iface_wd
  real(CUSTOM_REAL) :: tx, ty, tz
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
       file=prname(1:len_trim(prname))//'wavefield_discontinuity_points',&
       action='write', form='formatted')
  do iglob_wd = 1, nglob_wd
    iglob = boundary_to_iglob_wd(iglob_wd)
    write(IFILE_WAVEFIELD_DISCONTINUITY, '(4e20.5)') xstore_dummy(iglob), &
       ystore_dummy(iglob), zstore_dummy(iglob), mass_in_wd(iglob_wd)
  enddo
  close(IFILE_WAVEFIELD_DISCONTINUITY)
  open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
       file=prname(1:len_trim(prname))//'wavefield_discontinuity_faces',&
       action='write', form='formatted')
  do iface_wd = 1, nfaces_wd
    ispec = face_ispec_wd(iface_wd)
    do igll = 1, NGLLSQUARE
      i = face_ijk_wd(1,igll,iface_wd)
      j = face_ijk_wd(2,igll,iface_wd)
      k = face_ijk_wd(3,igll,iface_wd)
      iglob = ibool(i,j,k,ispec)
      tx = face_normal_wd(1,igll,iface_wd)
      ty = face_normal_wd(2,igll,iface_wd)
      tz = face_normal_wd(3,igll,iface_wd)
      write(IFILE_WAVEFIELD_DISCONTINUITY, '(7e20.5)') xstore_dummy(iglob), &
       ystore_dummy(iglob), zstore_dummy(iglob), tx, ty, tz, &
       face_jacobian2Dw_wd(igll, iface_wd)
    enddo
  enddo
  end subroutine write_discontinuity_surface_file

  subroutine get_face_wd(ispec, iside, face_ijk_b, face_normal_b, &
                         face_jacobian2Dw_b)
  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,NGNOD2D_FOUR_CORNERS, &
                       CUSTOM_REAL, NGLLSQUARE
  use generate_databases_par, only: NGNOD2D, ibool, NSPEC_AB, NGLOB_AB, &
                                    xstore, ystore, zstore, myrank
  use create_regions_mesh_ext_par, only: xstore_dummy, ystore_dummy, &
                                         zstore_dummy, nglob_dummy, &
            dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
            wgllwgll_xy,wgllwgll_xz,wgllwgll_yz
  implicit none
  integer, intent(in) :: ispec, iside
  integer :: iface, igll, i, j
  integer :: ijk_face(NDIM, NGLLX, NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM, NGLLX, NGLLY), &
                                         jacobian2Dw_face(NGLLX, NGLLY)
  real(kind=CUSTOM_REAL), dimension(NDIM) :: lnormal
  integer, intent(out) :: face_ijk_b(NDIM, NGLLSQUARE)
  real(kind=CUSTOM_REAL), intent(out) :: face_normal_b(NDIM, NGLLSQUARE), &
                                         face_jacobian2Dw_b(NGLLSQUARE)
  real(kind=CUSTOM_REAL), dimension(NGNOD2D_FOUR_CORNERS) :: &
                                         xcoord,ycoord,zcoord
  if (nglob_dummy /= NGLOB_AB) &
    call exit_MPI_without_rank('nglob_dummy not equal to NGLOB_AB')
  select case (iside)
    case (21)
      ! bottom
      xcoord(1) = xstore(1,1,1,ispec)
      xcoord(2) = xstore(1,NGLLY,1,ispec)
      xcoord(3) = xstore(NGLLX,NGLLY,1,ispec)
      xcoord(4) = xstore(NGLLX,1,1,ispec)
      ycoord(1) = ystore(1,1,1,ispec)
      ycoord(2) = ystore(1,NGLLY,1,ispec)
      ycoord(3) = ystore(NGLLX,NGLLY,1,ispec)
      ycoord(4) = ystore(NGLLX,1,1,ispec)
      zcoord(1) = zstore(1,1,1,ispec)
      zcoord(2) = zstore(1,NGLLY,1,ispec)
      zcoord(3) = zstore(NGLLX,NGLLY,1,ispec)
      zcoord(4) = zstore(NGLLX,1,1,ispec)
      ! sets face id of reference element associated with this face
      call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,NSPEC_AB,NGLOB_AB, &
                             xstore_dummy,ystore_dummy,zstore_dummy,iface)
      if (iface /= 5) &
        call exit_MPI_without_rank('incorrect iface')
      ! ijk indices of GLL points on face
      call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)
      ! weighted jacobian and normal
      call get_jacobian_boundary_face(myrank, NSPEC_AB, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)
      ! normal convention: points away from element
      ! switch normal direction if necessary
      do j = 1,NGLLY
        do i = 1,NGLLX
          lnormal(:) = normal_face(:,i,j)
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                  ibool,NSPEC_AB,nglob_dummy, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  lnormal)
          normal_face(:,i,j) = lnormal(:)
        enddo
      enddo
      ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
      igll = 0
      do j = 1,NGLLY
         do i = 1,NGLLX
           igll = igll+1
           face_ijk_b(:,igll) = ijk_face(:,i,j)
           face_jacobian2Dw_b(igll) = jacobian2Dw_face(i,j)
           face_normal_b(:,igll) = normal_face(:,i,j)
         enddo
      enddo
    case (25)
      ! left, xmin
      xcoord(1) = xstore(1,1,1,ispec)
      xcoord(2) = xstore(1,1,NGLLZ,ispec)
      xcoord(3) = xstore(1,NGLLY,NGLLZ,ispec)
      xcoord(4) = xstore(1,NGLLY,1,ispec)
      ycoord(1) = ystore(1,1,1,ispec)
      ycoord(2) = ystore(1,1,NGLLZ,ispec)
      ycoord(3) = ystore(1,NGLLY,NGLLZ,ispec)
      ycoord(4) = ystore(1,NGLLY,1,ispec)
      zcoord(1) = zstore(1,1,1,ispec)
      zcoord(2) = zstore(1,1,NGLLZ,ispec)
      zcoord(3) = zstore(1,NGLLY,NGLLZ,ispec)
      zcoord(4) = zstore(1,NGLLY,1,ispec)
      ! sets face id of reference element associated with this face
      call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,NSPEC_AB,NGLOB_AB, &
                             xstore_dummy,ystore_dummy,zstore_dummy,iface)
      if (iface /= 1) &
        call exit_MPI_without_rank('incorrect iface')
      ! ijk indices of GLL points on face
      call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)
      ! weighted jacobian and normal
      call get_jacobian_boundary_face(myrank, NSPEC_AB, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)
      ! normal convention: points away from element
      ! switch normal direction if necessary
      do j = 1,NGLLZ
        do i = 1,NGLLX
          lnormal(:) = normal_face(:,i,j)
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                  ibool,NSPEC_AB,nglob_dummy, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  lnormal)
          normal_face(:,i,j) = lnormal(:)
        enddo
      enddo
      ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
      igll = 0
      do j = 1,NGLLZ
         do i = 1,NGLLX
           igll = igll+1
           face_ijk_b(:,igll) = ijk_face(:,i,j)
           face_jacobian2Dw_b(igll) = jacobian2Dw_face(i,j)
           face_normal_b(:,igll) = normal_face(:,i,j)
         enddo
      enddo
    case (22)
      ! front, ymin
      xcoord(1) = xstore(1,1,1,ispec)
      xcoord(2) = xstore(NGLLX,1,1,ispec)
      xcoord(3) = xstore(NGLLX,1,NGLLZ,ispec)
      xcoord(4) = xstore(1,1,NGLLZ,ispec)
      ycoord(1) = ystore(1,1,1,ispec)
      ycoord(2) = ystore(NGLLX,1,1,ispec)
      ycoord(3) = ystore(NGLLX,1,NGLLZ,ispec)
      ycoord(4) = ystore(1,1,NGLLZ,ispec)
      zcoord(1) = zstore(1,1,1,ispec)
      zcoord(2) = zstore(NGLLX,1,1,ispec)
      zcoord(3) = zstore(NGLLX,1,NGLLZ,ispec)
      zcoord(4) = zstore(1,1,NGLLZ,ispec)
      ! sets face id of reference element associated with this face
      call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,NSPEC_AB,NGLOB_AB, &
                             xstore_dummy,ystore_dummy,zstore_dummy,iface)
      if (iface /= 3) &
        call exit_MPI_without_rank('incorrect iface')
      ! ijk indices of GLL points on face
      call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)
      ! weighted jacobian and normal
      call get_jacobian_boundary_face(myrank, NSPEC_AB, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)
      ! normal convention: points away from element
      ! switch normal direction if necessary
      do j = 1,NGLLZ
        do i = 1,NGLLY
          lnormal(:) = normal_face(:,i,j)
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                  ibool,NSPEC_AB,nglob_dummy, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  lnormal)
          normal_face(:,i,j) = lnormal(:)
        enddo
      enddo
      ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
      igll = 0
      do j = 1,NGLLZ
         do i = 1,NGLLY
           igll = igll+1
           face_ijk_b(:,igll) = ijk_face(:,i,j)
           face_jacobian2Dw_b(igll) = jacobian2Dw_face(i,j)
           face_normal_b(:,igll) = normal_face(:,i,j)
         enddo
      enddo
    case (23)
      ! right, xmax
      xcoord(1) = xstore(NGLLX,1,1,ispec)
      xcoord(2) = xstore(NGLLX,NGLLY,1,ispec)
      xcoord(3) = xstore(NGLLX,NGLLY,NGLLZ,ispec)
      xcoord(4) = xstore(NGLLX,1,NGLLZ,ispec)
      ycoord(1) = ystore(NGLLX,1,1,ispec)
      ycoord(2) = ystore(NGLLX,NGLLY,1,ispec)
      ycoord(3) = ystore(NGLLX,NGLLY,NGLLZ,ispec)
      ycoord(4) = ystore(NGLLX,1,NGLLZ,ispec)
      zcoord(1) = zstore(NGLLX,1,1,ispec)
      zcoord(2) = zstore(NGLLX,NGLLY,1,ispec)
      zcoord(3) = zstore(NGLLX,NGLLY,NGLLZ,ispec)
      zcoord(4) = zstore(NGLLX,1,NGLLZ,ispec)
      ! sets face id of reference element associated with this face
      call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,NSPEC_AB,NGLOB_AB, &
                             xstore_dummy,ystore_dummy,zstore_dummy,iface)
      if (iface /= 2) &
        call exit_MPI_without_rank('incorrect iface')
      ! ijk indices of GLL points on face
      call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)
      ! weighted jacobian and normal
      call get_jacobian_boundary_face(myrank, NSPEC_AB, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)
      ! normal convention: points away from element
      ! switch normal direction if necessary
      do j = 1,NGLLZ
        do i = 1,NGLLX
          lnormal(:) = normal_face(:,i,j)
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                  ibool,NSPEC_AB,nglob_dummy, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  lnormal)
          normal_face(:,i,j) = lnormal(:)
        enddo
      enddo
      ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
      igll = 0
      do j = 1,NGLLZ
         do i = 1,NGLLX
           igll = igll+1
           face_ijk_b(:,igll) = ijk_face(:,i,j)
           face_jacobian2Dw_b(igll) = jacobian2Dw_face(i,j)
           face_normal_b(:,igll) = normal_face(:,i,j)
         enddo
      enddo
    case (24)
      ! back, ymax
      xcoord(1) = xstore(1,NGLLY,1,ispec)
      xcoord(2) = xstore(1,NGLLY,NGLLZ,ispec)
      xcoord(3) = xstore(NGLLX,NGLLY,NGLLZ,ispec)
      xcoord(4) = xstore(NGLLX,NGLLY,1,ispec)
      ycoord(1) = ystore(1,NGLLY,1,ispec)
      ycoord(2) = ystore(1,NGLLY,NGLLZ,ispec)
      ycoord(3) = ystore(NGLLX,NGLLY,NGLLZ,ispec)
      ycoord(4) = ystore(NGLLX,NGLLY,1,ispec)
      zcoord(1) = zstore(1,NGLLY,1,ispec)
      zcoord(2) = zstore(1,NGLLY,NGLLZ,ispec)
      zcoord(3) = zstore(NGLLX,NGLLY,NGLLZ,ispec)
      zcoord(4) = zstore(NGLLX,NGLLY,1,ispec)
      ! sets face id of reference element associated with this face
      call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,NSPEC_AB,NGLOB_AB, &
                             xstore_dummy,ystore_dummy,zstore_dummy,iface)
      if (iface /= 4) &
        call exit_MPI_without_rank('incorrect iface')
      ! ijk indices of GLL points on face
      call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)
      ! weighted jacobian and normal
      call get_jacobian_boundary_face(myrank, NSPEC_AB, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ,NGNOD2D)
      ! normal convention: points away from element
      ! switch normal direction if necessary
      do j = 1,NGLLZ
        do i = 1,NGLLY
          lnormal(:) = normal_face(:,i,j)
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                  ibool,NSPEC_AB,nglob_dummy, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  lnormal)
          normal_face(:,i,j) = lnormal(:)
        enddo
      enddo
      ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
      igll = 0
      do j = 1,NGLLZ
         do i = 1,NGLLY
           igll = igll+1
           face_ijk_b(:,igll) = ijk_face(:,i,j)
           face_jacobian2Dw_b(igll) = jacobian2Dw_face(i,j)
           face_normal_b(:,igll) = normal_face(:,i,j)
         enddo
      enddo
    case (26)
      ! top, zmax
      xcoord(1) = xstore(1,1,NGLLZ,ispec)
      xcoord(2) = xstore(NGLLX,1,NGLLZ,ispec)
      xcoord(3) = xstore(NGLLX,NGLLY,NGLLZ,ispec)
      xcoord(4) = xstore(1,NGLLY,NGLLZ,ispec)
      ycoord(1) = ystore(1,1,NGLLZ,ispec)
      ycoord(2) = ystore(NGLLX,1,NGLLZ,ispec)
      ycoord(3) = ystore(NGLLX,NGLLY,NGLLZ,ispec)
      ycoord(4) = ystore(1,NGLLY,NGLLZ,ispec)
      zcoord(1) = zstore(1,1,NGLLZ,ispec)
      zcoord(2) = zstore(NGLLX,1,NGLLZ,ispec)
      zcoord(3) = zstore(NGLLX,NGLLY,NGLLZ,ispec)
      zcoord(4) = zstore(1,NGLLY,NGLLZ,ispec)
      ! sets face id of reference element associated with this face
      call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,NSPEC_AB,NGLOB_AB, &
                             xstore_dummy,ystore_dummy,zstore_dummy,iface)
      if (iface /= 6) &
        call exit_MPI_without_rank('incorrect iface')
      ! ijk indices of GLL points on face
      call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)
      ! weighted jacobian and normal
      call get_jacobian_boundary_face(myrank, NSPEC_AB, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY,NGNOD2D)
      ! normal convention: points away from element
      ! switch normal direction if necessary
      do j = 1,NGLLY
        do i = 1,NGLLX
          lnormal(:) = normal_face(:,i,j)
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                  ibool,NSPEC_AB,nglob_dummy, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  lnormal)
          normal_face(:,i,j) = lnormal(:)
        enddo
      enddo
      ! GLL points -- assuming NGLLX = NGLLY = NGLLZ
      igll = 0
      do j = 1,NGLLY
         do i = 1,NGLLX
           igll = igll+1
           face_ijk_b(:,igll) = ijk_face(:,i,j)
           face_jacobian2Dw_b(igll) = jacobian2Dw_face(i,j)
           face_normal_b(:,igll) = normal_face(:,i,j)
         enddo
      enddo
  end select
  end subroutine get_face_wd

  subroutine get_mass_wd(i, j, k, ispec, val_mass)
  use generate_databases_par, only: CUSTOM_REAL
  use create_regions_mesh_ext_par, only: irregular_element_number, &
                                wxgll, wygll, wzgll, jacobianstore, &
                                rhostore, jacobian_regular
  implicit none
  integer :: i, j, k, ispec, ispec_irreg
  real(kind=CUSTOM_REAL) :: val_mass, jacobianl
  double precision :: weight
  ispec_irreg = irregular_element_number(ispec)
  if (ispec_irreg == 0) then
    jacobianl = jacobian_regular
  else
    jacobianl = jacobianstore(i,j,k,ispec_irreg)
  endif
  weight = wxgll(i)*wygll(j)*wzgll(k)
  val_mass = real(dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)), &
                  kind=CUSTOM_REAL)
  end subroutine get_mass_wd

  subroutine find_point_wd(p, plist, len, ip)
  implicit none
  integer, intent(in) :: p, len
  integer, dimension(len), intent(in) :: plist
  integer, intent(out) :: ip
  integer :: i
  ip = 0
  do i = 1, len
    if (plist(i) == p) then
      ip = i
      exit
    endif
  enddo
  end subroutine find_point_wd

  subroutine get_points_boundary_wd(iside, i1, i2, j1, j2, k1, k2, is_face)
  use generate_databases_par, only: NGLLX, NGLLY, NGLLZ
  implicit none
  integer, intent(in) :: iside
  integer, intent(out) :: i1, i2, j1, j2, k1, k2
  logical, intent(out) :: is_face
  integer, parameter :: NGNOD = 27, THREE = 3
  integer, dimension(NGNOD) :: iaddx, iaddy, iaddz
  integer :: i, j, k
  call usual_hex_nodes(NGNOD, iaddx, iaddy, iaddz)
  iaddx = iaddx * 2 + 1
  iaddy = iaddy * 2 + 1
  iaddz = iaddz * 2 + 1
  i = iaddx(iside)
  j = iaddy(iside)
  k = iaddz(iside)
  is_face = .false.
  if ((iside >= 1) .and. (iside <= 8)) then
    ! a vertex is on the boundary
    i1 = i; i2 = i
    j1 = j; j2 = j
    k1 = k; k2 = k
  else if ((iside >= 9) .and. (iside <= 20)) then
    ! an edge is on the boundary
    if (i == THREE) then
      i1 = 1; i2 = NGLLX
      j1 = j; j2 = j
      k1 = k; k2 = k
    else if (j == THREE) then
      i1 = i; i2 = i
      j1 = 1; j2 = NGLLY
      k1 = k; k2 = k
    else if (k == THREE) then
      i1 = i; i2 = i
      j1 = j; j2 = j
      k1 = 1; k2 = NGLLZ
    else
      call exit_MPI_without_rank('incorrect wavefield discontinuity point')
    endif
  else if ((iside >= 21) .and. (iside <= 26)) then
    ! a face is on the boundary
    is_face = .true.
    if (i /= THREE) then
      i1 = i; i2 = i
      j1 = 1; j2 = NGLLY
      k1 = 1; k2 = NGLLZ
    else if (j /= THREE) then
      i1 = 1; i2 = NGLLX
      j1 = j; j2 = j
      k1 = 1; k2 = NGLLZ
    else if (k /= THREE) then
      i1 = 1; i2 = NGLLX
      j1 = 1; j2 = NGLLY
      k1 = k; k2 = k
    else
      call exit_MPI_without_rank('incorrect wavefield discontinuity point')
    endif
  else
    call exit_MPI_without_rank('incorrect wavefield discontinuity point')
  endif
  end subroutine get_points_boundary_wd
