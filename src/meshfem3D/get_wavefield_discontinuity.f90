!! Setting up the surface on which wavefield discontinuity condition is
!! imposed
!! Tianshi Liu, 2023.5

subroutine write_wavefield_discontinuity_database(IIN_DB)
  use wavefield_discontinuity_par
  implicit none
  integer :: IIN_DB
  write(IIN_DB) nb_wd
  write(IIN_DB) boundary_to_ispec_wd
  write(IIN_DB) side_wd
  !deallocate(boundary_to_ispec_wd, side_wd)
end subroutine write_wavefield_discontinuity_database

subroutine write_wavefield_discontinuity_file()
  use wavefield_discontinuity_par
  use meshfem3d_par, only: xstore,ystore,zstore
  implicit none
  double precision :: xelm(27), yelm(27), zelm(27)
  integer :: ib, iside, ispec
  open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
       file='MESH/wavefield_discontinuity_boundary', &
       form='formatted', action='write')
  do ib = 1, nb_wd
    ispec = boundary_to_ispec_wd(ib)
    iside = side_wd(ib)
    call hex8_to_hex27(xstore(:,:,:,ispec), ystore(:,:,:,ispec), &
                       zstore(:,:,:,ispec), xelm, yelm, zelm)
    write(IFILE_WAVEFIELD_DISCONTINUITY, '(2i10,3f15.4)') ispec, iside, &
         xelm(iside), yelm(iside), zelm(iside)
  enddo
  close(IFILE_WAVEFIELD_DISCONTINUITY)
end subroutine write_wavefield_discontinuity_file

subroutine find_wavefield_discontinuity_elements()
  use constants, only: CUSTOM_REAL
  use meshfem3d_par, only: xstore,ystore,zstore,nspec
  use wavefield_discontinuity_par
  implicit none
  integer :: boundary_to_ispec_wd_temp(6*nspec), side_wd_temp(6*nspec)
  integer :: ispec, iside
  logical :: is_boundary_wd
  logical :: covered(26)
  double precision :: x_min, x_max, y_min, y_max, z_min, z_max
  double precision :: dx, dy, dz, x_mid, y_mid, z_mid, ratio_small=1.0e-6
  double precision :: xelm(27), yelm(27), zelm(27)
  !nqdu
  integer :: ier 
  integer :: myrank

  !nqdu comment, move to par_file
  ! open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
  !      file='wavefield_discontinuity_box', &
  !      form='formatted', action='read')
  ! !print*, 
  ! read(IFILE_WAVEFIELD_DISCONTINUITY, *) x_min
  ! read(IFILE_WAVEFIELD_DISCONTINUITY, *) x_max
  ! read(IFILE_WAVEFIELD_DISCONTINUITY, *) y_min
  ! read(IFILE_WAVEFIELD_DISCONTINUITY, *) y_max
  ! read(IFILE_WAVEFIELD_DISCONTINUITY, *) z_min
  ! read(IFILE_WAVEFIELD_DISCONTINUITY, *) z_max
  ! close(IFILE_WAVEFIELD_DISCONTINUITY)

  !nqdu read from par_file
  call world_rank(myrank)
  if(myrank == 0) then 
    call open_parameter_file(ier)
    call read_value_double_precision(x_min,'WAVEFIELD_DISCON_BOX_XMIN',ier)
    call read_value_double_precision(y_min,'WAVEFIELD_DISCON_BOX_YMIN',ier)
    call read_value_double_precision(z_min,'WAVEFIELD_DISCON_BOX_ZMIN',ier)
    call read_value_double_precision(x_max,'WAVEFIELD_DISCON_BOX_XMAX',ier)
    call read_value_double_precision(y_max,'WAVEFIELD_DISCON_BOX_YMAX',ier)
    call read_value_double_precision(z_max,'WAVEFIELD_DISCON_BOX_ZMAX',ier)

    ! CLOSE
    call close_parameter_file()
  endif
  call bcast_all_singledp(x_min)
  call bcast_all_singledp(y_min)
  call bcast_all_singledp(z_min)
  call bcast_all_singledp(x_max)
  call bcast_all_singledp(y_max)
  call bcast_all_singledp(z_max)

  !
  nb_wd = 0
  do ispec = 1, nspec
    covered(:) = .false.
    call hex8_to_hex27(xstore(:,:,:,ispec), ystore(:,:,:,ispec), &
                       zstore(:,:,:,ispec), xelm, yelm, zelm)
    x_mid = xelm(27)
    y_mid = yelm(27)
    z_mid = zelm(27)
    dx = ratio_small * (maxval(xelm(1:27)) - minval(xelm(1:27)))
    dy = ratio_small * (maxval(yelm(1:27)) - minval(yelm(1:27)))
    dz = ratio_small * (maxval(zelm(1:27)) - minval(zelm(1:27)))
    !! bottom
    iside = 21
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/9,10,11,12,1,2,3,4,21/)) = .true.
      endif
    endif
    !! front
    iside = 22
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/9,13,14,17,1,2,5,6,22/)) = .true.
      endif
    endif
    !! right
    iside = 23
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/10,14,15,18,2,3,6,7,23/)) = .true.
      endif
    endif
    !! back
    iside = 24
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/11,15,16,19,3,4,7,8,24/)) = .true.
      endif
    endif
    !! left
    iside = 25
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/12,13,16,20,1,4,5,8,25/)) = .true.
      endif
    endif
    !! top
    iside = 26
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/17,18,19,20,5,6,7,8,26/)) = .true.
      endif
    endif
    !! front - bottom
    iside = 9
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/1,2,9/)) = .true.
      endif
    endif
    !! right - bottom
    iside = 10
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/2,3,10/)) = .true.
      endif
    endif
    !! back - bottom
    iside = 11
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/3,4,11/)) = .true.
      endif
    endif
    !! left - bottom
    iside = 12
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/1,4,12/)) = .true.
      endif
    endif
    !! left - front
    iside = 13
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/1,5,13/)) = .true.
      endif
    endif
    !! right - front
    iside = 14
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/2,6,14/)) = .true.
      endif
    endif
    !! right - back
    iside = 15
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/3,7,15/)) = .true.
      endif
    endif
    !! left - front
    iside = 16
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/4,8,16/)) = .true.
      endif
    endif
    !! front - top
    iside = 17
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/5,6,17/)) = .true.
      endif
    endif
    !! right - top
    iside = 18
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/6,7,18/)) = .true.
      endif
    endif
    !! back - top
    iside = 19
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/7,8,19/)) = .true.
      endif
    endif
    !! left - bottom
    iside = 20
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered((/5,8,20/)) = .true.
      endif
    endif
    !! left - front - bottom
    iside = 1
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered(1) = .true.
      endif
    endif
    !! right - front - bottom
    iside = 2
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered(2) = .true.
      endif
    endif
    !! right - back - bottom
    iside = 3
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered(3) = .true.
      endif
    endif
    !! left - front - bottom
    iside = 4
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered(4) = .true.
      endif
    endif
    !! left - front - top
    iside = 5
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered(5) = .true.
      endif
    endif
    !! right - front - top
    iside = 6
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered(6) = .true.
      endif
    endif
    !! right - back - top
    iside = 7
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered(7) = .true.
      endif
    endif
    !! left - front - top
    iside = 8
    if (.not. covered(iside)) then
      if (is_boundary_wd(xelm(iside), yelm(iside), &
                       zelm(iside), x_min, x_max, y_min, y_max, &
                       z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                       IS_TOP_WAVEFIELD_DISCONTINUITY)) then
        nb_wd = nb_wd + 1
        boundary_to_ispec_wd_temp(nb_wd) = ispec
        side_wd_temp(nb_wd) = iside
        covered(8) = .true.
      endif
    endif
  enddo
  allocate(boundary_to_ispec_wd(nb_wd), side_wd(nb_wd))
  boundary_to_ispec_wd(1:nb_wd) = boundary_to_ispec_wd_temp(1:nb_wd)
  side_wd(1:nb_wd) = side_wd_temp(1:nb_wd)
end subroutine find_wavefield_discontinuity_elements

subroutine hex8_to_hex27(xstore, ystore, zstore, xelm, yelm, zelm)
  implicit none
  double precision, intent(in) :: xstore(2,2,2), ystore(2,2,2), &
                                  zstore(2,2,2)
  double precision, intent(out) :: xelm(27), yelm(27), zelm(27)
  xelm(1)  = xstore(1,1,1)
  xelm(2)  = xstore(2,1,1)
  xelm(3)  = xstore(2,2,1)
  xelm(4)  = xstore(1,2,1)
  xelm(5)  = xstore(1,1,2)
  xelm(6)  = xstore(2,1,2)
  xelm(7)  = xstore(2,2,2)
  xelm(8)  = xstore(1,2,2)
  xelm(9)  = 0.5 * sum(xstore(1:2,1,1))
  xelm(10) = 0.5 * sum(xstore(2,1:2,1))
  xelm(11) = 0.5 * sum(xstore(1:2,2,1))
  xelm(12) = 0.5 * sum(xstore(1,1:2,1))
  xelm(13) = 0.5 * sum(xstore(1,1,1:2))
  xelm(14) = 0.5 * sum(xstore(2,1,1:2))
  xelm(15) = 0.5 * sum(xstore(2,2,1:2))
  xelm(16) = 0.5 * sum(xstore(1,2,1:2))
  xelm(17) = 0.5 * sum(xstore(1:2,1,2))
  xelm(18) = 0.5 * sum(xstore(2,1:2,2))
  xelm(19) = 0.5 * sum(xstore(1:2,2,2))
  xelm(20) = 0.5 * sum(xstore(1,1:2,2))
  xelm(21) = 0.25* sum(xstore(1:2,1:2,1))
  xelm(22) = 0.25* sum(xstore(1:2,1,1:2))
  xelm(23) = 0.25* sum(xstore(2,1:2,1:2))
  xelm(24) = 0.25* sum(xstore(1:2,2,1:2))
  xelm(25) = 0.25* sum(xstore(1,1:2,1:2))
  xelm(26) = 0.25* sum(xstore(1:2,1:2,2))
  xelm(27) = 0.125*sum(xstore(1:2,1:2,1:2))

  yelm(1)  = ystore(1,1,1)
  yelm(2)  = ystore(2,1,1)
  yelm(3)  = ystore(2,2,1)
  yelm(4)  = ystore(1,2,1)
  yelm(5)  = ystore(1,1,2)
  yelm(6)  = ystore(2,1,2)
  yelm(7)  = ystore(2,2,2)
  yelm(8)  = ystore(1,2,2)
  yelm(9)  = 0.5 * sum(ystore(1:2,1,1))
  yelm(10) = 0.5 * sum(ystore(2,1:2,1))
  yelm(11) = 0.5 * sum(ystore(1:2,2,1))
  yelm(12) = 0.5 * sum(ystore(1,1:2,1))
  yelm(13) = 0.5 * sum(ystore(1,1,1:2))
  yelm(14) = 0.5 * sum(ystore(2,1,1:2))
  yelm(15) = 0.5 * sum(ystore(2,2,1:2))
  yelm(16) = 0.5 * sum(ystore(1,2,1:2))
  yelm(17) = 0.5 * sum(ystore(1:2,1,2))
  yelm(18) = 0.5 * sum(ystore(2,1:2,2))
  yelm(19) = 0.5 * sum(ystore(1:2,2,2))
  yelm(20) = 0.5 * sum(ystore(1,1:2,2))
  yelm(21) = 0.25* sum(ystore(1:2,1:2,1))
  yelm(22) = 0.25* sum(ystore(1:2,1,1:2))
  yelm(23) = 0.25* sum(ystore(2,1:2,1:2))
  yelm(24) = 0.25* sum(ystore(1:2,2,1:2))
  yelm(25) = 0.25* sum(ystore(1,1:2,1:2))
  yelm(26) = 0.25* sum(ystore(1:2,1:2,2))
  yelm(27) = 0.125*sum(ystore(1:2,1:2,1:2))

  zelm(1)  = zstore(1,1,1)
  zelm(2)  = zstore(2,1,1)
  zelm(3)  = zstore(2,2,1)
  zelm(4)  = zstore(1,2,1)
  zelm(5)  = zstore(1,1,2)
  zelm(6)  = zstore(2,1,2)
  zelm(7)  = zstore(2,2,2)
  zelm(8)  = zstore(1,2,2)
  zelm(9)  = 0.5 * sum(zstore(1:2,1,1))
  zelm(10) = 0.5 * sum(zstore(2,1:2,1))
  zelm(11) = 0.5 * sum(zstore(1:2,2,1))
  zelm(12) = 0.5 * sum(zstore(1,1:2,1))
  zelm(13) = 0.5 * sum(zstore(1,1,1:2))
  zelm(14) = 0.5 * sum(zstore(2,1,1:2))
  zelm(15) = 0.5 * sum(zstore(2,2,1:2))
  zelm(16) = 0.5 * sum(zstore(1,2,1:2))
  zelm(17) = 0.5 * sum(zstore(1:2,1,2))
  zelm(18) = 0.5 * sum(zstore(2,1:2,2))
  zelm(19) = 0.5 * sum(zstore(1:2,2,2))
  zelm(20) = 0.5 * sum(zstore(1,1:2,2))
  zelm(21) = 0.25* sum(zstore(1:2,1:2,1))
  zelm(22) = 0.25* sum(zstore(1:2,1,1:2))
  zelm(23) = 0.25* sum(zstore(2,1:2,1:2))
  zelm(24) = 0.25* sum(zstore(1:2,2,1:2))
  zelm(25) = 0.25* sum(zstore(1,1:2,1:2))
  zelm(26) = 0.25* sum(zstore(1:2,1:2,2))
  zelm(27) = 0.125*sum(zstore(1:2,1:2,1:2))
end subroutine hex8_to_hex27

logical function is_boundary_wd(x, y, z, x_min, x_max, y_min, y_max, &
                                z_min, z_max, x_mid, y_mid, z_mid, dx, dy, dz, &
                                IS_TOP_WAVEFIELD_DISCONTINUITY)
  implicit none
  double precision :: x_min, x_max, y_min, y_max, z_min, z_max
  double precision :: x, y, z, x_mid, y_mid, z_mid, dx, dy, dz
  logical :: IS_TOP_WAVEFIELD_DISCONTINUITY
  is_boundary_wd = .false.
  if (((x < x_min + dx) .or. (x > x_max - dx) .or. &
       (y < y_min + dy) .or. (y > y_max - dy) .or. &
       (z < z_min + dz) .or. ((z > z_max - dz) .and. &
       IS_TOP_WAVEFIELD_DISCONTINUITY)) .and. &
       (x_mid > x_min) .and. (x_mid < x_max) .and. &
       (y_mid > y_min) .and. (y_mid < y_max) .and. &
       (z_mid > z_min) .and. ((.not. IS_TOP_WAVEFIELD_DISCONTINUITY) .or. &
       (z_mid < z_max))) &
    is_boundary_wd = .true.
end function is_boundary_wd
