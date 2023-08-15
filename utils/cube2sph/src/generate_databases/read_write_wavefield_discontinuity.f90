subroutine read_partition_files_wavefield_discontinuity()
  use constants, only: IIN
  use wavefield_discontinuity_cube2sph, only: &
                              nb_wd, boundary_to_ispec_wd, side_wd
  implicit none
  !integer :: ier
  read(IIN) nb_wd
  allocate(boundary_to_ispec_wd(nb_wd), side_wd(nb_wd))
  read(IIN) boundary_to_ispec_wd
  read(IIN) side_wd
end subroutine read_partition_files_wavefield_discontinuity

subroutine write_partition_files_wavefield_discontinuity()
  use constants, only: IOUT
  use wavefield_discontinuity_cube2sph, only: &
                              nb_wd, boundary_to_ispec_wd, side_wd
  implicit none
  !integer :: ier
  write(IOUT) nb_wd
  write(IOUT) boundary_to_ispec_wd
  write(IOUT) side_wd
end subroutine write_partition_files_wavefield_discontinuity
