!! Solving the wavefield discontinuity problem with a non-split-node
!! scheme
!! Tianshi Liu, 2023.5
module wavefield_discontinuity_solver
  use constants, only: CUSTOM_REAL
  use wavefield_discontinuity_par, only: IS_WAVEFIELD_DISCONTINUITY, &
                                         IFILE_WAVEFIELD_DISCONTINUITY, &
                                         read_wavefield_discontinuity_switch
  !! ispec_to_elem_wd(NSPEC_AB)
  !! ispec_to_elem_wd(ispec) = ispec_wd (0 if element not belong to boundary)
  !! read from solver database
  integer, dimension(:), allocatable :: ispec_to_elem_wd

  !! number of distinct gll points on the boundary
  !! read from solver database
  integer :: nglob_wd

  !! number of elements on the inner side of the boundary
  !! read from solver database
  integer :: nspec_wd

  !! ibool_wd(NGLLX, NGLLY, NGLLZ, nspec_wd)
  !! ibool_wd(i,j,k,ispec_wd) = iglob_wd (0 if point not on boundary)
  !! read from solver database
  integer, dimension(:,:,:,:), allocatable :: ibool_wd

  !! boundary_to_iglob_wd(nglob_wd)
  !! boundary_to_iglob_wd(iglob_wd) = iglob
  !! read from solver database
  integer, dimension(:), allocatable :: boundary_to_iglob_wd

  !! mass_in_wd(nglob_wd)
  !! mass matrix on the inner side of the boundary
  !! note that it is not assembled over processors
  !! read from solver database
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mass_in_wd

  !! number of faces on the boundary
  !! read from solver database
  integer :: nfaces_wd

  !! face_ijk_wd(NDIM, NGLLSQUARE, nfaces_wd)
  !! read from solver database
  integer, dimension(:,:,:), allocatable :: face_ijk_wd

  !! face_ispec_wd(nfaces_wd)
  !! read from solver database
  integer, dimension(:), allocatable :: face_ispec_wd

  !! face_normal_wd(NDIM, NGLLSQUARE, nfaces_wd)
  !! read from solver database
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: face_normal_wd

  !! face_jacobian2Dw_wd(NGLLSQUARE, nfaces_wd)
  !! read from solver database
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: face_jacobian2dw_wd

  !! displ_wd(NDIM, nglob_wd)
  !! displacement discontinuity condition at current time step
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ_wd

  !! accel_wd(NDIM, nglob_wd)
  !! acceleration discontinuity condition at current time step
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_wd

  !! traction_wd(NDIM, NGLLSQUARE, nfaces_wd)
  !! traction discontinuity condition at current time step
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: traction_wd
end module wavefield_discontinuity_solver
