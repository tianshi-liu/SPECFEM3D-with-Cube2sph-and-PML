!! Setting up the surface on which wavefield discontinuity condition is
!! imposed
!! Tianshi Liu, 2023.5

module wavefield_discontinuity_generate_databases
  use constants, only: CUSTOM_REAL
  use wavefield_discontinuity_par
  
  !! ispec_to_elem_wd(NSPEC_AB)
  !! ispec_to_elem_wd(ispec) = ispec_wd (0 if element not belong to boundary)
  !! written in solver database and used in solver
  integer, dimension(:), allocatable :: ispec_to_elem_wd

  !! number of distinct gll points on the boundary
  !! written in solver database and used in solver
  integer :: nglob_wd

  !! number of elements on the inner side of the boundary
  !! written in solver database and used in solver
  integer :: nspec_wd

  !! ibool_wd(NGLLX, NGLLY, NGLLZ, nspec_wd)
  !! ibool_wd(i,j,k,ispec_wd) = iglob_wd (0 if point not on boundary)
  !! written in solver database and used in solver
  integer, dimension(:,:,:,:), allocatable :: ibool_wd

  !! boundary_to_iglob_wd(nglob_wd)
  !! boundary_to_iglob_wd(iglob_wd) = iglob
  !! written in solver database and used in solver
  integer, dimension(:), allocatable :: boundary_to_iglob_wd

  !! mass_in_wd(nglob_wd)
  !! mass matrix on the inner side of the boundary
  !! note that it is not assembled over processors
  !! written in solver database and used in solver
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mass_in_wd

  !! number of faces on the boundary
  !! written in solver database and used in solver
  integer :: nfaces_wd

  !! face_ijk_wd(NDIM, NGLLSQUARE, nfaces_wd)
  !! written in solver database and used in solver
  integer, dimension(:,:,:), allocatable :: face_ijk_wd

  !! face_ispec_wd(nfaces_wd)
  !! written in solver database and used in solver
  integer, dimension(:), allocatable :: face_ispec_wd
  
  !! face_normal_wd(NDIM, NGLLSQUARE, nfaces_wd)
  !! written in solver database and used in solver
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: face_normal_wd

  !! face_jacobian2Dw_wd(NGLLSQUARE, nfaces_wd)
  !! written in solver database and used in solver
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: face_jacobian2dw_wd
end module wavefield_discontinuity_generate_databases
