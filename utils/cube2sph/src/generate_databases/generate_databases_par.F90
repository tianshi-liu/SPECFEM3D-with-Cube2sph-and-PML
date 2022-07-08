!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  module generate_databases_par

!  use constants, only: NGLLX,NGLLY,NGLLZ,NGLLSQUARE,NDIM,NDIM2D,NGNOD2D_FOUR_CORNERS,N_SLS, &
!    CUSTOM_REAL,SIZE_REAL,SIZE_DOUBLE, &
!    IMAIN,IIN,IOUT,ISTANDARD_OUTPUT, &
!    ZERO,ONE,TWO,FOUR_THIRDS,PI,HUGEVAL,GAUSSALPHA,GAUSSBETA, &
!    SMALLVAL_TOL,TINYVAL,HUGEVAL,R_EARTH, &
!    MAX_STRING_LEN,ATTENUATION_COMP_MAXIMUM, &
!    MINIMUM_THICKNESS_3D_OCEANS,RHO_APPROXIMATE_OCEAN_LOAD, &
!    CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
!    CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ, &
!    NPOWER,CPML_Rcoef, &
!    IMODEL_DEFAULT,IMODEL_GLL,IMODEL_1D_PREM,IMODEL_1D_CASCADIA,IMODEL_1D_SOCAL, &
!    IMODEL_SALTON_TROUGH,IMODEL_TOMO,IMODEL_USER_EXTERNAL,IMODEL_IPATI,IMODEL_IPATI_WATER, &
!    IMODEL_1D_PREM_PB,IMODEL_SEP,IMODEL_COUPLED, &
!    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC, &
!    OUTPUT_FILES, &
!    NX_TOPO_FILE,NY_TOPO_FILE, &
!    USE_MESH_COLORING_GPU,MAX_NUMBER_OF_COLORS, &
!    ADIOS_TRANSPORT_METHOD

!  use shared_parameters
  use constants, only: MAX_STRING_LEN, CUSTOM_REAL
  implicit none

! number of spectral elements in each block
!  integer npointot

! local to global indexing array
!  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
!  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
!  integer :: myrank,sizeprocs,ier

!  integer :: NX_TOPO,NY_TOPO
!  integer, dimension(:,:), allocatable :: itopo_bathy

! timer MPI
!  double precision :: time_start,tCPU

! memory size that will be needed by the solver
!  double precision :: max_memory_size,max_memory_size_request
  logical :: SAVE_MOHO_MESH

! this for all the regions
  integer NSPEC_AB,NGLOB_AB

  integer NSPEC2D_BOTTOM,NSPEC2D_TOP

!  double precision min_elevation,max_elevation
!  double precision min_elevation_all,max_elevation_all

! for Databases of external meshes
  double precision, dimension(:,:), allocatable :: nodes_coords

!  integer :: dummy_node
!  integer :: dummy_elmnt

  integer :: ispec, inode, num_interface,ie,imat,iface,icorner
  integer :: npts, nelmnts_ext_mesh
  integer  :: num_interfaces_ext_mesh
  integer  :: max_interface_size_ext_mesh
  integer  :: nmat_ext_mesh, nundefMat_ext_mesh
  integer, dimension(:), allocatable  :: my_neighbors_ext_mesh
  integer, dimension(:), allocatable  :: my_nelmnts_neighbors_ext_mesh

  integer, dimension(:,:,:), allocatable  :: my_interfaces_ext_mesh
  integer, dimension(:,:), allocatable  :: ibool_interfaces_ext_mesh
  integer, dimension(:), allocatable  :: nibool_interfaces_ext_mesh

  integer, dimension(:,:), allocatable :: elmnts_ext_mesh
  integer, dimension(:,:), allocatable :: mat_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh

  character(len=MAX_STRING_LEN) :: prname

! boundaries and materials
  double precision, dimension(:,:), allocatable :: materials_ext_mesh

  integer :: ispec2D, boundary_number
  integer :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, nspec2D_bottom_ext, nspec2D_top_ext

  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax, &
              ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top
  integer, dimension(:,:), allocatable :: nodes_ibelm_xmin,nodes_ibelm_xmax, &
              nodes_ibelm_ymin, nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top

  character(len=MAX_STRING_LEN), dimension(:,:), allocatable :: undef_mat_prop

! C-PML absorbing boundary conditions

  ! local number of C-PML spectral elements
  integer :: nspec_cpml

  ! global number of C-PML spectral elements
  integer :: nspec_cpml_tot

  ! C-PML spectral elements global indexing
  integer, dimension(:), allocatable :: CPML_to_spec

  ! C-PML regions
  integer, dimension(:), allocatable :: CPML_regions

  ! mask of C-PML elements for the global mesh
  logical, dimension(:), allocatable :: is_CPML

  ! thickness of C-PML layers in each direction
  real(kind=CUSTOM_REAL) :: CPML_width_x,CPML_width_y,CPML_width_z

  ! C-PML damping profile arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: d_store_x, d_store_y, d_store_z

  ! auxiliary parameters arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: K_store_x, K_store_y, K_store_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: alpha_store_x,alpha_store_y,alpha_store_z

  ! minimum distance between parameters of CPML to avoid the singularities
  real(kind=CUSTOM_REAL) :: min_distance_between_CPML_parameter

  ! array recording the points on interface shared by PML and interior computational domain
  logical, dimension(:), allocatable :: mask_ibool_interior_domain
  integer :: nglob_interface_PML_acoustic,nglob_interface_PML_elastic
  integer, dimension(:), allocatable :: points_interface_PML_acoustic, points_interface_PML_elastic
  
  !! TL: add these variables for ADE-PML
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: pml_d,&
          pml_kappa,pml_beta
  integer :: num_pml_physical
  integer, dimension(:), allocatable :: pml_physical_ispec
  integer, dimension(:,:,:), allocatable :: pml_physical_ijk
  integer, dimension(:,:,:,:), allocatable :: ibool_CPML
  integer, dimension(:), allocatable :: CPML_to_glob
  integer :: nglob_CPML
  integer :: num_interfaces_PML,max_nibool_interfaces_PML
  integer, dimension(:), allocatable :: nibool_interfaces_PML,my_neighbors_PML
  integer, dimension(:,:), allocatable :: ibool_interfaces_PML
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! moho (optional)
  integer :: nspec2D_moho_ext
  integer, dimension(:), allocatable  :: ibelm_moho
  integer, dimension(:,:), allocatable  :: nodes_ibelm_moho

  integer :: nspec_total
  integer :: nspec_irregular

! this can overflow if more than 2 Gigapoints in the whole mesh, thus replaced with double precision version
! integer :: nglob_total
  double precision :: nglob_total

! mesh surface
  logical,dimension(:),allocatable :: ispec_is_surface_external_mesh,iglob_is_surface_external_mesh
  integer :: nfaces_surface,nfaces_surface_glob_ext_mesh
  integer :: NPROC

  end module generate_databases_par

!
!-------------------------------------------------------------------------------------------------
!

