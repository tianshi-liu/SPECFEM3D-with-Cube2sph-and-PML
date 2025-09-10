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


  ! nqdu added  
  !> time step for wavefield discontinuity solver
  real(kind=CUSTOM_REAL) :: DT_wd 
  integer :: NSTEP_wd  !> number of time steps for wavefield discontinuity solver
  logical :: SAVE_DOWNSAMPLED_WD

  !> discontinuity field for all time steps
  real(kind=CUSTOM_REAL), allocatable :: &
                              field_d_wd(:,:,:), field_a_wd(:,:,:),field_t_wd(:,:,:,:)

  ! public:: lanczos_resample_to_grid
! contains 

!   !> sinc function 
!   pure function sinc(a) result(b)
!     real(kind=CUSTOM_REAL), intent(in) :: a
!     real(kind=CUSTOM_REAL) :: b
    
!     if (abs(a) > 1.0e-12_CUSTOM_REAL) then
!       b = 1.0_CUSTOM_REAL
!     else
!       b = sin(pi * a) / (pi * a)
!     end if
!   end function sinc

!   !> lanczos resampler 
!   pure real(CUSTOM_REAL) function lanczos_kernel(x, a) result(w)
    
!     use constants, only: dp => CUSTOM_REAL
!     implicit none
!     real(kind=dp), parameter :: pi = 3.14159265358979323846_dp
!     real(dp), intent(in) :: x
!     integer,  intent(in) :: a
!     real(dp) :: ax
!     ax = abs(x)
!     if (ax >= real(a, dp)) then
!       w = 0.0_dp
!     else
!       w = sinc(pi*x) * sinc(pi*x/real(a,dp))
!     end if
!   end function lanczos_kernel
end module wavefield_discontinuity_solver
