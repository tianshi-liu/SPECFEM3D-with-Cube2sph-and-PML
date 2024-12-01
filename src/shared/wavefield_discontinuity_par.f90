!! Setting up the surface on which wavefield discontinuity condition is
!! imposed
!! Tianshi Liu, 2023.5

!nqdu move everything in Par_file
module wavefield_discontinuity_par
  logical :: IS_WAVEFIELD_DISCONTINUITY=.false.

  logical :: IS_TOP_WAVEFIELD_DISCONTINUITY=.false.

  integer, parameter :: IFILE_WAVEFIELD_DISCONTINUITY=527

  !! boundary of wavefield discontinuity, read from database file
  integer :: nb_wd

  !! boundary_to_ispec_wd(nb_wd)
  !! the element the boundary belongs to, read from database file
  !! each point on the boundary belongs to two sides of the boundary
  !! here the element must be on the inner side of the boundary
  integer, dimension(:), allocatable :: boundary_to_ispec_wd

  !! side_wd(nb_wd)
  !! integers specifying which side the boundary is in the element
  !! read from database file
  !! side_wd = 1--8: only one vertex is on the boundary
  !! side_wd = 9--20: only one edge is on the boundary
  !! side_wd = 21--26: one face is on the boundary
  integer, dimension(:), allocatable :: side_wd

  contains

  subroutine read_wavefield_discontinuity_switch()
    use constants, only: INJECTION_TECHNIQUE_IS_WAVEDISCON
    implicit none
    integer :: ier
    integer :: injection_tech_type 
    
    ! opens file Par_file
    call open_parameter_file(ier)
    
    ! read INJECTION_TECH_TYPE
    call read_value_integer(injection_tech_type,'INJECTION_TECHNIQUE_TYPE',ier)
    if (injection_tech_type == INJECTION_TECHNIQUE_IS_WAVEDISCON) &
        IS_WAVEFIELD_DISCONTINUITY = .true.
    
    if(IS_WAVEFIELD_DISCONTINUITY) then 
      ! IS_TOP_WAVEFIELD_DISCONTINUITY
      call read_value_logical(IS_TOP_WAVEFIELD_DISCONTINUITY,'WAVEFIELD_DISCON_AT_TOP',ier)
      if(ier /=0) then 
        print*,'WAVEFIELD_DISCON_AT_TOP is not in in Par_file, set = ',IS_TOP_WAVEFIELD_DISCONTINUITY
      endif
    endif

    ! close file 
    call close_parameter_file()

  end subroutine read_wavefield_discontinuity_switch
end module wavefield_discontinuity_par


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! original one 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module wavefield_discontinuity_par
!   logical :: IS_WAVEFIELD_DISCONTINUITY=.false.

!   logical :: IS_TOP_WAVEFIELD_DISCONTINUITY=.true.

!   integer, parameter :: IFILE_WAVEFIELD_DISCONTINUITY=527

!   !! boundary of wavefield discontinuity, read from database file
!   integer :: nb_wd

!   !! boundary_to_ispec_wd(nb_wd)
!   !! the element the boundary belongs to, read from database file
!   !! each point on the boundary belongs to two sides of the boundary
!   !! here the element must be on the inner side of the boundary
!   integer, dimension(:), allocatable :: boundary_to_ispec_wd

!   !! side_wd(nb_wd)
!   !! integers specifying which side the boundary is in the element
!   !! read from database file
!   !! side_wd = 1--8: only one vertex is on the boundary
!   !! side_wd = 9--20: only one edge is on the boundary
!   !! side_wd = 21--26: one face is on the boundary
!   integer, dimension(:), allocatable :: side_wd

!   contains

!   subroutine read_wavefield_discontinuity_switch()
!     implicit none
!     integer :: ier
!     !integer :: myrank
!     !call world_rank(myrank)
!     open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
!        file='wavefield_discontinuity_switch', &
!        form='formatted', action='read', iostat=ier)
!     if (ier /= 0) then
!       !if (myrank == 0) print *, 'cannot find switch file for wavefield discontinuity, skip'
!       !print *, 'cannot find switch file for wavefield discontinuity, skip'
!       !print*,IS_TOP_WAVEFIELD_DISCONTINUITY
!       return 
!     else
!       read(IFILE_WAVEFIELD_DISCONTINUITY, *) IS_WAVEFIELD_DISCONTINUITY
!       read(IFILE_WAVEFIELD_DISCONTINUITY, *) IS_TOP_WAVEFIELD_DISCONTINUITY
!       close(IFILE_WAVEFIELD_DISCONTINUITY)
!       !if (myrank == 0) then
!       !  print *, 'found switch file for wavefield discontinuity'
!       !  print *, 'IS_WAVEFIELD_DISCONTINUITY = ', IS_WAVEFIELD_DISCONTINUITY
!       !  print *, 'IS_TOP_WAVEFIELD_DISCONTINUITY = ', IS_TOP_WAVEFIELD_DISCONTINUITY
!       !endif
!     endif
!   end subroutine read_wavefield_discontinuity_switch
! end module wavefield_discontinuity_par
