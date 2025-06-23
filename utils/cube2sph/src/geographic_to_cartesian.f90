subroutine geographic_to_cartesian(lat_in, lon_in, elevation, burial, &
                                   TOPOGRAPHY, ELLIPTICITY, x, y, z, nu)
  use constants_solver, only: NDIM, &
    DEGREES_TO_RADIANS,R_UNIT_SPHERE, &
    R_EARTH,R_EARTH_KM,NX_BATHY,NY_BATHY
  use specfem_par, only: &
    ONE_CRUST, rspl, espl, espl2, nspl, ibathy_topo

  implicit none
  logical, intent(in) :: TOPOGRAPHY, ELLIPTICITY
  ! elevation, burial in meters
  double precision, intent(in) :: lat_in, lon_in, elevation, burial
  double precision, intent(out) :: x, y, z, nu(NDIM, NDIM)
  !! local variables
  integer :: iorientation 
  double precision :: stazi, stdip, r, r0, topo
  double precision :: lat, lon, theta, phi, thetan, phin, n(3), &
                      sint, cost, sinp, cosp
  double precision :: ell, p20

  lat = lat_in
  lon = lon_in
  ONE_CRUST = .true.
  ! limits longitude to [0.0,360.0]
  if (lon < 0.d0 ) lon = lon + 360.d0
  if (lon > 360.d0 ) lon = lon - 360.d0

  ! converts geographic latitude stlat (degrees) to geocentric colatitude
  ! theta (radians)
  call lat_2_geocentric_colat_dble(lat,theta)

  phi = lon*DEGREES_TO_RADIANS
  call reduce(theta,phi)

  ! record three components for each station
  do iorientation = 1,3
    !     North
    if (iorientation == 1) then
      stazi = 0.d0
      stdip = 0.d0
    !     East
    else if (iorientation == 2) then
      stazi = 90.d0
      stdip = 0.d0
    !     Vertical
    else if (iorientation == 3) then
      stazi = 0.d0
      stdip = - 90.d0
    endif

    !     get the orientation of the seismometer
    thetan=(90.0d0+stdip)*DEGREES_TO_RADIANS
    phin=stazi*DEGREES_TO_RADIANS

    ! we use the same convention as in Harvard normal modes for the
    ! orientation

    !     vertical component
    n(1) = cos(thetan)
    !     N-S component
    n(2) = - sin(thetan)*cos(phin)
    !     E-W component
    n(3) = sin(thetan)*sin(phin)

    !     get the Cartesian components of n in the model: nu
    sint = sin(theta)
    cost = cos(theta)
    sinp = sin(phi)
    cosp = cos(phi)
    nu(iorientation,1) = n(1)*sint*cosp+n(2)*cost*cosp-n(3)*sinp
    nu(iorientation,2) = n(1)*sint*sinp+n(2)*cost*sinp+n(3)*cosp
    nu(iorientation,3) = n(1)*cost-n(2)*sint
  enddo 
  
  ! normalized receiver radius
  r0 = R_UNIT_SPHERE

  ! finds elevation of receiver
  topo = 0.0
  if (TOPOGRAPHY) then
     call get_topo_bathy(lat,lon,topo,ibathy_topo)
  endif
  !print *, topo, elevation
  r0 = r0 + (elevation + topo)/R_EARTH
  ! ellipticity
  if (ELLIPTICITY) then
    cost=cos(theta)
    ! this is the Legendre polynomial of degree two, P2(cos(theta)),
    ! see the discussion above eq (14.4) in Dahlen and Tromp (1998)
    p20=0.5d0*(3.0d0*cost*cost-1.0d0)
    ! get ellipticity using spline evaluation
    call spline_evaluation(rspl,espl,espl2,nspl,r0,ell)
    ! this is eq (14.4) in Dahlen and Tromp (1998)
    r0=r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
  endif

  ! subtract station burial depth (in meters)
  r0 = r0 - burial/R_EARTH

  ! compute the Cartesian position of the receiver
  x = r0*sin(theta)*cos(phi)
  y = r0*sin(theta)*sin(phi)
  z = r0*cos(theta)

  x = x * R_EARTH
  y = y * R_EARTH
  z = z * R_EARTH 

end subroutine geographic_to_cartesian

subroutine prepare_topography_ellipticity(TOPOGRAPHY)
  use constants_solver, only: &
    CUSTOM_REAL, NX_BATHY, NY_BATHY
  use specfem_par, only: &
    ONE_CRUST, ibathy_topo, nspl, rspl, espl, espl2
  implicit none
  !! local variable
  logical, intent(in):: TOPOGRAPHY
  integer :: ier
  ONE_CRUST = .true.
  call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  if (TOPOGRAPHY) then
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    ibathy_topo(:,:) = 0
    call read_topo_bathy_file(ibathy_topo)
  endif
end subroutine prepare_topography_ellipticity

subroutine deallocate_topography()
  use specfem_par, only: ibathy_topo
  implicit none
  deallocate(ibathy_topo)
end subroutine deallocate_topography
