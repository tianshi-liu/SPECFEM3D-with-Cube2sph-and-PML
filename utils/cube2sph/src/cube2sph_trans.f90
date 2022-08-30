subroutine cube2sph_trans(nodes_coords,nodes_coords_sph,npts,&
                r_earth,center_lat,center_lon,rotation_azi)
  use constants, only: NDIM, RADIANS_TO_DEGREES,PI
  integer :: npts,ipt
  double precision, dimension(NDIM,npts) :: nodes_coords, nodes_coords_sph
  double precision :: r_earth,center_lat,center_lon,rotation_azi
  double precision :: &
          alpha, beta, gamma, cosa, sina, cosb, sinb, cosg, sing, t1, t2, &
          M11, M12, M13, M21, M22, M23, M31, M32, M33, xi, eta, x, y, z
  alpha = center_lon / RADIANS_TO_DEGREES
  beta = PI / 2.0 - center_lat / RADIANS_TO_DEGREES
  gamma = rotation_azi / RADIANS_TO_DEGREES
  cosa = cos(alpha)
  sina = sin(alpha)
  cosb = cos(beta)
  sinb = sin(beta)
  cosg = cos(gamma)
  sing = sin(gamma)
  M11 = cosg * cosb * cosa - sing * sina
  M12 = - sing* cosb * cosa - cosg * sina
  M13 = sinb * cosa
  M21 = cosg * cosb * sina + sing * cosa
  M22 = - sing * cosb * sina + cosg * cosa
  M23 = sinb * sina
  M31 = - cosg * sinb
  M32 = sing * sinb
  M33 = cosb
  do ipt = 1, npts
    x = nodes_coords(1, ipt)
    y = nodes_coords(2, ipt)
    z = nodes_coords(3, ipt)
    if (r_earth > 0.0) then
      xi = x / r_earth
      eta = y / r_earth
      t1 = tan(xi)
      t2 = tan(eta)
      z = (1.0+z/r_earth)/sqrt(1.0+t1*t1+t2*t2)
      x = -z*t2*r_earth
      y = z*t1*r_earth
      z = z*r_earth
    endif
    nodes_coords_sph(1, ipt) = M11 * x + M12 * y + M13 * z
    nodes_coords_sph(2, ipt) = M21 * x + M22 * y + M23 * z
    nodes_coords_sph(3, ipt) = M31 * x + M32 * y + M33 * z
  enddo
end subroutine cube2sph_trans

subroutine cube2sph_trans_inv(nodes_coords,nodes_coords_sph,npts,&
                r_earth,center_lat,center_lon,rotation_azi)
  use constants, only: NDIM, RADIANS_TO_DEGREES,PI
  integer :: npts,ipt
  double precision, dimension(NDIM,npts) :: nodes_coords, nodes_coords_sph
  double precision :: r_earth,center_lat,center_lon,rotation_azi
  double precision :: &
          alpha, beta, gamma, cosa, sina, cosb, sinb, cosg, sing, t1, t2, &
          M11, M12, M13, M21, M22, M23, M31, M32, M33, xi, eta, x, y, z
  alpha = center_lon / RADIANS_TO_DEGREES
  beta = PI / 2.0 - center_lat / RADIANS_TO_DEGREES
  gamma = rotation_azi / RADIANS_TO_DEGREES
  cosa = cos(alpha)
  sina = sin(alpha)
  cosb = cos(beta)
  sinb = sin(beta)
  cosg = cos(gamma)
  sing = sin(gamma)
  M11 = cosg * cosb * cosa - sing * sina
  M12 = - sing* cosb * cosa - cosg * sina
  M13 = sinb * cosa
  M21 = cosg * cosb * sina + sing * cosa
  M22 = - sing * cosb * sina + cosg * cosa
  M23 = sinb * sina
  M31 = - cosg * sinb
  M32 = sing * sinb
  M33 = cosb
  do ipt = 1, npts
    !x = nodes_coords(1, ipt)
    !y = nodes_coords(2, ipt)
    !z = nodes_coords(3, ipt)
    !if (r_earth > 0.0) then
    !  xi = x / r_earth
    !  eta = y / r_earth
    !  t1 = tan(xi)
    !  t2 = tan(eta)
    !  z = (1.0+z/r_earth)/sqrt(1.0+t1*t1+t2*t2)
    !  x = -z*t2*r_earth
    !  y = z*t1*r_earth
    !  z = z*r_earth
    !endif
    !nodes_coords_sph(1, ipt) = M11 * x + M12 * y + M13 * z
    !nodes_coords_sph(2, ipt) = M21 * x + M22 * y + M23 * z
    !nodes_coords_sph(3, ipt) = M31 * x + M32 * y + M33 * z

    x = M11*nodes_coords_sph(1, ipt)+M21*nodes_coords_sph(2, ipt)+&
        M31*nodes_coords_sph(3, ipt)
    y = M12*nodes_coords_sph(1, ipt)+M22*nodes_coords_sph(2, ipt)+&
        M32*nodes_coords_sph(3, ipt)
    z = M13*nodes_coords_sph(1, ipt)+M23*nodes_coords_sph(2, ipt)+&
        M33*nodes_coords_sph(3, ipt)
    z = z / r_earth
    t1 = y / z / r_earth
    t2 = - x / z / r_earth
    z = (z * sqrt(1.0+t1*t1+t2*t2) - 1.0) * r_earth
    xi = atan(t1)
    eta = atan(t2)
    nodes_coords(1, ipt) = xi * r_earth
    nodes_coords(2, ipt) = eta * r_earth
    nodes_coords(3, ipt) = z
  enddo
end subroutine cube2sph_trans_inv
