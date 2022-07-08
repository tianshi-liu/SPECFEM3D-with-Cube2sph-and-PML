
subroutine rotate_aniso_tensor(theta,phi,d11,d12,d13,d14,d15,d16, &
                           d22,d23,d24,d25,d26, &
                           d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                           c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  implicit none

  double precision theta,phi
  double precision c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                   c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  double precision d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                   d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision costheta,sintheta,cosphi,sinphi
  double precision costhetasq,sinthetasq,cosphisq,sinphisq
  double precision costwotheta,sintwotheta,costwophi,sintwophi
  double precision cosfourtheta,sinfourtheta
  double precision costhetafour,sinthetafour,cosphifour,sinphifour
  double precision sintwophisq,sintwothetasq

  costheta = dcos(theta)
  sintheta = dsin(theta)
  cosphi = dcos(phi)
  sinphi = dsin(phi)

  costhetasq = costheta * costheta
  sinthetasq = sintheta * sintheta
  cosphisq = cosphi * cosphi
  sinphisq = sinphi * sinphi

  costhetafour = costhetasq * costhetasq
  sinthetafour = sinthetasq * sinthetasq
  cosphifour = cosphisq * cosphisq
  sinphifour = sinphisq * sinphisq

  costwotheta = dcos(2.d0*theta)
  sintwotheta = dsin(2.d0*theta)
  costwophi = dcos(2.d0*phi)
  sintwophi = dsin(2.d0*phi)

  cosfourtheta = dcos(4.d0*theta)
  sinfourtheta = dsin(4.d0*theta)
  sintwothetasq = sintwotheta * sintwotheta
  sintwophisq = sintwophi * sintwophi

  ! recompute 21 anisotropic coefficients for full anisotropic model using Mathematica

  c11 = d22*sinphifour - 2.*sintwophi*sinphisq*(d26*costheta + d24*sintheta) - &
      2.*cosphisq*sintwophi*(d16*costhetasq*costheta + &
      (d14 + 2*d56)*costhetasq*sintheta + &
      (d36 + 2*d45)*costheta*sinthetasq + d34*sintheta*sinthetasq) + &
      cosphifour*(d11*costhetafour + 2.*d15*costhetasq*sintwotheta + &
      (d13 + 2.*d55)*sintwothetasq/2. + &
      2.*d35*sintwotheta*sinthetasq + d33*sinthetafour) + &
      (sintwophisq/4.)*(d12 + d23 + 2.*(d44 + d66) + &
      (d12 - d23 - 2.*d44 + 2.*d66)*costwotheta + &
      2.*(d25 + 2.*d46)*sintwotheta)

  c12 = -((sintwophi/2.)*sinphisq*((3.*d16 - 4.*d26 + d36 + 2.*d45)*costheta + &
      (d16 - d36 - 2.*d45)*(4.*costhetasq*costheta - 3.*costheta) + &
      2.*(d14 - 2.*d24 + d34 + 2.*d56 + &
      (d14 - d34 + 2.*d56)*costwotheta)*sintheta))/2. + &
      cosphisq*sintwophi*(d16*costhetasq*costheta - d24*sintheta + &
      (d14 + 2.*d56)*costhetasq*sintheta + d34*sintheta*sinthetasq + &
      costheta*(-d26 + (d36 + 2.*d45)*sinthetasq)) + &
      (sintwophisq/4.)*(d22 + d11*costhetafour + &
      2.*d15*costhetasq*sintwotheta - 4.*d44*sinthetasq + &
      d33*sinthetafour + costhetasq*(-4.*d66 + &
      2.*(d13 + 2.*d55)*sinthetasq) + &
      costheta*(-8.*d46*sintheta + 4.*d35*sintheta*sinthetasq)) + &
      (cosphifour + sinphifour)*(d12*costhetasq + &
      d23*sinthetasq + d25*sintwotheta)

  c13 = sinphisq*(d23*costhetasq - d25*sintwotheta + d12*sinthetasq) - &
      sintwophi*(d36*costhetasq*costheta + &
      (d34 - 2.*d56)*costhetasq*sintheta + &
      (d16 - 2.*d45)*costheta*sinthetasq + d14*sintheta*sinthetasq) + &
      (cosphisq*(d11 + 6.*d13 + d33 - 4.*d55 - &
      (d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta + &
      4.*(-d15 + d35)*sinfourtheta))/8.

  c14 = (-4.*cosphi*sinphisq*((-d14 - 2.*d24 + d34 + 2.*d56)*costheta + &
      (d14 - d34 + 2.*d56)*(4.*costhetasq*costheta - 3.*costheta) + &
      2.*(-d16 + d26 + d36 + (-d16 + d36 + 2.*d45)*costwotheta)*sintheta) + &
      8.*cosphisq*cosphi*(d14*costhetasq*costheta - &
      (d16 - 2.*d45)*costhetasq*sintheta + &
      (d34 - 2.*d56)*costheta*sinthetasq - d36*sintheta*sinthetasq) + &
      4.*sinphi*sinphisq*(2.*d25*costwotheta + (-d12 + d23)*sintwotheta) + &
      cosphisq*sinphi*(4.*(d15 + d35 - 4*d46)*costwotheta + &
      4.*(d15 - d35)*cosfourtheta - &
      2.*(d11 - d33 + 4.*d44 - 4.*d66 + &
      (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sintwotheta))/8.

  c15 = (8.*sinphi*sinphisq*(-(d24*costheta) + d26*sintheta) + &
      4.*cosphi*sinphisq*(2.*(d25 + 2.*d46)*costwotheta + &
      (-d12 + d23 + 2.*d44 - 2.*d66)*sintwotheta) + &
      cosphisq*cosphi*(4.*(d15 + d35)*costwotheta + &
      4.*(d15 - d35)*cosfourtheta - 2.*(d11 - d33 + &
      (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sintwotheta) - &
      2.*cosphisq*sinphi*((d14 + 3.*d34 + 2.*d56)*costheta + &
      3.*(d14 - d34 + 2.*d56)*(4.*costhetasq*costheta - 3.*costheta) - &
      (3.*d16 + d36 + 2.*d45)*sintheta + &
      3.*(-d16 + d36 + 2.*d45)*(-4.*sinthetasq*sintheta + 3.*sintheta)))/8.

  c16 = -(sinphifour*(d26*costheta + d24*sintheta)) - &
      (3.*(sintwophisq/4.)*((3.*d16 - 4.*d26 + d36 + 2.*d45)*costheta + &
      (d16 - d36 - 2.*d45)*(4.*costhetasq*costheta - 3.*costheta) + &
      2.*(d14 - 2.*d24 + d34 + 2.*d56 + &
      (d14 - d34 + 2.*d56)*costwotheta)*sintheta))/4. + &
      cosphifour*(d16*costhetasq*costheta + &
      (d14 + 2.*d56)*costhetasq*sintheta + &
      (d36 + 2.*d45)*costheta*sinthetasq + d34*sintheta*sinthetasq) + &
      (sintwophi/2.)*sinphisq*(-d22 + (d12 + 2.*d66)*costhetasq + &
      2.*d46*sintwotheta + (d23 + 2.*d44)*sinthetasq + d25*sintwotheta) + &
      cosphisq*(sintwophi/2.)*(d11*costhetafour + &
      2.*d15*costhetasq*sintwotheta - (d23 + 2.*d44)*sinthetasq + &
      d33*sinthetafour - costhetasq*(d12 + &
      2.*d66 - 2.*(d13 + 2.*d55)*sinthetasq) - &
      (d25 - d35 + 2.*d46 + d35*costwotheta)*sintwotheta)

  c22 = d22*cosphifour + 2.*cosphisq*sintwophi*(d26*costheta + d24*sintheta) + &
      2.*sintwophi*sinphisq*(d16*costhetasq*costheta + &
      (d14 + 2.*d56)*costhetasq*sintheta + &
      (d36 + 2.*d45)*costheta*sinthetasq + d34*sintheta*sinthetasq) + &
      sinphifour*(d11*costhetafour + 2.*d15*costhetasq*sintwotheta + &
      (d13 + 2.*d55)*sintwothetasq/2. + &
      2.*d35*sintwotheta*sinthetasq + d33*sinthetafour) + &
      (sintwophisq/4.)*(d12 + d23 + 2.*(d44 + d66) + &
      (d12 - d23 - 2.*d44 + 2.*d66)*costwotheta + &
      2.*(d25 + 2.*d46)*sintwotheta)

  c23 = d13*costhetafour*sinphisq + &
      sintheta*sinthetasq*(d14*sintwophi + d13*sinphisq*sintheta) + &
      costheta*sinthetasq*((d16 - 2.*d45)*sintwophi + &
      2.*(d15 - d35)*sinphisq*sintheta) + &
      costhetasq*costheta*(d36*sintwophi + &
      2.*(-d15 + d35)*sinphisq*sintheta) + &
      costhetasq*sintheta*((d34 - 2.*d56)*sintwophi + &
      (d11 + d33 - 4.*d55)*sinphisq*sintheta) + &
      cosphisq*(d23*costhetasq - d25*sintwotheta + d12*sinthetasq)

  c24 = (8.*cosphisq*cosphi*(d24*costheta - d26*sintheta) + &
      4.*cosphisq*sinphi*(2.*(d25 + 2.*d46)*costwotheta + &
      (-d12 + d23 + 2.*d44 - 2.*d66)*sintwotheta) + &
      sinphi*sinphisq*(4.*(d15 + d35)*costwotheta + &
      4.*(d15 - d35)*cosfourtheta - &
      2.*(d11 - d33 + (d11 - 2.*d13 + &
      d33 - 4.*d55)*costwotheta)*sintwotheta) + &
      2.*cosphi*sinphisq*((d14 + 3.*d34 + 2.*d56)*costheta + &
      3.*(d14 - d34 + 2.*d56)*(4.*costhetasq*costheta - 3.*costheta) - &
      (3.*d16 + d36 + 2.*d45)*sintheta + &
      3.*(-d16 + d36 + 2.*d45)*(-4.*sinthetasq*sintheta + 3.*sintheta)))/8.

  c25 = (4.*cosphisq*sinphi*((-d14 - 2.*d24 + d34 + 2.*d56)*costheta + &
      (d14 - d34 + 2.*d56)*(4.*costhetasq*costheta - 3.*costheta) + &
      2.*(-d16 + d26 + d36 + (-d16 + d36 + 2.*d45)*costwotheta)*sintheta) - &
      8.*sinphi*sinphisq*(d14*costhetasq*costheta - &
      (d16 - 2.*d45)*costhetasq*sintheta + &
      (d34 - 2.*d56)*costheta*sinthetasq - d36*sintheta*sinthetasq) + &
      4.*cosphisq*cosphi*(2.*d25*costwotheta + (-d12 + d23)*sintwotheta) + &
      cosphi*sinphisq*(4.*(d15 + d35 - 4.*d46)*costwotheta + &
      4.*(d15 - d35)*cosfourtheta - 2.*(d11 - d33 + 4.*d44 - 4.*d66 + &
      (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sintwotheta))/8.

  c26 = cosphifour*(d26*costheta + d24*sintheta) + &
      (3.*(sintwophisq/4.)*((3.*d16 - 4.*d26 + d36 + 2.*d45)*costheta + &
      (d16 - d36 - 2.*d45)*(4.*costhetasq*costheta - 3.*costheta) + &
      2.*(d14 - 2.*d24 + d34 + 2.*d56 + &
      (d14 - d34 + 2.*d56)*costwotheta)*sintheta))/4. - &
      sinphifour*(d16*costhetasq*costheta + &
      (d14 + 2.*d56)*costhetasq*sintheta + &
      (d36 + 2.*d45)*costheta*sinthetasq + d34*sintheta*sinthetasq) + &
      cosphisq*(sintwophi/2.)*(-d22 + (d12 + 2.*d66)*costhetasq + &
      2.*d46*sintwotheta + (d23 + 2.*d44)*sinthetasq + &
      d25*sintwotheta) + (sintwophi/2.)*sinphisq*(d11*costhetafour + &
      2.*d15*costhetasq*sintwotheta - (d23 + 2.*d44)*sinthetasq + &
      d33*sinthetafour - costhetasq*(d12 + &
      2.*d66 - 2.*(d13 + 2.*d55)*sinthetasq) - &
      (d25 - d35 + 2.*d46 + d35*costwotheta)*sintwotheta)

  c33 = d33*costhetafour - 2.*d35*costhetasq*sintwotheta + &
      (d13 + 2.*d55)*sintwothetasq/2. - &
      2.*d15*sintwotheta*sinthetasq + d11*sinthetafour

  c34 = cosphi*(d34*costhetasq*costheta - (d36 + 2.*d45)*costhetasq*sintheta + &
      (d14 + 2.*d56)*costheta*sinthetasq - d16*sintheta*sinthetasq) + &
      (sinphi*(4.*(d15 + d35)*costwotheta + 4.*(-d15 + d35)*cosfourtheta + &
      2.*(-d11 + d33)*sintwotheta + &
      (d11 - 2.*d13 + d33 - 4.*d55)*sinfourtheta))/8.

  c35 = sinphi*(-(d34*costhetasq*costheta) + &
      (d36 + 2.*d45)*costhetasq*sintheta - &
      (d14 + 2.*d56)*costheta*sinthetasq + d16*sintheta*sinthetasq) + &
      (cosphi*(4.*(d15 + d35)*costwotheta + 4.*(-d15 + d35)*cosfourtheta + &
      2.*(-d11 + d33)*sintwotheta + &
      (d11 - 2.*d13 + d33 - 4.*d55)*sinfourtheta))/8.

  c36 = (4.*costwophi*((d16 + 3.*d36 - 2.*d45)*costheta + &
      (-d16 + d36 + 2.*d45)*(4.*costhetasq*costheta - 3.*costheta) + &
      (3.*d14 + d34 - 2.*d56)*sintheta + &
      (-d14 + d34 - 2.*d56)*(-4.*sinthetasq*sintheta + 3.*sintheta)) + &
      sintwophi*(d11 - 4.*d12 + 6.*d13 - 4.*d23 + d33 - 4.*d55 + &
      4.*(d12 - d23)*costwotheta - &
      (d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta + &
      8.*d25*sintwotheta + 4.*(-d15 + d35)*sinfourtheta))/16.

  c44 = (d11 - 2.*d13 + d33 + 4.*(d44 + d55 + d66) - &
      (d11 - 2.*d13 + d33 - 4.*(d44 - d55 + d66))*costwophi + &
      4.*sintwophi*((d16 - d36 + 2.*d45)*costheta + &
      (-d16 + d36 + 2.*d45)*(4.*costhetasq*costheta - 3.*costheta) - &
      2.*(d14 - d34 + (d14 - d34 + 2.*d56)*costwotheta)*sintheta) + &
      8.*cosphisq*((d44 - d66)*costwotheta - 2.*d46*sintwotheta) + &
      2.*sinphisq*(-((d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta) + &
      4.*(-d15 + d35)*sinfourtheta))/16.

  c45 = (4.*costwophi*((d16 - d36 + 2.*d45)*costheta + &
      (-d16 + d36 + 2.*d45)*(4.*costhetasq*costheta - 3.*costheta) - &
      2.*(d14 - d34 + (d14 - d34 + 2.*d56)*costwotheta)*sintheta) + &
      sintwophi*(d11 - 2.*d13 + d33 - 4.*(d44 - d55 + d66) + &
      4.*(-d44 + d66)*costwotheta - &
      (d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta + 8.*d46*sintwotheta + &
      4.*(-d15 + d35)*sinfourtheta))/16.

  c46 = (-2.*sinphi*sinphisq*((-d14 + d34 + 2.*d56)*costheta + &
      (d14 - d34 + 2.*d56)*(4.*costhetasq*costheta - 3.*costheta) + &
      2.*(-d16 + d36 + (-d16 + d36 + 2.*d45)*costwotheta)*sintheta) + &
      4.*cosphisq*cosphi*(2.*d46*costwotheta + (d44 - d66)*sintwotheta) + &
      cosphi*sinphisq*(4.*(d15 - 2.*d25 + d35 - 2.*d46)*costwotheta + &
      4.*(d15 - d35)*cosfourtheta - &
      2.*(d11 - 2.*d12 + 2.*d23 - d33 + 2.*d44 - 2.*d66 + &
      (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sintwotheta) + &
      4.*cosphisq*sinphi*((d14 - 2.*d24 + d34)*costheta + &
      (d14 - d34 + 2.*d56)*(4.*costhetasq*costheta - 3.*costheta) - &
      (d16 - 2.*d26 + d36)*sintheta + &
      (-d16 + d36 + 2.*d45)*(-4.*sinthetasq*sintheta + 3.*sintheta)))/8.

  c55 = d66*sinphisq*sinthetasq + (sintwotheta/2.)*(-2.*d46*sinphisq + &
      (d36 + d45)*sintwophi*sintheta) + &
      costhetasq*(d44*sinphisq + (d14 + d56)*sintwophi*sintheta) - &
      sintwophi*(d45*costhetasq*costheta + d34*costhetasq*sintheta + &
      d16*costheta*sinthetasq + d56*sintheta*sinthetasq) + &
      (cosphisq*(d11 - 2.*d13 + d33 + 4.*d55 - &
      (d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta + &
      4.*(-d15 + d35)*sinfourtheta))/8.

  c56 = (8.*cosphisq*cosphi*(d56*costhetasq*costheta - &
      (d16 - d36 - d45)*costhetasq*sintheta - &
      (d14 - d34 + d56)*costheta*sinthetasq - d45*sintheta*sinthetasq) + &
      4.*sinphi*sinphisq*(2.*d46*costwotheta + (d44 - d66)*sintwotheta) + &
      cosphisq*sinphi*(4.*(d15 - 2.*d25 + d35 - 2.*d46)*costwotheta + &
      4.*(d15 - d35)*cosfourtheta - &
      2.*(d11 - 2.*d12 + 2.*d23 - d33 + 2.*d44 - 2.*d66 + &
      (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sintwotheta) - &
      4.*cosphi*sinphisq*((d14 - 2.*d24 + d34)*costheta + &
      (d14 - d34 + 2.*d56)*(4.*costhetasq*costheta - 3.*costheta) - &
      (d16 - 2.*d26 + d36)*sintheta + &
      (-d16 + d36 + 2.*d45)*(-4.*sinthetasq*sintheta + 3.*sintheta)))/8.

  c66 = -((sintwophi/2.)*sinphisq*((3.*d16 - 4.*d26 + d36 + 2.*d45)*costheta + &
      (d16 - d36 - 2.*d45)*(4.*costhetasq*costheta - 3.*costheta) + &
      2.*(d14 - 2.*d24 + d34 + 2.*d56 + &
      (d14 - d34 + 2.*d56)*costwotheta)*sintheta))/2. + &
      cosphisq*sintwophi*(d16*costhetasq*costheta - d24*sintheta + &
      (d14 + 2.*d56)*costhetasq*sintheta + d34*sintheta*sinthetasq + &
      costheta*(-d26 + (d36 + 2.*d45)*sinthetasq)) + &
      (sintwophisq/4.)*(d22 + d11*costhetafour + &
      2.*d15*costhetasq*sintwotheta - 2.*(d23 + d44)*sinthetasq + &
      d33*sinthetafour - 2.*sintwotheta*(d25 + d46 - d35*sinthetasq) - &
      2.*costhetasq*(d12 + d66 - (d13 + 2.*d55)*sinthetasq)) + &
      (cosphifour + sinphifour)*(d66*costhetasq + &
      d44*sinthetasq + d46*sintwotheta)

end subroutine rotate_aniso_tensor



subroutine rotate_kernels_dble(cij_kl,cij_kll,theta_in,phi_in)

! Purpose : compute the kernels in r,theta,phi (cij_kll)
! from the kernels in x,y,z (cij_kl) (x,y,z to r,theta,phi)
! At r,theta,phi fixed
! theta and phi are in radians

! Coeff from Min's routine rotate_anisotropic_tensor
! with the help of Collect[Expand[cij],{dij}] in Mathematica

! Definition of the output array cij_kll :
! cij_kll(1) = C11 ; cij_kll(2) = C12 ; cij_kll(3) = C13
! cij_kll(4) = C14 ; cij_kll(5) = C15 ; cij_kll(6) = C16
! cij_kll(7) = C22 ; cij_kll(8) = C23 ; cij_kll(9) = C24
! cij_kll(10) = C25 ; cij_kll(11) = C26 ; cij_kll(12) = C33
! cij_kll(13) = C34 ; cij_kll(14) = C35 ; cij_kll(15) = C36
! cij_kll(16) = C44 ; cij_kll(17) = C45 ; cij_kll(18) = C46
! cij_kll(19) = C55 ; cij_kll(20) = C56 ; cij_kll(21) = C66
! where the Cij (Voigt's notation) are defined as function of
! the components of the elastic tensor in spherical coordinates
! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)

  use constants, only: CUSTOM_REAL
  

  implicit none
  double precision, parameter :: ONE_HALF = 0.5d0, ONE_FOURTH = 0.25d0, &
                              ONE_EIGHTH = 0.125d0, ONE_SIXTEENTH = 0.0625d0
  real(kind=CUSTOM_REAL), intent(in) :: theta_in,phi_in

  real(kind=CUSTOM_REAL), dimension(21), intent(in) :: cij_kl
  real(kind=CUSTOM_REAL), dimension(21), intent(out) :: cij_kll

  double precision :: theta,phi
  double precision :: costheta,sintheta,cosphi,sinphi
  double precision :: costhetasq,sinthetasq,cosphisq,sinphisq
  double precision :: costwotheta,sintwotheta,costwophi,sintwophi
  double precision :: cosfourtheta,sinfourtheta,cosfourphi,sinfourphi
  double precision :: costhetafour,sinthetafour,cosphifour,sinphifour
  double precision :: sintwophisq,sintwothetasq
  double precision :: costhreetheta,sinthreetheta,costhreephi,sinthreephi

  theta = dble(theta_in)
  phi = dble(phi_in)

  costheta = dcos(theta)
  sintheta = dsin(theta)
  cosphi = dcos(phi)
  sinphi = dsin(phi)

  costhetasq = costheta * costheta
  sinthetasq = sintheta * sintheta
  cosphisq = cosphi * cosphi
  sinphisq = sinphi * sinphi

  costhetafour = costhetasq * costhetasq
  sinthetafour = sinthetasq * sinthetasq
  cosphifour = cosphisq * cosphisq
  sinphifour = sinphisq * sinphisq

  costwotheta = dcos(2.d0*theta)
  sintwotheta = dsin(2.d0*theta)
  costwophi = dcos(2.d0*phi)
  sintwophi = dsin(2.d0*phi)

  costhreetheta=dcos(3.d0*theta)
  sinthreetheta=dsin(3.d0*theta)
  costhreephi=dcos(3.d0*phi)
  sinthreephi=dsin(3.d0*phi)

  cosfourtheta = dcos(4.d0*theta)
  sinfourtheta = dsin(4.d0*theta)
  cosfourphi = dcos(4.d0*phi)
  sinfourphi = dsin(4.d0*phi)
  sintwothetasq = sintwotheta * sintwotheta
  sintwophisq = sintwophi * sintwophi

  cij_kll(1) = ONE_SIXTEENTH * (cij_kl(16) - cij_kl(16)* costwophi + &
     16.d0* cosphi*cosphisq* costhetafour* (cij_kl(1)* cosphi + cij_kl(6)* sinphi) + &
     2.d0* (cij_kl(15) + cij_kl(17))* sintwophi* sintwothetasq - &
     2.d0* (cij_kl(16)* cosfourtheta* sinphisq + &
     2.d0* costhetafour* (-4* cij_kl(7)* sinphifour - &
     (cij_kl(2) + cij_kl(21))* sintwophisq) + &
     8.d0* cij_kl(5)* cosphi*cosphisq* costheta*costhetasq* sintheta - &
     8.d0* cij_kl(8)* costhetasq* sinphisq* sinthetasq - &
     8.d0* cij_kl(12)* sinthetafour + &
     8.d0* cosphisq* costhetasq* sintheta* ((cij_kl(4) + &
     cij_kl(20))* costheta* sinphi - &
     (cij_kl(3) + cij_kl(19))*sintheta) + &
     8.d0* cosphi* costheta* (-cij_kl(11)* costheta*costhetasq* &
     sinphi*sinphisq + (cij_kl(10) + cij_kl(18))* costhetasq* sinphisq* sintheta + &
     cij_kl(14)* sintheta*sinthetasq) + 2.d0* sinphi* (cij_kl(13) + &
     cij_kl(9)* sinphisq)* sintwotheta + &
     sinphi* (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta))

  cij_kll(2) = ONE_FOURTH * (costhetasq* (cij_kl(1) + 3.d0* cij_kl(2) + cij_kl(7) - &
      cij_kl(21) + (-cij_kl(1) + cij_kl(2) - cij_kl(7) + &
      cij_kl(21))* cosfourphi + (-cij_kl(6) + cij_kl(11))* sinfourphi) + &
      4.d0* (cij_kl(8)* cosphisq - cij_kl(15)* cosphi* sinphi + &
      cij_kl(3)* sinphisq)* sinthetasq - &
      2.d0* (cij_kl(10)* cosphisq*cosphi + &
      (cij_kl(9) - cij_kl(20))* cosphisq* sinphi + &
      (cij_kl(5) - cij_kl(18))* cosphi* sinphisq + &
      cij_kl(4)* sinphisq*sinphi)* sintwotheta)

  cij_kll(3) = ONE_EIGHTH * (sintwophi* (3.d0* cij_kl(15) - cij_kl(17) + &
     4.d0* (cij_kl(2) + cij_kl(21))* costhetasq* sintwophi* sinthetasq) + &
     4.d0* cij_kl(12)* sintwothetasq + 4.d0* cij_kl(1)* cosphifour* sintwothetasq + &
     2.d0* cosphi*cosphisq* (8.d0* cij_kl(6)* costhetasq* sinphi* sinthetasq + &
     cij_kl(5)* sinfourtheta) + 2.d0* cosphisq* (3.d0* cij_kl(3) -  cij_kl(19) + &
     (cij_kl(3) + cij_kl(19))* cosfourtheta + &
     (cij_kl(4) + cij_kl(20))* sinphi* sinfourtheta) + &
     2.d0* sinphi* (sinphi* (3.d0* cij_kl(8) - &
     cij_kl(16) + (cij_kl(8) + cij_kl(16))* cosfourtheta + &
     2.d0* cij_kl(7)* sinphisq* sintwothetasq)+ &
     (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta)+ &
     2.d0* cosphi* ((cij_kl(15) + cij_kl(17))* cosfourtheta* sinphi + &
     8.d0* cij_kl(11)* costhetasq* sinphi*sinphisq* sinthetasq + &
     (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq)*sinfourtheta))

  cij_kll(4) = ONE_EIGHTH * (cosphi* costheta *(5.d0* cij_kl(4) - &
     cij_kl(9) + 4.d0* cij_kl(13) - &
     3.d0* cij_kl(20) + (cij_kl(4) + 3.d0* cij_kl(9) - &
     4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + &
     ONE_HALF* (cij_kl(4) - cij_kl(9) + &
     cij_kl(20))* costhreephi * (costheta + 3.d0* costhreetheta) - &
     costheta* (-cij_kl(5) + 5.d0* cij_kl(10) + &
     4.d0* cij_kl(14) - 3.d0* cij_kl(18) + &
     (3.d0* cij_kl(5) + cij_kl(10) - &
     4.d0* cij_kl(14) + cij_kl(18))* costwotheta)* sinphi - &
     ONE_HALF* (cij_kl(5) - cij_kl(10) - cij_kl(18))* (costheta + &
     3.d0* costhreetheta)* sinthreephi + &
     4.d0* (cij_kl(6) - cij_kl(11))* cosfourphi* costhetasq* sintheta - &
     4.d0* (cij_kl(1) + cij_kl(3) - cij_kl(7) - cij_kl(8) + cij_kl(16) - cij_kl(19) + &
     (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + &
     cij_kl(16) - cij_kl(19))* costwotheta)* sintwophi* sintheta - &
     4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
     cij_kl(21))* costhetasq* sinfourphi* sintheta + &
     costwophi* ((cij_kl(6) + cij_kl(11) + 6.d0* cij_kl(15) - &
     2.d0* cij_kl(17))* sintheta + &
     (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + cij_kl(17)))* sinthreetheta))
 
  cij_kll(5) = ONE_FOURTH * (2.d0* (cij_kl(4) + &
     cij_kl(20))* cosphisq* (costwotheta + cosfourtheta)* sinphi + &
     2.d0* cij_kl(9)* (costwotheta + cosfourtheta)* sinphi*sinphisq + &
     16.d0* cij_kl(1)* cosphifour* costheta*costhetasq* sintheta + &
     4.d0* costheta*costhetasq* (-2.d0* cij_kl(8)* sinphisq + &
     4.d0* cij_kl(7)* sinphifour + &
     (cij_kl(2) + cij_kl(21))* sintwophisq)* sintheta + &
     4.d0* cij_kl(13)* (1.d0 + 2.d0* costwotheta)* sinphi* sinthetasq + &
     8.d0* costheta* (-2.d0* cij_kl(12) + cij_kl(8)* sinphisq)* sintheta*sinthetasq + &
     2.d0* cosphi*cosphisq* (cij_kl(5)* (costwotheta + cosfourtheta) + &
     8.d0* cij_kl(6)* costheta*costhetasq* sinphi* sintheta) + &
     2.d0* cosphi* (cosfourtheta* (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
     costwotheta* (cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
     8.d0* cij_kl(11)* costheta*costhetasq* sinphi*sinphisq* sintheta) - &
     (cij_kl(3) + cij_kl(16) + cij_kl(19) + &
     (cij_kl(3) - cij_kl(16) + cij_kl(19))* costwophi + &
     (cij_kl(15) + cij_kl(17))* sintwophi)* sinfourtheta)

  cij_kll(6) = ONE_HALF * costheta*costhetasq* ((cij_kl(6) + cij_kl(11))* costwophi + &
      (cij_kl(6) - cij_kl(11))* cosfourphi + 2.d0* (-cij_kl(1) + cij_kl(7))* sintwophi + &
      (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* sinfourphi) + &
      ONE_FOURTH* costhetasq* (-(cij_kl(4) + 3* cij_kl(9) + cij_kl(20))* cosphi - &
      3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi + &
      (3.d0* cij_kl(5) + cij_kl(10) + cij_kl(18))* sinphi + &
      3.d0* (cij_kl(5) - cij_kl(10) - cij_kl(18))* sinthreephi)* sintheta + &
      costheta* ((cij_kl(15) + cij_kl(17))* costwophi + &
      (-cij_kl(3) + cij_kl(8) + cij_kl(16) - cij_kl(19))* sintwophi)* sinthetasq + &
      (-cij_kl(13)* cosphi + cij_kl(14)* sinphi)* sintheta*sinthetasq

  cij_kll(7) = cij_kl(7) * cosphifour - cij_kl(11)* cosphi*cosphisq* sinphi + &
      (cij_kl(2) + cij_kl(21))* cosphisq* sinphisq - &
      cij_kl(6)* cosphi* sinphi*sinphisq + &
      cij_kl(1)* sinphifour

  cij_kll(8) = ONE_HALF * (2.d0* costhetasq* sinphi* (-cij_kl(15)* cosphi + &
      cij_kl(3)* sinphi) + 2.d0* cij_kl(2)* cosphifour* sinthetasq + &
      (2.d0* cij_kl(2)* sinphifour + &
      (cij_kl(1) + cij_kl(7) - cij_kl(21))* sintwophisq)* sinthetasq + &
      cij_kl(4)* sinphi*sinphisq* sintwotheta + &
      cosphi*cosphisq* (2.d0* (-cij_kl(6) + cij_kl(11))* sinphi* sinthetasq + &
      cij_kl(10)* sintwotheta) + cosphi* sinphisq* (2.d0* (cij_kl(6) - &
      cij_kl(11))* sinphi* sinthetasq + &
      (cij_kl(5) - cij_kl(18))* sintwotheta) + &
      cosphisq* (2.d0* cij_kl(8)* costhetasq + &
      (cij_kl(9) - cij_kl(20))* sinphi* sintwotheta))

  cij_kll(9) = cij_kl(11) * cosphifour* sintheta - sinphi*sinphisq* (cij_kl(5)* costheta + &
      cij_kl(6)* sinphi* sintheta) +  cosphisq* sinphi* (-(cij_kl(10) + &
      cij_kl(18))* costheta + &
      3.d0* (cij_kl(6) - cij_kl(11))* sinphi* sintheta) + &
      cosphi* sinphisq* ((cij_kl(4) + cij_kl(20))* costheta + &
      2.d0* (-2.d0* cij_kl(1) + cij_kl(2) + cij_kl(21))* sinphi* sintheta) + &
      cosphi*cosphisq* (cij_kl(9)* costheta - 2.d0* (cij_kl(2) - 2.d0* cij_kl(7) + &
      cij_kl(21))* sinphi* sintheta)

  cij_kll(10) = ONE_FOURTH * (4.d0* costwotheta* (cij_kl(10)* cosphi*cosphisq +&
      (cij_kl(9) - cij_kl(20))* cosphisq* sinphi + &
      (cij_kl(5) - cij_kl(18))* cosphi* sinphisq + &
      cij_kl(4)* sinphi*sinphisq) + (cij_kl(1) + 3.d0* cij_kl(2) - &
      2.d0* cij_kl(3) + cij_kl(7) - &
      2.d0* cij_kl(8) - cij_kl(21) + 2.d0* (cij_kl(3) - cij_kl(8))* costwophi + &
      (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
      2.d0* cij_kl(15)* sintwophi + &
      (-cij_kl(6) + cij_kl(11))* sinfourphi)* sintwotheta)

  cij_kll(11) = ONE_FOURTH * (2.d0* costheta* ((cij_kl(6) + cij_kl(11))* costwophi + &
      (-cij_kl(6) + cij_kl(11))* cosfourphi + &
      2.d0* (-cij_kl(1) + cij_kl(7))* sintwophi + &
      (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(21))* sinfourphi) + &
      (-(cij_kl(4) + 3.d0* cij_kl(9) + cij_kl(20))* cosphi + &
      (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi + &
      (3.d0* cij_kl(5) + cij_kl(10) + cij_kl(18))* sinphi + &
      (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sintheta)

  cij_kll(12) = ONE_SIXTEENTH * (cij_kl(16) - 2.d0* cij_kl(16)* cosfourtheta* sinphisq + &
      costwophi* (-cij_kl(16) + 8.d0* costheta* sinthetasq* ((cij_kl(3) - &
      cij_kl(8) + cij_kl(19))* costheta + &
      (cij_kl(5) - cij_kl(10) - cij_kl(18))* cosphi* sintheta)) + &
      2.d0* (cij_kl(15) + cij_kl(17))* sintwophi* sintwothetasq + &
      2.d0* (8.d0* cij_kl(12)* costhetafour + &
      8.d0* cij_kl(14)* cosphi* costheta*costhetasq* sintheta + &
      4.d0* cosphi* costheta* (cij_kl(5) + cij_kl(10) + cij_kl(18) + &
      (cij_kl(4) + cij_kl(20))* sintwophi)* &
      sintheta*sinthetasq + 8.d0* cij_kl(1)* cosphifour* sinthetafour + &
      8.d0* cij_kl(6)* cosphi*cosphisq* sinphi* sinthetafour + &
      8.d0* cij_kl(11)* cosphi* sinphi*sinphisq* sinthetafour + &
      8.d0* cij_kl(7)* sinphifour* sinthetafour + &
      2.d0* cij_kl(2)* sintwophisq* sinthetafour + &
      2.d0* cij_kl(21)* sintwophisq* sinthetafour + &
      2.d0* cij_kl(13)* sinphi* sintwotheta + &
      2.d0* cij_kl(9)* sinphi*sinphisq* sintwotheta + &
      cij_kl(3)* sintwothetasq + cij_kl(8)* sintwothetasq + &
      cij_kl(19)* sintwothetasq + cij_kl(13)* sinphi* sinfourtheta - &
      cij_kl(9)* sinphi*sinphisq* sinfourtheta))

  cij_kll(13) = ONE_EIGHTH * (cosphi* costheta* (cij_kl(4) + 3.d0* cij_kl(9) + &
      4.d0* cij_kl(13) + cij_kl(20) - (cij_kl(4) + 3.d0* cij_kl(9) - &
      4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + 4.d0* (-cij_kl(1) - &
      cij_kl(3) + cij_kl(7) + cij_kl(8) + cij_kl(16) - cij_kl(19) + &
      (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
      cij_kl(19))* costwotheta)* sintwophi* sintheta + &
      4.d0* (cij_kl(6) - cij_kl(11))* cosfourphi* sinthetasq*sintheta - &
      4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* sinfourphi* sinthetasq*sintheta + &
      costheta* ((-3.d0* cij_kl(5) - cij_kl(10) - 4.d0* cij_kl(14) - &
      cij_kl(18) + (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + &
      cij_kl(18))* costwotheta)* sinphi + 6.d0* ((cij_kl(4) - cij_kl(9) + &
      cij_kl(20))* costhreephi + (-cij_kl(5) + cij_kl(10) + &
      cij_kl(18))* sinthreephi)* sinthetasq) + costwophi* ((3* cij_kl(6) + &
      3.d0* cij_kl(11) + 2.d0* (cij_kl(15) + cij_kl(17)))* sintheta - &
      (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + &
      cij_kl(17)))* sinthreetheta))

  cij_kll(14) = ONE_FOURTH * (2.d0* cij_kl(13)* (costwotheta + cosfourtheta)* sinphi + &
      8.d0* costheta*costhetasq* (-2.d0* cij_kl(12) + cij_kl(8)* sinphisq)* sintheta + &
      4.d0* (cij_kl(4) + cij_kl(20))* cosphisq* (1.d0 + &
      2.d0* costwotheta)* sinphi* sinthetasq + &
      4.d0* cij_kl(9)* (1.d0 + 2.d0* costwotheta)* sinphi*sinphisq* sinthetasq + &
      16.d0* cij_kl(1)* cosphifour* costheta* sintheta*sinthetasq + &
      4.d0* costheta* (-2.d0* cij_kl(8)* sinphisq + 4.d0* cij_kl(7)* sinphifour + &
      (cij_kl(2) + cij_kl(21))* sintwophisq)* sintheta*sinthetasq + &
      4.d0* cosphi*cosphisq* sinthetasq* (cij_kl(5) + 2.d0* cij_kl(5)* costwotheta + &
      4.d0* cij_kl(6)* costheta* sinphi* sintheta) + &
      2.d0* cosphi* (cosfourtheta* (cij_kl(14) - (cij_kl(10) + cij_kl(18))* sinphisq) + &
      costwotheta* (cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
      8.d0* cij_kl(11)* costheta* sinphi*sinphisq* sintheta*sinthetasq) + &
      (cij_kl(3) + cij_kl(16) + cij_kl(19) + (cij_kl(3) - cij_kl(16) + &
      cij_kl(19))* costwophi + (cij_kl(15) + cij_kl(17))* sintwophi)* sinfourtheta)

  cij_kll(15) = costwophi * costheta* (-cij_kl(17) + (cij_kl(15) + cij_kl(17))* costhetasq) + &
       ONE_SIXTEENTH* (-((11.d0* cij_kl(4) + cij_kl(9) + 4.d0* cij_kl(13) - &
       5.d0* cij_kl(20))* cosphi + (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (cij_kl(5) + 11.d0* cij_kl(10) + 4.d0* cij_kl(14) - &
       5.d0* cij_kl(18))* sinphi + (-cij_kl(5) + cij_kl(10) + &
       cij_kl(18))* sinthreephi)* sintheta + &
       8.d0* costheta* ((-cij_kl(1) - cij_kl(3) + cij_kl(7) + cij_kl(8) - cij_kl(16) +&
       cij_kl(19) + (cij_kl(1) - cij_kl(3) - &
       cij_kl(7) + cij_kl(8) + cij_kl(16) - cij_kl(19))* costwotheta)* sintwophi +&
       ((cij_kl(6) + cij_kl(11))* costwophi + &
       (cij_kl(6) - cij_kl(11))* cosfourphi + (-cij_kl(1) + cij_kl(2) - cij_kl(7) +&
       cij_kl(21))* sinfourphi)* sinthetasq) +&
       ((cij_kl(4) + 3.d0* cij_kl(9) - 4.d0* cij_kl(13) + cij_kl(20))* cosphi + &
       3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))* sinphi + &
       3.d0* (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sinthreetheta)

  cij_kll(16) = ONE_FOURTH *(cij_kl(1) - cij_kl(2) + cij_kl(7) + cij_kl(16) + &
       cij_kl(19) + cij_kl(21) + 2.d0*(cij_kl(16) - cij_kl(19))*costwophi* costhetasq + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(16) + &
       cij_kl(19) - cij_kl(21))*costwotheta - 2.d0* cij_kl(17)* costhetasq* sintwophi + &
       2.d0* ((-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
       (-cij_kl(6) + cij_kl(11))* sinfourphi)* sinthetasq + ((cij_kl(5) - cij_kl(10) +&
       cij_kl(18))* cosphi + (-cij_kl(5) + cij_kl(10) + cij_kl(18))* costhreephi +&
       (-cij_kl(4) + cij_kl(9) + cij_kl(20))* sinphi - &
       (cij_kl(4) - cij_kl(9) + cij_kl(20))* sinthreephi)* sintwotheta)

  cij_kll(17) = ONE_EIGHTH * (4.d0* costwophi* costheta* (cij_kl(6) + cij_kl(11) - &
       2.d0* cij_kl(15) - (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + &
       cij_kl(17)))* costwotheta) - (2.d0* cosphi* (-3.d0* cij_kl(4) +&
       cij_kl(9) + 2.d0* cij_kl(13) + cij_kl(20) + (cij_kl(4) - cij_kl(9) + &
       cij_kl(20))* costwophi) - (cij_kl(5) - 5.d0* cij_kl(10) + &
       4.d0* cij_kl(14) + 3.d0* cij_kl(18))* sinphi + (-cij_kl(5) + cij_kl(10) + &
       cij_kl(18))* sinthreephi)* sintheta + &
       8.d0* costheta* ((-cij_kl(1) + cij_kl(3) + cij_kl(7) - cij_kl(8) + &
       (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
       cij_kl(19))* costwotheta)* sintwophi + ((cij_kl(6) - cij_kl(11))* cosfourphi + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* sinfourphi)* sinthetasq) +&
       ((cij_kl(4) + 3.d0* cij_kl(9) - 4.d0* cij_kl(13) + cij_kl(20))* cosphi + &
       3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))* sinphi + &
       3.d0* (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sinthreetheta)

  cij_kll(18) = ONE_HALF * ((cij_kl(5) - cij_kl(10) + cij_kl(18))* cosphi* costwotheta - &
       (cij_kl(5) - cij_kl(10) - cij_kl(18))* costhreephi* costwotheta - &
       2.d0* (cij_kl(4) - cij_kl(9) + &
       (cij_kl(4) - cij_kl(9) + cij_kl(20))* costwophi)* costwotheta* sinphi + &
       (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(16) - cij_kl(19) + cij_kl(21) + &
       (-cij_kl(16) + cij_kl(19))* costwophi + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
       cij_kl(17)* sintwophi + &
       (-cij_kl(6) + cij_kl(11))* sinfourphi)* sintwotheta)

  cij_kll(19) = ONE_FOURTH * (cij_kl(16) - cij_kl(16)* costwophi + &
      (-cij_kl(15) + cij_kl(17))* sintwophi + &
      4.d0* cij_kl(12)* sintwothetasq + &
      2.d0* (2.d0* cij_kl(1)* cosphifour* sintwothetasq + &
      cosphi*cosphisq* (8.d0* cij_kl(6)* costhetasq* sinphi* sinthetasq + &
      cij_kl(5)* sinfourtheta) + cosphisq* (-cij_kl(3) + cij_kl(19) + (cij_kl(3) +&
      cij_kl(19))* cosfourtheta + (cij_kl(4) + cij_kl(20))* sinphi* sinfourtheta) + &
      sinphi* (cosfourtheta* ((cij_kl(15) + cij_kl(17))* cosphi + &
      cij_kl(16)* sinphi) + (cij_kl(2) + cij_kl(7) - 2.d0* cij_kl(8) + cij_kl(21) + &
      (cij_kl(2) - cij_kl(7) + cij_kl(21))* costwophi)* sinphi* sintwothetasq + &
      (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta) + &
      cosphi* (8.d0* cij_kl(11)* costhetasq* sinphi*sinphisq* sinthetasq + &
      (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq)* sinfourtheta)))

  cij_kll(20) = ONE_EIGHTH * (2.d0* cosphi* costheta* (-3.d0* cij_kl(4) - cij_kl(9) + &
      4.d0* cij_kl(13) + cij_kl(20) + (cij_kl(4) + 3.d0* cij_kl(9) - &
      4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + &
      (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi* (costheta + &
      3.d0* costhreetheta) - &
      2.d0* costheta* (-cij_kl(5) - 3.d0* cij_kl(10) + 4.d0* cij_kl(14) + &
      cij_kl(18) + (3.d0* cij_kl(5) + &
      cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))*costwotheta)* sinphi - &
      (cij_kl(5) - cij_kl(10) - cij_kl(18))* &
      (costheta + 3.d0* costhreetheta)* sinthreephi + 8.d0* (cij_kl(6) - &
      cij_kl(11))* cosfourphi* costhetasq* sintheta - 8.d0* (cij_kl(1) - &
      cij_kl(3) - cij_kl(7) + cij_kl(8) + &
      (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
      cij_kl(19))* costwotheta)* sintwophi* sintheta - &
      8.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* costhetasq* sinfourphi* sintheta + &
      2.d0* costwophi* ((cij_kl(6) + cij_kl(11) - 2.d0* cij_kl(15) + &
      2.d0* cij_kl(17))* sintheta + &
      (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + cij_kl(17)))* sinthreetheta))

  cij_kll(21) = ONE_FOURTH * (cij_kl(1) - cij_kl(2) + cij_kl(7) + cij_kl(16) + &
      cij_kl(19) + cij_kl(21) - 2.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* cosfourphi* costhetasq + &
      (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(16) - cij_kl(19) + &
      cij_kl(21))* costwotheta + &
      2.d0* (-cij_kl(6) + cij_kl(11))* costhetasq* sinfourphi - &
      2.d0* ((-cij_kl(16) + cij_kl(19))* costwophi + cij_kl(17)* sintwophi)* sinthetasq - &
      ((cij_kl(5) - cij_kl(10) + cij_kl(18))* cosphi + (-cij_kl(5) + cij_kl(10) +&
      cij_kl(18))* costhreephi + &
      (-cij_kl(4) + cij_kl(9) + cij_kl(20))* sinphi - (cij_kl(4) - cij_kl(9) + &
      cij_kl(20))* sinthreephi)* sintwotheta)

  end subroutine rotate_kernels_dble
