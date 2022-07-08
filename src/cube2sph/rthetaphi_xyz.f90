!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

  subroutine xyz_2_rthetaphi(x,y,z,r,theta,phi)

! convert x y z to r theta phi, single precision call

  use constants, only: CUSTOM_REAL,ZERO

  implicit none
  double precision, parameter :: SMALL_VAL_ANGLE = 1.d-10

  real(kind=CUSTOM_REAL), intent(in) :: x,y,z
  real(kind=CUSTOM_REAL), intent(out) :: r,theta,phi

  double precision xmesh,ymesh,zmesh

  xmesh = dble(x)
  ymesh = dble(y)
  zmesh = dble(z)

  if (zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
  if (zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE

  theta = real(datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh), kind=CUSTOM_REAL)

  if (xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
  if (xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE

  phi = real(datan2(ymesh,xmesh), kind=CUSTOM_REAL)

  r = real(dsqrt(xmesh*xmesh + ymesh*ymesh + zmesh*zmesh), kind=CUSTOM_REAL)

  end subroutine xyz_2_rthetaphi


