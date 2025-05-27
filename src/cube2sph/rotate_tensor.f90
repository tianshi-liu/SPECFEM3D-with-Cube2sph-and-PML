subroutine rotate_tensor(cij_in,rotmat,cij_out)

  ! rotates from (6,6)-tensor cij_in to cij_out using the (3,3)-rotation matrix rotmat
  !
  ! cij (Voigt) tensors only need upper triangle set

  implicit none

  double precision,dimension(3,3),intent(in) :: rotmat
  double precision,dimension(6,6),intent(in) :: cij_in
  double precision,dimension(6,6),intent(out) :: cij_out

  ! local parameters
  double precision,dimension(6,6) :: bond,bond_t
  double precision,dimension(6,6) :: tensor_tmp
  integer :: i,j,k

  ! creates Bond matrix (see e.g. Auld, 1973)
  ! First column
  bond(1,1) = rotmat(1,1)*rotmat(1,1)
  bond(2,1) = rotmat(2,1)*rotmat(2,1)
  bond(3,1) = rotmat(3,1)*rotmat(3,1)
  bond(4,1) = rotmat(2,1)*rotmat(3,1)
  bond(5,1) = rotmat(3,1)*rotmat(1,1)
  bond(6,1) = rotmat(1,1)*rotmat(2,1)

  ! Second column
  bond(1,2) = rotmat(1,2)*rotmat(1,2)
  bond(2,2) = rotmat(2,2)*rotmat(2,2)
  bond(3,2) = rotmat(3,2)*rotmat(3,2)
  bond(4,2) = rotmat(2,2)*rotmat(3,2)
  bond(5,2) = rotmat(3,2)*rotmat(1,2)
  bond(6,2) = rotmat(1,2)*rotmat(2,2)

  ! Third column
  bond(1,3) = rotmat(1,3)*rotmat(1,3)
  bond(2,3) = rotmat(2,3)*rotmat(2,3)
  bond(3,3) = rotmat(3,3)*rotmat(3,3)
  bond(4,3) = rotmat(2,3)*rotmat(3,3)
  bond(5,3) = rotmat(3,3)*rotmat(1,3)
  bond(6,3) = rotmat(1,3)*rotmat(2,3)

  ! Fourth column
  bond(1,4) = 2.d0*rotmat(1,2)*rotmat(1,3)
  bond(2,4) = 2.d0*rotmat(2,2)*rotmat(2,3)
  bond(3,4) = 2.d0*rotmat(3,2)*rotmat(3,3)
  bond(4,4) = rotmat(2,2)*rotmat(3,3) + rotmat(2,3)*rotmat(3,2)
  bond(5,4) = rotmat(1,2)*rotmat(3,3) + rotmat(1,3)*rotmat(3,2)
  bond(6,4) = rotmat(1,2)*rotmat(2,3) + rotmat(1,3)*rotmat(2,2)

  ! Fifth column
  bond(1,5) = 2.d0*rotmat(1,3)*rotmat(1,1)
  bond(2,5) = 2.d0*rotmat(2,3)*rotmat(2,1)
  bond(3,5) = 2.d0*rotmat(3,3)*rotmat(3,1)
  bond(4,5) = rotmat(2,1)*rotmat(3,3) + rotmat(2,3)*rotmat(3,1)
  bond(5,5) = rotmat(1,3)*rotmat(3,1) + rotmat(1,1)*rotmat(3,3)
  bond(6,5) = rotmat(1,3)*rotmat(2,1) + rotmat(1,1)*rotmat(2,3)

  ! Sixth column
  bond(1,6) = 2.d0*rotmat(1,1)*rotmat(1,2)
  bond(2,6) = 2.d0*rotmat(2,1)*rotmat(2,2)
  bond(3,6) = 2.d0*rotmat(3,1)*rotmat(3,2)
  bond(4,6) = rotmat(2,2)*rotmat(3,1) + rotmat(2,1)*rotmat(3,2)
  bond(5,6) = rotmat(1,1)*rotmat(3,2) + rotmat(1,2)*rotmat(3,1)
  bond(6,6) = rotmat(1,1)*rotmat(2,2) + rotmat(1,2)*rotmat(2,1)

  bond_t = transpose(bond)

  ! rotates Cij
  ! First compute C M^t
  tensor_tmp(:,:) = 0.d0
  do j = 1,6
    do k = 1,6
      do i = 1,6
        tensor_tmp(i,j) = tensor_tmp(i,j) + cij_in(i,k) * bond_t(k,j)
      enddo
    enddo
  enddo
  ! Second compute M * (C M^t)
  cij_out(:,:) = 0.d0
  do j = 1,6
    do k = 1,6
      do i = 1,j ! half only
        cij_out(i,j) = cij_out(i,j) + bond(i,k) * tensor_tmp(k,j)
      enddo
    enddo
  enddo

  end subroutine rotate_tensor

subroutine rotate_tensor_global_to_radial(theta,phi, &
  d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
  d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
  c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  ! rotates from global (c_ij) to local (d_ij) anisotropic parameters.
  ! The c_ij are the coefficients in the global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(out) :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
  d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  double precision,intent(in) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! Purpose : compute the kernels in r,theta,phi (cij_kl_spherical)
  ! from the kernels in x,y,z (cij_kl) (x,y,z to r,theta,phi)
  ! At r,theta,phi fixed
  ! theta and phi are in radians

  ! local parameters
  double precision,dimension(3,3) :: rotmat
  double precision,dimension(6,6) :: cij,dij
  integer :: i,j

  ! rotation matrix
  ! rotates pole (Cartesian) to spherical (radial) position
  ! First column
  rotmat(1,1) =  cos(phi) * cos(theta)
  rotmat(2,1) = -sin(phi)
  rotmat(3,1) =  cos(phi) * sin(theta)

  ! Second column
  rotmat(1,2) =  sin(phi) * cos(theta)
  rotmat(2,2) =  cos(phi)
  rotmat(3,2) =  sin(phi) * sin(theta)

  ! Third column
  rotmat(1,3) = -sin(theta)
  rotmat(2,3) =  0.d0
  rotmat(3,3) =  cos(theta)

  ! Cij Voigt notation
  cij(1,1) = c11
  cij(1,2) = c12
  cij(1,3) = c13
  cij(1,4) = c14
  cij(1,5) = c15
  cij(1,6) = c16
  cij(2,2) = c22
  cij(2,3) = c23
  cij(2,4) = c24
  cij(2,5) = c25
  cij(2,6) = c26
  cij(3,3) = c33
  cij(3,4) = c34
  cij(3,5) = c35
  cij(3,6) = c36
  cij(4,4) = c44
  cij(4,5) = c45
  cij(4,6) = c46
  cij(5,5) = c55
  cij(5,6) = c56
  cij(6,6) = c66

  ! fills lower-triangle, for example: C(2,1) = C(1,2) <->  c21 = c12
  !                                    C(3,1) = C(1,3) <->  c31 = c13
  !                                    C(3,2) = C(2,3) <->  c32 = c23
    do j = 2,6
      do i = 1,j - 1
        cij(j,i) = cij(i,j)
      enddo
    enddo

  call rotate_tensor(cij,rotmat,dij)

  ! returns local dij
  d11 = dij(1,1)
  d12 = dij(1,2)
  d13 = dij(1,3)
  d14 = dij(1,4)
  d15 = dij(1,5)
  d16 = dij(1,6)
  d22 = dij(2,2)
  d23 = dij(2,3)
  d24 = dij(2,4)
  d25 = dij(2,5)
  d26 = dij(2,6)
  d33 = dij(3,3)
  d34 = dij(3,4)
  d35 = dij(3,5)
  d36 = dij(3,6)
  d44 = dij(4,4)
  d45 = dij(4,5)
  d46 = dij(4,6)
  d55 = dij(5,5)
  d56 = dij(5,6)
  d66 = dij(6,6)


end subroutine rotate_tensor_global_to_radial

subroutine rotate_tensor_radial_to_global(theta,phi, &
  d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
  d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
  c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  ! rotates from local (d_ij) to global (c_ij) anisotropic parameters.
  ! The c_ij are the coefficients in the global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(in) :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
  d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision,dimension(3,3) :: rotmat
  double precision,dimension(6,6) :: cij,dij
  integer :: i,j
  ! original
  ! double precision :: costheta,sintheta,cosphi,sinphi
  ! double precision :: costhetasq,sinthetasq,cosphisq,sinphisq
  ! double precision :: costwotheta,sintwotheta,costwophi,sintwophi
  ! double precision :: cosfourtheta,sinfourtheta
  ! double precision :: costhetafour,sinthetafour,cosphifour,sinphifour
  ! double precision :: sintwophisq,sintwothetasq

  ! by default, we use Bond's matrix approach which is also used for the backward rotation
  !logical, parameter :: USE_BOND_MATRIX_ROTATION = .true.

  ! tensor rotation
  ! rotation matrix
  ! rotates from local (radial) to pole (Cartesian) reference
  ! First column
  rotmat(1,1) =  cos(phi) * cos(theta)
  rotmat(2,1) =  sin(phi) * cos(theta)
  rotmat(3,1) = -sin(theta)

  ! Second column
  rotmat(1,2) = -sin(phi)
  rotmat(2,2) =  cos(phi)
  rotmat(3,2) =  0.d0

  ! Third column
  rotmat(1,3) =  cos(phi)*sin(theta)
  rotmat(2,3) =  sin(phi)*sin(theta)
  rotmat(3,3) =  cos(theta)

  ! Cij Voigt notation
  dij(1,1) = d11
  dij(1,2) = d12
  dij(1,3) = d13
  dij(1,4) = d14
  dij(1,5) = d15
  dij(1,6) = d16
  dij(2,2) = d22
  dij(2,3) = d23
  dij(2,4) = d24
  dij(2,5) = d25
  dij(2,6) = d26
  dij(3,3) = d33
  dij(3,4) = d34
  dij(3,5) = d35
  dij(3,6) = d36
  dij(4,4) = d44
  dij(4,5) = d45
  dij(4,6) = d46
  dij(5,5) = d55
  dij(5,6) = d56
  dij(6,6) = d66
  ! fills lower-triangle, for example: C(2,1) = C(1,2) <->  c21 = c12
  !                                    C(3,1) = C(1,3) <->  c31 = c13
  !                                    C(3,2) = C(2,3) <->  c32 = c23
  do j = 2,6
    do i = 1,j - 1
      dij(j,i) = dij(i,j)
    enddo
  enddo

  call rotate_tensor(dij,rotmat,cij)

  ! returns global cij (SPECFEM x/y/z reference)
  c11 = cij(1,1)
  c12 = cij(1,2)
  c13 = cij(1,3)
  c14 = cij(1,4)
  c15 = cij(1,5)
  c16 = cij(1,6)
  c22 = cij(2,2)
  c23 = cij(2,3)
  c24 = cij(2,4)
  c25 = cij(2,5)
  c26 = cij(2,6)
  c33 = cij(3,3)
  c34 = cij(3,4)
  c35 = cij(3,5)
  c36 = cij(3,6)
  c44 = cij(4,4)
  c45 = cij(4,5)
  c46 = cij(4,6)
  c55 = cij(5,5)
  c56 = cij(5,6)
  c66 = cij(6,6)

end subroutine rotate_tensor_radial_to_global