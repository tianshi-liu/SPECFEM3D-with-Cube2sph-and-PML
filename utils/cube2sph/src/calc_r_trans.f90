subroutine calc_r_trans_cube2sph(xstore,ystore,zstore,nspec,nspec_cpml,&
              r_trans,r_trans_inv,r_earth,center_lat,center_lon,rotation_azi)
  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM, RADIANS_TO_DEGREES, PI
  use generate_databases_par, only: CPML_to_spec
  integer :: nspec,nspec_cpml,ispec,ispec_cpml,i,j,k
  double precision :: r_earth,center_lat,center_lon,rotation_azi
  double precision :: &
          alpha, beta, gamma, cosa, sina, cosb, sinb, cosg, sing, t1, t2, &
          xi, eta, x, y, z, c1,c2
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
  real, dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,nspec_cpml) :: r_trans,&
          r_trans_inv
  double precision, dimension(NDIM,NDIM) :: J_trans, J_trans_inv, rot_mat,temp
  double precision :: l1,l2,l3
  double precision :: tempd,tempdsqrt,tempz
  alpha = center_lon / RADIANS_TO_DEGREES
  beta = PI / 2.0 - center_lat / RADIANS_TO_DEGREES
  gamma = rotation_azi / RADIANS_TO_DEGREES
  cosa = cos(alpha)
  sina = sin(alpha)
  cosb = cos(beta)
  sinb = sin(beta)
  cosg = cos(gamma)
  sing = sin(gamma)
  rot_mat(1,1) = cosg * cosb * cosa - sing * sina
  rot_mat(2,1) = - sing* cosb * cosa - cosg * sina
  rot_mat(3,1) = sinb * cosa
  rot_mat(1,2) = cosg * cosb * sina + sing * cosa
  rot_mat(2,2) = - sing * cosb * sina + cosg * cosa
  rot_mat(3,2) = sinb * sina
  rot_mat(1,3) = - cosg * sinb
  rot_mat(2,3) = sing * sinb
  rot_mat(3,3) = cosb
  do ispec_cpml = 1, nspec_cpml
    ispec = CPML_to_spec(ispec_cpml)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          x = xstore(i,j,k,ispec)
          y = ystore(i,j,k,ispec)
          z = zstore(i,j,k,ispec)
          xi = x / r_earth
          eta = y / r_earth
          t1 = tan(xi)
          t2 = tan(eta)
          c1 = cos(xi)
          c2 = cos(eta)
          c1 = c1 * c1
          c2 = c2 * c2
          tempd = 1.0 + t1*t1+t2*t2
          tempdsqrt = sqrt(tempd)
          tempz = 1.0 + z / r_earth
          J_trans_inv(1,1)=tempz*t1*t2/c1/tempd/tempdsqrt
          J_trans_inv(2,1)=-tempz/tempd/c1/c2/tempdsqrt
          J_trans_inv(3,1)=-t2/tempdsqrt
          J_trans_inv(1,2)=tempz/tempd/c1/c2/tempdsqrt
          J_trans_inv(2,2)=-tempz*t1*t2/c2/tempd/tempdsqrt
          J_trans_inv(3,2)=t1/tempdsqrt
          J_trans_inv(1,3)=-tempz*t1/tempd/c1/tempdsqrt
          J_trans_inv(2,3)=-tempz*t2/tempd/c2/tempdsqrt
          J_trans_inv(3,3)=1.0/tempdsqrt
          call calc_mat_mul(J_trans_inv,rot_mat,temp,NDIM)
          J_trans_inv(:,:)=temp(:,:)
          !call calc_mat_inv3(J_trans_inv,J_trans,NDIM)
          l1 = dsqrt(J_trans_inv(1,1)**2+J_trans_inv(1,2)**2+&
                  J_trans_inv(1,3)**2)
          l2 = dsqrt(J_trans_inv(2,1)**2+J_trans_inv(2,2)**2+&
                  J_trans_inv(2,3)**2)
          l3 = dsqrt(J_trans_inv(3,1)**2+J_trans_inv(3,2)**2+&
                  J_trans_inv(3,3)**2)
          r_trans_inv(1,:,i,j,k,ispec_cpml) = real(J_trans_inv(1,:) / l1)
          r_trans_inv(2,:,i,j,k,ispec_cpml) = real(J_trans_inv(2,:) / l2)
          r_trans_inv(3,:,i,j,k,ispec_cpml) = real(J_trans_inv(3,:) / l3)
          call calc_mat_inv3(dble(r_trans_inv(:,:,i,j,k,ispec_cpml)),temp,NDIM)
          r_trans(:,1,i,j,k,ispec_cpml) = real(temp(1,:))
          r_trans(:,2,i,j,k,ispec_cpml) = real(temp(2,:))
          r_trans(:,3,i,j,k,ispec_cpml) = real(temp(3,:))
        enddo
      enddo
    enddo
  enddo

end subroutine calc_r_trans_cube2sph



subroutine calc_r_trans(nodes_coords,nodes_coords_new,npts,nspec_cpml,&
                        r_trans,r_trans_inv)
  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,NDIM,GAUSSALPHA,GAUSSBETA
  use generate_databases_par, only: elmnts_ext_mesh,CPML_to_spec
  implicit none
  integer :: npts,nspec_cpml,ispec,ispec_cpml,i,j,k,ia
  double precision, dimension(NDIM,npts) :: nodes_coords,nodes_coords_new
  real, dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,nspec_cpml) :: r_trans,&
          r_trans_inv
  double precision :: shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision, dimension(NGNOD) :: xelm,yelm,zelm,&
                                        xelm_new,yelm_new,zelm_new
  double precision, dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ) :: J_local,&
          J_local_inv, J_local_new,J_local_inv_new
  double precision, dimension(NDIM,NDIM) :: J_trans, J_trans_inv, temp
  !double precision, dimension(NDIM,NDIM) :: temp_check
  double precision :: l1,l2,l3
  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! get the 3-D shape functions
  call get_shape3D(shape3D,dershape3D,xigll,yigll,zigll)

  do ispec_cpml = 1, nspec_cpml
    ispec = CPML_to_spec(ispec_cpml)
    do ia = 1, NGNOD
      xelm(ia) = nodes_coords(1,elmnts_ext_mesh(ia,ispec))
      yelm(ia) = nodes_coords(2,elmnts_ext_mesh(ia,ispec))
      zelm(ia) = nodes_coords(3,elmnts_ext_mesh(ia,ispec))
      xelm_new(ia) = nodes_coords_new(1,elmnts_ext_mesh(ia,ispec))
      yelm_new(ia) = nodes_coords_new(2,elmnts_ext_mesh(ia,ispec))
      zelm_new(ia) = nodes_coords_new(3,elmnts_ext_mesh(ia,ispec))
    enddo
    call get_J_local(xelm,yelm,zelm,dershape3D,J_local,J_local_inv)
    call get_J_local(xelm_new,yelm_new,zelm_new,dershape3D,&
            J_local_new,J_local_inv_new)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          call calc_mat_mul(J_local(:,:,i,j,k),J_local_inv_new(:,:,i,j,k),&
                  J_trans,NDIM)
          call calc_mat_mul(J_local_new(:,:,i,j,k), J_local_inv(:,:,i,j,k),&
                  J_trans_inv,NDIM)
          !call calc_mat_mul(J_trans,J_trans_inv,temp_check,NDIM)
          l1 = dsqrt(J_trans_inv(1,1)**2+J_trans_inv(2,1)**2+&
                  J_trans_inv(3,1)**2)
          l2 = dsqrt(J_trans_inv(1,2)**2+J_trans_inv(2,2)**2+&
                  J_trans_inv(3,2)**2)
          l3 = dsqrt(J_trans_inv(1,3)**2+J_trans_inv(2,3)**2+&
                  J_trans_inv(3,3)**2)
          r_trans(1,:,i,j,k,ispec_cpml) = real(J_trans(1,:) * l1)
          r_trans(2,:,i,j,k,ispec_cpml) = real(J_trans(2,:) * l2)
          r_trans(3,:,i,j,k,ispec_cpml) = real(J_trans(3,:) * l3)
          call calc_mat_inv3(dble(r_trans(:,:,i,j,k,ispec_cpml)),temp,NDIM)
          r_trans_inv(:,1,i,j,k,ispec_cpml) = real(temp(1,:))
          r_trans_inv(:,2,i,j,k,ispec_cpml) = real(temp(2,:))
          r_trans_inv(:,3,i,j,k,ispec_cpml) = real(temp(3,:))
        enddo
      enddo
    enddo
  enddo
  
end subroutine calc_r_trans

subroutine get_J_local(xelm,yelm,zelm,dershape3D,J_local,J_local_inv)
  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,NDIM,ZERO
  implicit none
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ) :: J_local,&
          J_local_inv
  double precision, dimension(NDIM,NDIM) :: temp, temp_inv
  integer i,j,k,ia
  
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

        temp(:,:) = 0.0
        temp_inv(:,:) = 0.0

        do ia=1,NGNOD
          temp(1,1) = temp(1,1) + dershape3D(1,ia,i,j,k)*xelm(ia)
          temp(1,2) = temp(1,2) + dershape3D(2,ia,i,j,k)*xelm(ia)
          temp(1,3) = temp(1,3) + dershape3D(3,ia,i,j,k)*xelm(ia)

          temp(2,1) = temp(2,1) + dershape3D(1,ia,i,j,k)*yelm(ia)
          temp(2,2) = temp(2,2) + dershape3D(2,ia,i,j,k)*yelm(ia)
          temp(2,3) = temp(2,3) + dershape3D(3,ia,i,j,k)*yelm(ia)

          temp(3,1) = temp(3,1) + dershape3D(1,ia,i,j,k)*zelm(ia)
          temp(3,2) = temp(3,2) + dershape3D(2,ia,i,j,k)*zelm(ia)
          temp(3,3) = temp(3,3) + dershape3D(3,ia,i,j,k)*zelm(ia)
        enddo

        call calc_mat_inv3(temp,temp_inv,NDIM)
        J_local(:,:,i,j,k) = temp(:,:)
        J_local_inv(:,:,i,j,k) = temp_inv(:,:)
      enddo
    enddo
  enddo
end subroutine get_J_local

subroutine calc_mat_inv3(A,B,n)
  integer :: n
  double precision, dimension(n,n) :: A,B
  double precision :: detinv
  if (.not. (n==3)) then
    print *, 'wrong dimension'
    return
  endif
  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

end subroutine calc_mat_inv3

subroutine calc_mat_mul(A,B,C,n)
  integer :: n
  double precision, dimension(n,n) :: A,B,C
  integer :: i,j,k
  C(:,:) = 0.0
  do j = 1, n
    do i = 1, n
      do k = 1, n
        C(i,j) = C(i,j) + A(i,k) * B(k,j)
      enddo
    enddo 
  enddo
end subroutine calc_mat_mul
