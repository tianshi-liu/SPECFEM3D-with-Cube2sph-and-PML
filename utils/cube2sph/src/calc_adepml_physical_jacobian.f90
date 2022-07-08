subroutine calc_adepml_physical_jacobian(num_face,xstore,ystore,&
                zstore,nspec,face_ispec,face_ijk,face_normal,face_jacobianw)
  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,GAUSSALPHA,GAUSSBETA,&
          NGLLX,NGLLY,NGNOD2D,NDIM2D,NGLLSQUARE
  use generate_databases_par, only: CPML_to_spec
  !use generate_databases_par, only: elmnts_ext_mesh
  implicit none
  !! local variables
  integer :: npts, nspec, ia, ispec, iface,i,j,k,igll
  integer :: num_face, iboun
  integer :: face_ispec(num_face), face_ijk(NDIM,NGLLSQUARE,num_face)
  real :: face_normal(NDIM,NGLLSQUARE,num_face), &
          face_jacobianw(NGLLSQUARE,num_face)
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,&
          ystore,zstore
  double precision :: shape2D(NGNOD2D,NGLLX,NGLLY)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision :: xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)
  double precision :: jacobian2Dw(NGLLSQUARE), normal(NDIM,NGLLSQUARE)
  double precision :: xref,yref,zref
  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)

  ! get the 2-D shape functions
  call get_shape2D(shape2D,dershape2D,xigll,yigll,NGLLX,NGLLY)

  do iface = 1, num_face
    i = face_ijk(1,(NGLLSQUARE+1)/2,iface)
    j = face_ijk(2,(NGLLSQUARE+1)/2,iface)
    k = face_ijk(3,(NGLLSQUARE+1)/2,iface)
    ispec = CPML_to_spec(face_ispec(iface))
    if (i .eq. 1) then !xmin
      iboun = 1
      xelm(1)=xstore(1,1,1,ispec)
      yelm(1)=ystore(1,1,1,ispec)
      zelm(1)=zstore(1,1,1,ispec)
      xelm(2)=xstore(1,NGLLY,1,ispec)
      yelm(2)=ystore(1,NGLLY,1,ispec)
      zelm(2)=zstore(1,NGLLY,1,ispec)
      xelm(3)=xstore(1,NGLLY,NGLLZ,ispec)
      yelm(3)=ystore(1,NGLLY,NGLLZ,ispec)
      zelm(3)=zstore(1,NGLLY,NGLLZ,ispec)
      xelm(4)=xstore(1,1,NGLLZ,ispec)
      yelm(4)=ystore(1,1,NGLLZ,ispec)
      zelm(4)=zstore(1,1,NGLLZ,ispec)
      xelm(5)=xstore(1,(NGLLY+1)/2,1,ispec)
      yelm(5)=ystore(1,(NGLLY+1)/2,1,ispec)
      zelm(5)=zstore(1,(NGLLY+1)/2,1,ispec)
      xelm(6)=xstore(1,NGLLY,(NGLLZ+1)/2,ispec)
      yelm(6)=ystore(1,NGLLY,(NGLLZ+1)/2,ispec)
      zelm(6)=zstore(1,NGLLY,(NGLLZ+1)/2,ispec)
      xelm(7)=xstore(1,(NGLLY+1)/2,NGLLZ,ispec)
      yelm(7)=ystore(1,(NGLLY+1)/2,NGLLZ,ispec)
      zelm(7)=zstore(1,(NGLLY+1)/2,NGLLZ,ispec)
      xelm(8)=xstore(1,1,(NGLLZ+1)/2,ispec)
      yelm(8)=ystore(1,1,(NGLLZ+1)/2,ispec)
      zelm(8)=zstore(1,1,(NGLLZ+1)/2,ispec)
      xelm(9)=xstore(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
      yelm(9)=ystore(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
      zelm(9)=zstore(1,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
    else if (i .eq. NGLLX) then  !xmax
      iboun = 2
      xelm(1)=xstore(NGLLX,1,1,ispec)
      yelm(1)=ystore(NGLLX,1,1,ispec)
      zelm(1)=zstore(NGLLX,1,1,ispec)
      xelm(2)=xstore(NGLLX,NGLLY,1,ispec)
      yelm(2)=ystore(NGLLX,NGLLY,1,ispec)
      zelm(2)=zstore(NGLLX,NGLLY,1,ispec)
      xelm(3)=xstore(NGLLX,NGLLY,NGLLZ,ispec)
      yelm(3)=ystore(NGLLX,NGLLY,NGLLZ,ispec)
      zelm(3)=zstore(NGLLX,NGLLY,NGLLZ,ispec)
      xelm(4)=xstore(NGLLX,1,NGLLZ,ispec)
      yelm(4)=ystore(NGLLX,1,NGLLZ,ispec)
      zelm(4)=zstore(NGLLX,1,NGLLZ,ispec)
      xelm(5)=xstore(NGLLX,(NGLLY+1)/2,1,ispec)
      yelm(5)=ystore(NGLLX,(NGLLY+1)/2,1,ispec)
      zelm(5)=zstore(NGLLX,(NGLLY+1)/2,1,ispec)
      xelm(6)=xstore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
      yelm(6)=ystore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
      zelm(6)=zstore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
      xelm(7)=xstore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
      yelm(7)=ystore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
      zelm(7)=zstore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
      xelm(8)=xstore(NGLLX,1,(NGLLZ+1)/2,ispec)
      yelm(8)=ystore(NGLLX,1,(NGLLZ+1)/2,ispec)
      zelm(8)=zstore(NGLLX,1,(NGLLZ+1)/2,ispec)
      xelm(9)=xstore(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
      yelm(9)=ystore(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
      zelm(9)=zstore(NGLLX,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
    else if (j .eq. 1) then !ymin
      iboun = 3
      xelm(1)=xstore(1,1,1,ispec)
      yelm(1)=ystore(1,1,1,ispec)
      zelm(1)=zstore(1,1,1,ispec)
      xelm(2)=xstore(NGLLX,1,1,ispec)
      yelm(2)=ystore(NGLLX,1,1,ispec)
      zelm(2)=zstore(NGLLX,1,1,ispec)
      xelm(3)=xstore(NGLLX,1,NGLLZ,ispec)
      yelm(3)=ystore(NGLLX,1,NGLLZ,ispec)
      zelm(3)=zstore(NGLLX,1,NGLLZ,ispec)
      xelm(4)=xstore(1,1,NGLLZ,ispec)
      yelm(4)=ystore(1,1,NGLLZ,ispec)
      zelm(4)=zstore(1,1,NGLLZ,ispec)
      xelm(5)=xstore((NGLLX+1)/2,1,1,ispec)
      yelm(5)=ystore((NGLLX+1)/2,1,1,ispec)
      zelm(5)=zstore((NGLLX+1)/2,1,1,ispec)
      xelm(6)=xstore(NGLLX,1,(NGLLZ+1)/2,ispec)
      yelm(6)=ystore(NGLLX,1,(NGLLZ+1)/2,ispec)
      zelm(6)=zstore(NGLLX,1,(NGLLZ+1)/2,ispec)
      xelm(7)=xstore((NGLLX+1)/2,1,NGLLZ,ispec)
      yelm(7)=ystore((NGLLX+1)/2,1,NGLLZ,ispec)
      zelm(7)=zstore((NGLLX+1)/2,1,NGLLZ,ispec)
      xelm(8)=xstore(1,1,(NGLLZ+1)/2,ispec)
      yelm(8)=ystore(1,1,(NGLLZ+1)/2,ispec)
      zelm(8)=zstore(1,1,(NGLLZ+1)/2,ispec)
      xelm(9)=xstore((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec)
      yelm(9)=ystore((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec)
      zelm(9)=zstore((NGLLX+1)/2,1,(NGLLZ+1)/2,ispec)
    else if (j .eq. NGLLX) then !ymax
      iboun = 4
      xelm(1)=xstore(1,NGLLY,1,ispec)
      yelm(1)=ystore(1,NGLLY,1,ispec)
      zelm(1)=zstore(1,NGLLY,1,ispec)
      xelm(2)=xstore(NGLLX,NGLLY,1,ispec)
      yelm(2)=ystore(NGLLX,NGLLY,1,ispec)
      zelm(2)=zstore(NGLLX,NGLLY,1,ispec)
      xelm(3)=xstore(NGLLX,NGLLY,NGLLZ,ispec)
      yelm(3)=ystore(NGLLX,NGLLY,NGLLZ,ispec)
      zelm(3)=zstore(NGLLX,NGLLY,NGLLZ,ispec)
      xelm(4)=xstore(1,NGLLY,NGLLZ,ispec)
      yelm(4)=ystore(1,NGLLY,NGLLZ,ispec)
      zelm(4)=zstore(1,NGLLY,NGLLZ,ispec)
      xelm(5)=xstore((NGLLX+1)/2,NGLLY,1,ispec)
      yelm(5)=ystore((NGLLX+1)/2,NGLLY,1,ispec)
      zelm(5)=zstore((NGLLX+1)/2,NGLLY,1,ispec)
      xelm(6)=xstore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
      yelm(6)=ystore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
      zelm(6)=zstore(NGLLX,NGLLY,(NGLLZ+1)/2,ispec)
      xelm(7)=xstore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
      yelm(7)=ystore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
      zelm(7)=zstore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
      xelm(8)=xstore(1,NGLLY,(NGLLZ+1)/2,ispec)
      yelm(8)=ystore(1,NGLLY,(NGLLZ+1)/2,ispec)
      zelm(8)=zstore(1,NGLLY,(NGLLZ+1)/2,ispec)
      xelm(9)=xstore((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec)
      yelm(9)=ystore((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec)
      zelm(9)=zstore((NGLLX+1)/2,NGLLY,(NGLLZ+1)/2,ispec)
    else if (k .eq. 1) then  !zmin
      iboun = 5
      xelm(1)=xstore(1,1,1,ispec)
      yelm(1)=ystore(1,1,1,ispec)
      zelm(1)=zstore(1,1,1,ispec)
      xelm(2)=xstore(NGLLX,1,1,ispec)
      yelm(2)=ystore(NGLLX,1,1,ispec)
      zelm(2)=zstore(NGLLX,1,1,ispec)
      xelm(3)=xstore(NGLLX,NGLLY,1,ispec)
      yelm(3)=ystore(NGLLX,NGLLY,1,ispec)
      zelm(3)=zstore(NGLLX,NGLLY,1,ispec)
      xelm(4)=xstore(1,NGLLY,1,ispec)
      yelm(4)=ystore(1,NGLLY,1,ispec)
      zelm(4)=zstore(1,NGLLY,1,ispec)
      xelm(5)=xstore((NGLLX+1)/2,1,1,ispec)
      yelm(5)=ystore((NGLLX+1)/2,1,1,ispec)
      zelm(5)=zstore((NGLLX+1)/2,1,1,ispec)
      xelm(6)=xstore(NGLLX,(NGLLY+1)/2,1,ispec)
      yelm(6)=ystore(NGLLX,(NGLLY+1)/2,1,ispec)
      zelm(6)=zstore(NGLLX,(NGLLY+1)/2,1,ispec)
      xelm(7)=xstore((NGLLX+1)/2,NGLLY,1,ispec)
      yelm(7)=ystore((NGLLX+1)/2,NGLLY,1,ispec)
      zelm(7)=zstore((NGLLX+1)/2,NGLLY,1,ispec)
      xelm(8)=xstore(1,(NGLLY+1)/2,1,ispec)
      yelm(8)=ystore(1,(NGLLY+1)/2,1,ispec)
      zelm(8)=zstore(1,(NGLLY+1)/2,1,ispec)
      xelm(9)=xstore((NGLLX+1)/2,(NGLLY+1)/2,1,ispec)
      yelm(9)=ystore((NGLLX+1)/2,(NGLLY+1)/2,1,ispec)
      zelm(9)=zstore((NGLLX+1)/2,(NGLLY+1)/2,1,ispec)
    else if (k .eq. NGLLX) then !zmax
      iboun = 6
      xelm(1)=xstore(1,1,NGLLZ,ispec)
      yelm(1)=ystore(1,1,NGLLZ,ispec)
      zelm(1)=zstore(1,1,NGLLZ,ispec)
      xelm(2)=xstore(NGLLX,1,NGLLZ,ispec)
      yelm(2)=ystore(NGLLX,1,NGLLZ,ispec)
      zelm(2)=zstore(NGLLX,1,NGLLZ,ispec)
      xelm(3)=xstore(NGLLX,NGLLY,NGLLZ,ispec)
      yelm(3)=ystore(NGLLX,NGLLY,NGLLZ,ispec)
      zelm(3)=zstore(NGLLX,NGLLY,NGLLZ,ispec)
      xelm(4)=xstore(1,NGLLY,NGLLZ,ispec)
      yelm(4)=ystore(1,NGLLY,NGLLZ,ispec)
      zelm(4)=zstore(1,NGLLY,NGLLZ,ispec)
      xelm(5)=xstore((NGLLX+1)/2,1,NGLLZ,ispec)
      yelm(5)=ystore((NGLLX+1)/2,1,NGLLZ,ispec)
      zelm(5)=zstore((NGLLX+1)/2,1,NGLLZ,ispec)
      xelm(6)=xstore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
      yelm(6)=ystore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
      zelm(6)=zstore(NGLLX,(NGLLY+1)/2,NGLLZ,ispec)
      xelm(7)=xstore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
      yelm(7)=ystore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
      zelm(7)=zstore((NGLLX+1)/2,NGLLY,NGLLZ,ispec)
      xelm(8)=xstore(1,(NGLLY+1)/2,NGLLZ,ispec)
      yelm(8)=ystore(1,(NGLLY+1)/2,NGLLZ,ispec)
      zelm(8)=zstore(1,(NGLLY+1)/2,NGLLZ,ispec)
      xelm(9)=xstore((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec)
      yelm(9)=ystore((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec)
      zelm(9)=zstore((NGLLX+1)/2,(NGLLY+1)/2,NGLLZ,ispec)
    endif
    !call compute_jacobian_2D(xelm,yelm,zelm,dershape2D,jacobian2Dw,normal,&
    !        face_ijk(:,:,face_ispec(iface)),iboun)
    call compute_jacobian_2D(xelm,yelm,zelm,dershape2D,jacobian2Dw,normal,&
            face_ijk(:,:,iface),iboun)
    do igll = 1, NGLLSQUARE
      i = face_ijk(1,igll,iface)
      j = face_ijk(2,igll,iface)
      k = face_ijk(3,igll,iface)
      if (iboun .eq. 1) then
        xref = xstore(i,j,k,ispec) - xstore(i+1,j,k,ispec)
        yref = ystore(i,j,k,ispec) - ystore(i+1,j,k,ispec)
        zref = zstore(i,j,k,ispec) - zstore(i+1,j,k,ispec)
        jacobian2Dw(igll) = jacobian2Dw(igll) * wxgll(j) * wygll(k)
      elseif (iboun .eq. 2) then
        xref = xstore(i,j,k,ispec) - xstore(i-1,j,k,ispec)
        yref = ystore(i,j,k,ispec) - ystore(i-1,j,k,ispec)
        zref = zstore(i,j,k,ispec) - zstore(i-1,j,k,ispec)
        jacobian2Dw(igll) = jacobian2Dw(igll) * wxgll(j) * wygll(k)
      elseif (iboun .eq. 3) then
        xref = xstore(i,j,k,ispec) - xstore(i,j+1,k,ispec)
        yref = ystore(i,j,k,ispec) - ystore(i,j+1,k,ispec)
        zref = zstore(i,j,k,ispec) - zstore(i,j+1,k,ispec)
        jacobian2Dw(igll) = jacobian2Dw(igll) * wxgll(i) * wygll(k)
      elseif (iboun .eq. 4) then
        xref = xstore(i,j,k,ispec) - xstore(i,j-1,k,ispec)
        yref = ystore(i,j,k,ispec) - ystore(i,j-1,k,ispec)
        zref = zstore(i,j,k,ispec) - zstore(i,j-1,k,ispec)
        jacobian2Dw(igll) = jacobian2Dw(igll) * wxgll(i) * wygll(k)
      elseif (iboun .eq. 5) then
        xref = xstore(i,j,k,ispec) - xstore(i,j,k+1,ispec)
        yref = ystore(i,j,k,ispec) - ystore(i,j,k+1,ispec)
        zref = zstore(i,j,k,ispec) - zstore(i,j,k+1,ispec)
        jacobian2Dw(igll) = jacobian2Dw(igll) * wxgll(i) * wygll(j)
      elseif (iboun .eq. 6) then
        xref = xstore(i,j,k,ispec) - xstore(i,j,k-1,ispec)
        yref = ystore(i,j,k,ispec) - ystore(i,j,k-1,ispec)
        zref = zstore(i,j,k,ispec) - zstore(i,j,k-1,ispec)
        jacobian2Dw(igll) = jacobian2Dw(igll) * wxgll(i) * wygll(j)
      endif
      if ((xref*normal(1,igll)+yref*normal(2,igll)+zref*normal(3,igll))<0) &
              print *, iboun, ' wrong normal direction'
    enddo
    face_normal(:,:,iface) = real(normal(:,:))
    face_jacobianw(:,iface) = real(jacobian2Dw(:))
  enddo
end subroutine calc_adepml_physical_jacobian


subroutine get_shape2D(shape2D,dershape2D,xigll,yigll,NGLLA,NGLLB)

  use constants

  implicit none

! generic routine that accepts any polynomial degree in each direction

  integer NGLLA,NGLLB

  double precision xigll(NGLLA)
  double precision yigll(NGLLB)

! 2D shape functions and their derivatives
  double precision shape2D(NGNOD2D,NGLLA,NGLLB)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLA,NGLLB)

  integer i,j,ia

! location of the nodes of the 2D quadrilateral elements
  double precision xi,eta
  double precision xi_map,eta_map

! for checking the 2D shape functions
  !double precision sumshape,sumdershapexi,sumdershapeeta

! check that the parameter file is correct
  !if (NGNOD /= 8 .and. NGNOD /= 27) &
  !  call exit_MPI(myrank,'volume elements should have 8 or 27 control nodes')
  !if (NGNOD2D /= 4 .and. NGNOD2D /= 9) &
  !  call exit_MPI(myrank,'surface elements should have 4 or 9 control nodes')

  if (NGNOD2D == 4) then

    ! generate the 2D shape functions and their derivatives (4 nodes)
    do i = 1,NGLLA

      xi = xigll(i)

      do j = 1,NGLLB

        eta = yigll(j)

        ! map coordinates to [0,1]
        xi_map = 0.5d0 * (xi + 1.d0)
        eta_map = 0.5d0 * (eta + 1.)

        ! corner nodes
        shape2D(1,i,j) = (1.d0 - xi_map)*(1.d0 - eta_map)
        shape2D(2,i,j) = xi_map*(1.d0 - eta_map)
        shape2D(3,i,j) = xi_map*eta_map
        shape2D(4,i,j) = (1.d0 - xi_map)*eta_map

        dershape2D(1,1,i,j) = 0.25d0 * (eta - 1.d0)
        dershape2D(2,1,i,j) = 0.25d0 * (xi - 1.d0)

        dershape2D(1,2,i,j) = 0.25d0 * (1.d0 - eta)
        dershape2D(2,2,i,j) = 0.25d0 * (-1.d0 - xi)

        dershape2D(1,3,i,j) = 0.25d0 * (1.d0 + eta)
        dershape2D(2,3,i,j) = 0.25d0 * (1.d0 + xi)

        dershape2D(1,4,i,j) = 0.25d0 * (- 1.d0 - eta)
        dershape2D(2,4,i,j) = 0.25d0 * (1.d0 - xi)

      enddo
    enddo

  else

    ! note: put further initialization for ngnod2d == 9 into subroutine
    !       to avoid compilation errors in case ngnod2d == 4
    call get_shape2D_9(shape2D,dershape2D,xigll,yigll,NGLLA,NGLLB)

  endif


  end subroutine get_shape2D


subroutine get_shape2D_9(shape2D,dershape2D,xigll,yigll,NGLLA,NGLLB)

  use constants

  implicit none

! generic routine that accepts any polynomial degree in each direction

  !integer :: NGNOD2D
  integer :: NGLLA,NGLLB

  double precision xigll(NGLLA)
  double precision yigll(NGLLB)

! 2D shape functions and their derivatives
  double precision shape2D(NGNOD2D,NGLLA,NGLLB)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLA,NGLLB)

  integer i,j

! location of the nodes of the 2D quadrilateral elements
  double precision xi,eta
  double precision l1xi,l2xi,l3xi,l1eta,l2eta,l3eta
  double precision l1pxi,l2pxi,l3pxi,l1peta,l2peta,l3peta

  ! check that the parameter file is correct
  if (NGNOD2D /= 9) stop 'surface elements should have 9 control nodes'

  ! generate the 2D shape functions and their derivatives (9 nodes)
  do i = 1,NGLLA

    xi = xigll(i)

    l1xi = HALF*xi*(xi - ONE)
    l2xi = ONE - xi**2
    l3xi = HALF*xi*(xi + ONE)

    l1pxi = xi - HALF
    l2pxi = -TWO * xi
    l3pxi = xi + HALF

    do j = 1,NGLLB

      eta = yigll(j)

      l1eta = HALF*eta*(eta - ONE)
      l2eta = ONE - eta**2
      l3eta = HALF*eta*(eta + ONE)

      l1peta = eta - HALF
      l2peta = -TWO * eta
      l3peta = eta + HALF

!   corner nodes

      shape2D(1,i,j) = l1xi*l1eta
      shape2D(2,i,j) = l3xi*l1eta
      shape2D(3,i,j) = l3xi*l3eta
      shape2D(4,i,j) = l1xi*l3eta

      dershape2D(1,1,i,j) = l1pxi*l1eta
      dershape2D(1,2,i,j) = l3pxi*l1eta
      dershape2D(1,3,i,j) = l3pxi*l3eta
      dershape2D(1,4,i,j) = l1pxi*l3eta

      dershape2D(2,1,i,j) = l1xi*l1peta
      dershape2D(2,2,i,j) = l3xi*l1peta
      dershape2D(2,3,i,j) = l3xi*l3peta
      dershape2D(2,4,i,j) = l1xi*l3peta

!   midside nodes

      shape2D(5,i,j) = l2xi*l1eta
      shape2D(6,i,j) = l3xi*l2eta
      shape2D(7,i,j) = l2xi*l3eta
      shape2D(8,i,j) = l1xi*l2eta

      dershape2D(1,5,i,j) = l2pxi*l1eta
      dershape2D(1,6,i,j) = l3pxi*l2eta
      dershape2D(1,7,i,j) = l2pxi*l3eta
      dershape2D(1,8,i,j) = l1pxi*l2eta

      dershape2D(2,5,i,j) = l2xi*l1peta
      dershape2D(2,6,i,j) = l3xi*l2peta
      dershape2D(2,7,i,j) = l2xi*l3peta
      dershape2D(2,8,i,j) = l1xi*l2peta

!   center node

      shape2D(9,i,j) = l2xi*l2eta

      dershape2D(1,9,i,j) = l2pxi*l2eta
      dershape2D(2,9,i,j) = l2xi*l2peta

    enddo
  enddo

  end subroutine get_shape2D_9

subroutine compute_jacobian_2D(xelm,yelm,zelm,dershape2D,jacobian2D,normal,&
                this_face_ijk,iboun)

  use constants

  implicit none

! generic routine that accepts any polynomial degree in each direction

  !integer ispecb,NGLLA,NGLLB,NSPEC2DMAX_AB
  integer :: this_face_ijk(NDIM,NGLLSQUARE)

  double precision xelm(NGNOD2D),yelm(NGNOD2D),zelm(NGNOD2D)
  double precision dershape2D(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  double precision jacobian2D(NGLLSQUARE)
  double precision normal(NDIM,NGLLSQUARE)

  integer i,j,ia,igll,iboun
  double precision xxi,xeta,yxi,yeta,zxi,zeta
  double precision unx,uny,unz,jacobian

  do igll = 1, NGLLSQUARE
    if ((iboun .eq. 1) .or. (iboun .eq. 2)) then
      i = this_face_ijk(2,igll)
      j = this_face_ijk(3,igll)
    elseif ((iboun .eq. 3) .or. (iboun .eq. 4)) then
      i = this_face_ijk(1,igll)
      j = this_face_ijk(3,igll)
    elseif ((iboun .eq. 5) .or. (iboun .eq. 6)) then
      i = this_face_ijk(1,igll)
      j = this_face_ijk(2,igll)
    endif

    xxi=ZERO
    xeta=ZERO
    yxi=ZERO
    yeta=ZERO
    zxi=ZERO
    zeta=ZERO

    do ia = 1,NGNOD2D
      xxi=xxi+dershape2D(1,ia,i,j)*xelm(ia)
      xeta=xeta+dershape2D(2,ia,i,j)*xelm(ia)
      yxi=yxi+dershape2D(1,ia,i,j)*yelm(ia)
      yeta=yeta+dershape2D(2,ia,i,j)*yelm(ia)
      zxi=zxi+dershape2D(1,ia,i,j)*zelm(ia)
      zeta=zeta+dershape2D(2,ia,i,j)*zelm(ia)
    enddo

    !   calculate the unnormalized normal to the boundary
    unx=yxi*zeta-yeta*zxi
    uny=zxi*xeta-zeta*xxi
    unz=xxi*yeta-xeta*yxi
    jacobian=dsqrt(unx**2+uny**2+unz**2)

    !if (jacobian <= ZERO) call exit_MPI(myrank,'2D Jacobian undefined')

    !   normalize normal vector and store surface Jacobian
    jacobian2D(igll)=jacobian
    normal(1,igll)=unx/jacobian
    normal(2,igll)=uny/jacobian
    normal(3,igll)=unz/jacobian


  enddo
  if ((iboun .eq. 1) .or. (iboun .eq. 4) .or. (iboun .eq. 5)) &
            normal(:,:) = - normal(:,:)

  end subroutine compute_jacobian_2D
