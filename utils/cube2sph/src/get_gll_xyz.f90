subroutine get_gll_xyz(nodes_coords,npts,xstore,ystore,zstore,jacobian3D,nspec)
  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,NDIM,GAUSSALPHA,GAUSSBETA
  use generate_databases_par, only: elmnts_ext_mesh
  implicit none
  !! local variables
  integer :: npts, nspec, ia, ispec
  double precision, dimension(NDIM,npts) :: nodes_coords
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,&
          ystore,zstore,jacobian3D
  double precision :: shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! get the 3-D shape functions
  call get_shape3D(shape3D,dershape3D,xigll,yigll,zigll)

  ! point locations
  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0


  do ispec = 1, nspec
    do ia = 1,NGNOD
      xelm(ia) = nodes_coords(1,elmnts_ext_mesh(ia,ispec))
      yelm(ia) = nodes_coords(2,elmnts_ext_mesh(ia,ispec))
      zelm(ia) = nodes_coords(3,elmnts_ext_mesh(ia,ispec))
    enddo
    call calc_coords(xstore(:,:,:,ispec),ystore(:,:,:,ispec),&
            zstore(:,:,:,ispec), xelm,yelm,zelm,shape3D)
    call calc_jacobian(jacobian3D(:,:,:,ispec),xelm,yelm,zelm,dershape3D)
  enddo

end subroutine get_gll_xyz

!subroutine get_shape3D(shape3D,dershape3D,xigll,yigll,zigll)
!  use constants
!  implicit none
!
!! Gauss-Lobatto-Legendre points of integration
!  double precision xigll(NGLLX)
!  double precision yigll(NGLLY)
!  double precision zigll(NGLLZ)
!
!! 3D shape functions and their derivatives
!  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
!  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
!
!  integer i,j,k,ia
!
!! location of the nodes of the 3D hexahedra elements
!  double precision xi,eta,gamma
!  double precision ra1,ra2,rb1,rb2,rc1,rc2
!
!! for checking the 3D shape functions
!  !double precision sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma
!
!  !double precision, parameter :: ONE_EIGHTH = 0.125d0
!
!! ***
!! *** create 3D shape functions and jacobian
!! ***
!
!  do i=1,NGLLX
!    do j=1,NGLLY
!      do k=1,NGLLZ
!
!        xi = xigll(i)
!        eta = yigll(j)
!        gamma = zigll(k)
!
!        !--- case of a 3D 8-node element (Dhatt-Touzot p. 115)
!        if (NGNOD == 8) then
!
!          ra1 = one + xi
!          ra2 = one - xi
!
!          rb1 = one + eta
!          rb2 = one - eta
!
!          rc1 = one + gamma
!          rc2 = one - gamma
!
!          shape3D(1,i,j,k) = ONE_EIGHTH*ra2*rb2*rc2
!          shape3D(2,i,j,k) = ONE_EIGHTH*ra1*rb2*rc2
!          shape3D(3,i,j,k) = ONE_EIGHTH*ra1*rb1*rc2
!          shape3D(4,i,j,k) = ONE_EIGHTH*ra2*rb1*rc2
!          shape3D(5,i,j,k) = ONE_EIGHTH*ra2*rb2*rc1
!          shape3D(6,i,j,k) = ONE_EIGHTH*ra1*rb2*rc1
!          shape3D(7,i,j,k) = ONE_EIGHTH*ra1*rb1*rc1
!          shape3D(8,i,j,k) = ONE_EIGHTH*ra2*rb1*rc1
!          
!          dershape3D(1,1,i,j,k) = - ONE_EIGHTH*rb2*rc2
!          dershape3D(1,2,i,j,k) = ONE_EIGHTH*rb2*rc2
!          dershape3D(1,3,i,j,k) = ONE_EIGHTH*rb1*rc2
!          dershape3D(1,4,i,j,k) = - ONE_EIGHTH*rb1*rc2
!          dershape3D(1,5,i,j,k) = - ONE_EIGHTH*rb2*rc1
!          dershape3D(1,6,i,j,k) = ONE_EIGHTH*rb2*rc1
!          dershape3D(1,7,i,j,k) = ONE_EIGHTH*rb1*rc1
!          dershape3D(1,8,i,j,k) = - ONE_EIGHTH*rb1*rc1
!
!          dershape3D(2,1,i,j,k) = - ONE_EIGHTH*ra2*rc2
!          dershape3D(2,2,i,j,k) = - ONE_EIGHTH*ra1*rc2
!          dershape3D(2,3,i,j,k) = ONE_EIGHTH*ra1*rc2
!          dershape3D(2,4,i,j,k) = ONE_EIGHTH*ra2*rc2
!          dershape3D(2,5,i,j,k) = - ONE_EIGHTH*ra2*rc1
!          dershape3D(2,6,i,j,k) = - ONE_EIGHTH*ra1*rc1
!          dershape3D(2,7,i,j,k) = ONE_EIGHTH*ra1*rc1
!          dershape3D(2,8,i,j,k) = ONE_EIGHTH*ra2*rc1
!
!          dershape3D(3,1,i,j,k) = - ONE_EIGHTH*ra2*rb2
!          dershape3D(3,2,i,j,k) = - ONE_EIGHTH*ra1*rb2
!          dershape3D(3,3,i,j,k) = - ONE_EIGHTH*ra1*rb1
!          dershape3D(3,4,i,j,k) = - ONE_EIGHTH*ra2*rb1
!          dershape3D(3,5,i,j,k) = ONE_EIGHTH*ra2*rb2
!          dershape3D(3,6,i,j,k) = ONE_EIGHTH*ra1*rb2
!          dershape3D(3,7,i,j,k) = ONE_EIGHTH*ra1*rb1
!          dershape3D(3,8,i,j,k) = ONE_EIGHTH*ra2*rb1
!        else  
!          call get_shape3D_27(shape3D,dershape3D,xi,eta,gamma,i,j,k)
!        endif
!      enddo
!    enddo
!  enddo
!
!end subroutine get_shape3D
!
!!--- case of a 3D 27-node element
!
!  subroutine get_shape3D_27(shape3D,dershape3D,xi,eta,gamma,i,j,k)
!
!  use constants
!
!  implicit none
!
!  integer :: i,j,k
!
!! 3D shape functions and their derivatives
!  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
!  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
!
!! location of the nodes of the 3D hexahedra elements
!  double precision xi,eta,gamma
!  double precision l1xi,l2xi,l3xi,l1eta,l2eta,l3eta,l1gamma,l2gamma,l3gamma
!  double precision l1pxi,l2pxi,l3pxi,l1peta,l2peta,l3peta,l1pgamma,l2pgamma,l3pgamma
!
!  l1xi=HALF*xi*(xi-ONE)
!  l2xi=ONE-xi**2
!  l3xi=HALF*xi*(xi+ONE)
!
!  l1pxi=xi-HALF
!  l2pxi=-TWO*xi
!  l3pxi=xi+HALF
!
!  l1eta=HALF*eta*(eta-ONE)
!  l2eta=ONE-eta**2
!  l3eta=HALF*eta*(eta+ONE)
!
!  l1peta=eta-HALF
!  l2peta=-TWO*eta
!  l3peta=eta+HALF
!
!  l1gamma=HALF*gamma*(gamma-ONE)
!  l2gamma=ONE-gamma**2
!  l3gamma=HALF*gamma*(gamma+ONE)
!
!  l1pgamma=gamma-HALF
!  l2pgamma=-TWO*gamma
!  l3pgamma=gamma+HALF
!
!  ! corner nodes
!  shape3D(1,i,j,k)=l1xi*l1eta*l1gamma
!  shape3D(2,i,j,k)=l3xi*l1eta*l1gamma
!  shape3D(3,i,j,k)=l3xi*l3eta*l1gamma
!  shape3D(4,i,j,k)=l1xi*l3eta*l1gamma
!  shape3D(5,i,j,k)=l1xi*l1eta*l3gamma
!  shape3D(6,i,j,k)=l3xi*l1eta*l3gamma
!  shape3D(7,i,j,k)=l3xi*l3eta*l3gamma
!  shape3D(8,i,j,k)=l1xi*l3eta*l3gamma
!
!  dershape3D(1,1,i,j,k)=l1pxi*l1eta*l1gamma
!  dershape3D(1,2,i,j,k)=l3pxi*l1eta*l1gamma
!  dershape3D(1,3,i,j,k)=l3pxi*l3eta*l1gamma
!  dershape3D(1,4,i,j,k)=l1pxi*l3eta*l1gamma
!  dershape3D(1,5,i,j,k)=l1pxi*l1eta*l3gamma
!  dershape3D(1,6,i,j,k)=l3pxi*l1eta*l3gamma
!  dershape3D(1,7,i,j,k)=l3pxi*l3eta*l3gamma
!  dershape3D(1,8,i,j,k)=l1pxi*l3eta*l3gamma
!
!  dershape3D(2,1,i,j,k)=l1xi*l1peta*l1gamma
!  dershape3D(2,2,i,j,k)=l3xi*l1peta*l1gamma
!  dershape3D(2,3,i,j,k)=l3xi*l3peta*l1gamma
!  dershape3D(2,4,i,j,k)=l1xi*l3peta*l1gamma
!  dershape3D(2,5,i,j,k)=l1xi*l1peta*l3gamma
!  dershape3D(2,6,i,j,k)=l3xi*l1peta*l3gamma
!  dershape3D(2,7,i,j,k)=l3xi*l3peta*l3gamma
!  dershape3D(2,8,i,j,k)=l1xi*l3peta*l3gamma
!
!  dershape3D(3,1,i,j,k)=l1xi*l1eta*l1pgamma
!  dershape3D(3,2,i,j,k)=l3xi*l1eta*l1pgamma
!  dershape3D(3,3,i,j,k)=l3xi*l3eta*l1pgamma
!  dershape3D(3,4,i,j,k)=l1xi*l3eta*l1pgamma
!  dershape3D(3,5,i,j,k)=l1xi*l1eta*l3pgamma
!  dershape3D(3,6,i,j,k)=l3xi*l1eta*l3pgamma
!  dershape3D(3,7,i,j,k)=l3xi*l3eta*l3pgamma
!  dershape3D(3,8,i,j,k)=l1xi*l3eta*l3pgamma
!
!  ! midside nodes
!  shape3D(9,i,j,k)=l2xi*l1eta*l1gamma
!  shape3D(10,i,j,k)=l3xi*l2eta*l1gamma
!  shape3D(11,i,j,k)=l2xi*l3eta*l1gamma
!  shape3D(12,i,j,k)=l1xi*l2eta*l1gamma
!  shape3D(13,i,j,k)=l1xi*l1eta*l2gamma
!  shape3D(14,i,j,k)=l3xi*l1eta*l2gamma
!  shape3D(15,i,j,k)=l3xi*l3eta*l2gamma
!  shape3D(16,i,j,k)=l1xi*l3eta*l2gamma
!  shape3D(17,i,j,k)=l2xi*l1eta*l3gamma
!  shape3D(18,i,j,k)=l3xi*l2eta*l3gamma
!  shape3D(19,i,j,k)=l2xi*l3eta*l3gamma
!  shape3D(20,i,j,k)=l1xi*l2eta*l3gamma
!
!  dershape3D(1,9,i,j,k)=l2pxi*l1eta*l1gamma
!  dershape3D(1,10,i,j,k)=l3pxi*l2eta*l1gamma
!  dershape3D(1,11,i,j,k)=l2pxi*l3eta*l1gamma
!  dershape3D(1,12,i,j,k)=l1pxi*l2eta*l1gamma
!  dershape3D(1,13,i,j,k)=l1pxi*l1eta*l2gamma
!  dershape3D(1,14,i,j,k)=l3pxi*l1eta*l2gamma
!  dershape3D(1,15,i,j,k)=l3pxi*l3eta*l2gamma
!  dershape3D(1,16,i,j,k)=l1pxi*l3eta*l2gamma
!  dershape3D(1,17,i,j,k)=l2pxi*l1eta*l3gamma
!  dershape3D(1,18,i,j,k)=l3pxi*l2eta*l3gamma
!  dershape3D(1,19,i,j,k)=l2pxi*l3eta*l3gamma
!  dershape3D(1,20,i,j,k)=l1pxi*l2eta*l3gamma
!
!  dershape3D(2,9,i,j,k)=l2xi*l1peta*l1gamma
!  dershape3D(2,10,i,j,k)=l3xi*l2peta*l1gamma
!  dershape3D(2,11,i,j,k)=l2xi*l3peta*l1gamma
!  dershape3D(2,12,i,j,k)=l1xi*l2peta*l1gamma
!  dershape3D(2,13,i,j,k)=l1xi*l1peta*l2gamma
!  dershape3D(2,14,i,j,k)=l3xi*l1peta*l2gamma
!  dershape3D(2,15,i,j,k)=l3xi*l3peta*l2gamma
!  dershape3D(2,16,i,j,k)=l1xi*l3peta*l2gamma
!  dershape3D(2,17,i,j,k)=l2xi*l1peta*l3gamma
!  dershape3D(2,18,i,j,k)=l3xi*l2peta*l3gamma
!  dershape3D(2,19,i,j,k)=l2xi*l3peta*l3gamma
!  dershape3D(2,20,i,j,k)=l1xi*l2peta*l3gamma
!
!  dershape3D(3,9,i,j,k)=l2xi*l1eta*l1pgamma
!  dershape3D(3,10,i,j,k)=l3xi*l2eta*l1pgamma
!  dershape3D(3,11,i,j,k)=l2xi*l3eta*l1pgamma
!  dershape3D(3,12,i,j,k)=l1xi*l2eta*l1pgamma
!  dershape3D(3,13,i,j,k)=l1xi*l1eta*l2pgamma
!  dershape3D(3,14,i,j,k)=l3xi*l1eta*l2pgamma
!  dershape3D(3,15,i,j,k)=l3xi*l3eta*l2pgamma
!  dershape3D(3,16,i,j,k)=l1xi*l3eta*l2pgamma
!  dershape3D(3,17,i,j,k)=l2xi*l1eta*l3pgamma
!  dershape3D(3,18,i,j,k)=l3xi*l2eta*l3pgamma
!  dershape3D(3,19,i,j,k)=l2xi*l3eta*l3pgamma
!  dershape3D(3,20,i,j,k)=l1xi*l2eta*l3pgamma
!
!  ! side center nodes
!  shape3D(21,i,j,k)=l2xi*l2eta*l1gamma
!  shape3D(22,i,j,k)=l2xi*l1eta*l2gamma
!  shape3D(23,i,j,k)=l3xi*l2eta*l2gamma
!  shape3D(24,i,j,k)=l2xi*l3eta*l2gamma
!  shape3D(25,i,j,k)=l1xi*l2eta*l2gamma
!  shape3D(26,i,j,k)=l2xi*l2eta*l3gamma
!
!  dershape3D(1,21,i,j,k)=l2pxi*l2eta*l1gamma
!  dershape3D(1,22,i,j,k)=l2pxi*l1eta*l2gamma
!  dershape3D(1,23,i,j,k)=l3pxi*l2eta*l2gamma
!  dershape3D(1,24,i,j,k)=l2pxi*l3eta*l2gamma
!  dershape3D(1,25,i,j,k)=l1pxi*l2eta*l2gamma
!  dershape3D(1,26,i,j,k)=l2pxi*l2eta*l3gamma
!
!  dershape3D(2,21,i,j,k)=l2xi*l2peta*l1gamma
!  dershape3D(2,22,i,j,k)=l2xi*l1peta*l2gamma
!  dershape3D(2,23,i,j,k)=l3xi*l2peta*l2gamma
!  dershape3D(2,24,i,j,k)=l2xi*l3peta*l2gamma
!  dershape3D(2,25,i,j,k)=l1xi*l2peta*l2gamma
!  dershape3D(2,26,i,j,k)=l2xi*l2peta*l3gamma
!
!  dershape3D(3,21,i,j,k)=l2xi*l2eta*l1pgamma
!  dershape3D(3,22,i,j,k)=l2xi*l1eta*l2pgamma
!  dershape3D(3,23,i,j,k)=l3xi*l2eta*l2pgamma
!  dershape3D(3,24,i,j,k)=l2xi*l3eta*l2pgamma
!  dershape3D(3,25,i,j,k)=l1xi*l2eta*l2pgamma
!  dershape3D(3,26,i,j,k)=l2xi*l2eta*l3pgamma
!
!  ! center node
!  shape3D(27,i,j,k)=l2xi*l2eta*l2gamma
!
!  dershape3D(1,27,i,j,k)=l2pxi*l2eta*l2gamma
!  dershape3D(2,27,i,j,k)=l2xi*l2peta*l2gamma
!  dershape3D(3,27,i,j,k)=l2xi*l2eta*l2pgamma
!
!  end subroutine get_shape3D_27

  subroutine calc_jacobian(jacobian3D,xelm,yelm,zelm,dershape3D)
  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,NDIM,ZERO
  implicit none
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGNOD) :: xelm,yelm,zelm

  integer i,j,k,ia
  double precision xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision jacobian3D(NGLLX,NGLLY,NGLLZ)

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

      xxi = ZERO
      xeta = ZERO
      xgamma = ZERO
      yxi = ZERO
      yeta = ZERO
      ygamma = ZERO
      zxi = ZERO
      zeta = ZERO
      zgamma = ZERO

      do ia=1,NGNOD
        xxi = xxi + dershape3D(1,ia,i,j,k)*xelm(ia)
        xeta = xeta + dershape3D(2,ia,i,j,k)*xelm(ia)
        xgamma = xgamma + dershape3D(3,ia,i,j,k)*xelm(ia)

        yxi = yxi + dershape3D(1,ia,i,j,k)*yelm(ia)
        yeta = yeta + dershape3D(2,ia,i,j,k)*yelm(ia)
        ygamma = ygamma + dershape3D(3,ia,i,j,k)*yelm(ia)

        zxi = zxi + dershape3D(1,ia,i,j,k)*zelm(ia)
        zeta = zeta + dershape3D(2,ia,i,j,k)*zelm(ia)
        zgamma = zgamma + dershape3D(3,ia,i,j,k)*zelm(ia)
      enddo

      jacobian3D(i,j,k) = xxi*(yeta*zgamma-ygamma*zeta) - &
                          xeta*(yxi*zgamma-ygamma*zxi) + &
                          xgamma*(yxi*zeta-yeta*zxi)
      enddo
    enddo
  enddo
  end subroutine calc_jacobian

  subroutine calc_coords(x_elem,y_elem,z_elem,xelm,yelm,zelm,shape3D)

  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,ZERO

  implicit none

  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: x_elem,y_elem,z_elem

  !local
  integer i,j,k,ia
  double precision xmesh,ymesh,zmesh

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

      xmesh = ZERO
      ymesh = ZERO
      zmesh = ZERO

      do ia=1,NGNOD
        xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
        ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
        zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
      enddo

      x_elem(i,j,k) = xmesh
      y_elem(i,j,k) = ymesh
      z_elem(i,j,k) = zmesh

      enddo
    enddo
  enddo

  end subroutine calc_coords
