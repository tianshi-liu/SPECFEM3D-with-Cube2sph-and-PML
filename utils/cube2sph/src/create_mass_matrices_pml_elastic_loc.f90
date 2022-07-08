subroutine create_volume_matrices_pml_elastic_loc(nspec,jacobian3D,&
           nspec_cpml,rvolume_loc)
  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,GAUSSALPHA,GAUSSBETA
  use generate_databases_par, only: CPML_to_spec
  implicit none
  integer :: nspec,nspec_cpml,ispec_CPML,ispec, i,j,k
  !integer :: ibool_CPML(NGLLX,NGLLY,NGLLZ,nspec_cpml)
  real :: rvolume_loc(NGLLX,NGLLY,NGLLZ,nspec_cpml)
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision :: weight
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec_cpml) :: jacobian3D
  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  

  rvolume_loc(:,:,:,:) = 0.0
  
  do ispec_CPML = 1, nspec_cpml
    ispec = CPML_to_spec(ispec_CPML)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          !iglob_CPML = ibool_CPML(i,j,k,ispec_CPML)
          weight = wxgll(i)*wygll(j)*wzgll(k)
          rvolume_loc(i,j,k,ispec_CPML) = &
              real(jacobian3D(i,j,k,ispec)*weight)
        enddo
      enddo
    enddo
  enddo
end subroutine create_volume_matrices_pml_elastic_loc
