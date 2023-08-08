program write_slice_vtk
  use constants, only: R_EARTH, MAX_STRING_LEN, NDIM, ONE, TWO
  use specfem_par, only: ONE_CRUST, ELLIPTICITY,nspl,rspl,espl,espl2
  implicit none 
  character(len=MAX_STRING_LEN) :: xyzfile, vtkfile
  integer :: i, j, ier,nv = 8
  double precision, dimension(:,:), allocatable :: coords_cube, coords
  double precision :: xi_start, xi_end, eta_start, eta_end, dxi, deta, depth, &
                    center_lat = 62.5, center_lon = -151.0, rotation_azi = 20.0
  double precision :: r, x, y, z, theta, phi, cost, p20, ell, factor
  integer :: nxi, neta, npts, ipt, nelem, ielem
  integer , dimension(:,:), allocatable :: ipt_map
  ONE_CRUST = .true.
  ELLIPTICITY = .true.
  xyzfile = 'slice_90000m.xyz'
  vtkfile = 'slice_90000m.vtk'
  xi_start = -1223000.0
  xi_end = 1223000.0
  eta_start = -1223000.0
  eta_end = 1223000.0
  nxi = 100
  neta = 100
  dxi = (xi_end - xi_start) / (nxi * 1.0)
  deta = (eta_end - eta_start) / (neta * 1.0)
  depth = 90000.0
  nelem = nxi * neta
  allocate(ipt_map(0:nv-1,0:nelem-1))
  do i = 0, nxi-1
    do j = 0, neta-1
      ielem = j * nxi + i
      ipt_map(0,ielem) = j * (nxi + 1) + i
      ipt_map(1,ielem) = j * (nxi + 1) + i + 1
      ipt_map(2,ielem) = (j + 1) * (nxi + 1) + i + 1
      ipt_map(3,ielem) = (j + 1) * (nxi + 1) + i
      ipt_map(4,ielem) = j * nxi + i + (nxi + 1) * (neta + 1)
      ipt_map(5,ielem) = j * (nxi + 1) + i + 1 + (2 * nxi + 1) * (neta + 1)
      ipt_map(6,ielem) = (j + 1) * nxi + i + (nxi + 1) * (neta + 1)
      ipt_map(7,ielem) = j * (nxi + 1) + i + (2 * nxi + 1) * (neta + 1)
    enddo
  enddo
  npts = (nxi+1)*(neta+1) + (nxi+1)*neta + nxi*(neta+1)
  allocate(coords_cube(NDIM,0:npts-1))
  ipt = 0
  do j = 0, neta
    do i = 0, nxi
      coords_cube(1,ipt) = xi_start + i * dxi
      coords_cube(2,ipt) = eta_start + j * deta
      coords_cube(3,ipt) = -depth
      ipt = ipt + 1
    enddo
  enddo
  do j = 0, neta
    do i = 0, nxi-1
      coords_cube(1,ipt) = xi_start + (i + 0.5) * dxi
      coords_cube(2,ipt) = eta_start + j * deta
      coords_cube(3,ipt) = -depth
      ipt = ipt + 1
    enddo
  enddo
  do j = 0, neta-1
    do i = 0, nxi
      coords_cube(1,ipt) = xi_start + i * dxi
      coords_cube(2,ipt) = eta_start + (j + 0.5) * deta
      coords_cube(3,ipt) = -depth
      ipt = ipt + 1
    enddo
  enddo
  if (ipt /= npts) print *, 'incorrect indexing'
  allocate(coords(NDIM,0:npts-1))
  call cube2sph_trans(coords_cube, coords, npts, R_EARTH, &
                      center_lat,center_lon,rotation_azi)
  ! sets up spline coefficients for ellipticity
  if (ELLIPTICITY) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  if (ELLIPTICITY) then
    do ipt = 0, npts-1
      !call get_ellipticity_single_point(coords(1,ipt), coords(2,ipt), &
      !                                  coords(3,ipt),nspl,rspl,espl,espl2)
      x = coords(1,ipt) / R_EARTH
      y = coords(2,ipt) / R_EARTH
      z = coords(3,ipt) / R_EARTH
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      cost = dcos(theta)
      ! this is the Legendre polynomial of degree two, P2(cos(theta)), see the
      ! discussion above eq (14.4) in Dahlen and Tromp (1998)
      p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

      ! get ellipticity using spline evaluation
      call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

      ! this is eq (14.4) in Dahlen and Tromp (1998)
      factor = ONE-(TWO/3.0d0)*ell*p20
      coords(:,ipt) = coords(:,ipt) * factor
    enddo
  endif
  open(55,file=trim(xyzfile),status='unknown',iostat=ier, action='write')
  do ipt = 0, npts-1
    write(55, '(3e18.6)') coords(1,ipt), coords(2,ipt), coords(3,ipt)
  enddo
  close(55)
  open(77,file=trim(vtkfile),status='unknown',iostat=ier, action='write')
  write(77,'(a)') '# vtk DataFile Version 3.1'
  write(77,'(a)') 'material model VTK file'
  write(77,'(a)') 'ASCII'
  write(77,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(77,'(a,i15,a)') 'POINTS ', npts, ' float'
  do ipt = 0, npts-1
    write(77,'(3f18.6)') coords(1,ipt), coords(2,ipt), coords(3,ipt)
  enddo
  write(77,'(a)') ''
  !write(77,'(a,2i12)') 'CELLS', nelem, nelem * (nv+1)
  write(77,'(a,2i12)') 'CELLS', nelem, nelem * 5
  do ielem = 0, nelem-1
    !write(77,'(9i12)') 8, ipt_map(0:7,ielem)
    write(77,'(9i12)') 4, ipt_map(0:3,ielem)
  enddo
  write(77,'(a)') ''
  write(77,'(a,i12)') 'CELL_TYPES', nelem
  do ielem = 0, nelem-1
    write(77,'(i12)') 9
  enddo
  close(55)
  deallocate(ipt_map,coords,coords_cube)
end program write_slice_vtk
