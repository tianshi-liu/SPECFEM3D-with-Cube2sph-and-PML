program write_slice_vtk
  use constants, only: R_EARTH, MAX_STRING_LEN, NDIM
  implicit none 
  character(len=MAX_STRING_LEN) :: xyzfile, vtkfile
  integer :: i, j, ier,nv = 8
  double precision, dimension(:,:), allocatable :: coords_cube, coords
  double precision :: xi_start, xi_end, eta_start, eta_end, dxi, deta, depth, &
                      depth_start, depth_end, ddepth, &
                      center_lat = 90.0, center_lon = 0.0, rotation_azi = 0.0
  integer :: nhor, nver, npts, ipt, nelem, ielem
  integer , dimension(:,:), allocatable :: ipt_map
  xyzfile = 'slice_vert1.xyz'
  vtkfile = 'slice_vert1.vtk'
  xi_start = -1111949.25
  xi_end = 1111949.25
  eta_start = -1111949.25
  eta_end = 1111949.25
  depth_start = 0.0
  depth_end = 220000.0
  nhor = 320
  nver = 30
  dxi = (xi_end - xi_start) / (nhor * 1.0)
  deta = (eta_end - eta_start) / (nhor * 1.0)
  ddepth = (depth_end - depth_start) / (nver)
  !depth = 500.0
  nelem = nhor * nver
  allocate(ipt_map(0:nv-1,0:nelem-1))
  do i = 0, nhor-1
    do j = 0, nver-1
      ielem = j * nhor + i
      ipt_map(0,ielem) = j * (nhor + 1) + i
      ipt_map(1,ielem) = j * (nhor + 1) + i + 1
      ipt_map(2,ielem) = (j + 1) * (nhor + 1) + i + 1
      ipt_map(3,ielem) = (j + 1) * (nhor + 1) + i
      ipt_map(4,ielem) = j * nhor + i + (nhor + 1) * (nver + 1)
      ipt_map(5,ielem) = j * (nhor + 1) + i + 1 + (2 * nhor + 1) * (nver + 1)
      ipt_map(6,ielem) = (j + 1) * nhor + i + (nhor + 1) * (nver + 1)
      ipt_map(7,ielem) = j * (nhor + 1) + i + (2 * nhor + 1) * (nver + 1)
    enddo
  enddo
  npts = (nhor+1)*(nver+1) + (nhor+1)*nver + nhor*(nver+1)
  allocate(coords_cube(NDIM,0:npts-1))
  ipt = 0
  do j = 0, nver
    do i = 0, nhor
      coords_cube(1,ipt) = xi_start + i * dxi
      coords_cube(2,ipt) = eta_start + i * deta
      coords_cube(3,ipt) = -(depth_start + j * ddepth)
      ipt = ipt + 1
    enddo
  enddo
  do j = 0, nver
    do i = 0, nhor-1
      coords_cube(1,ipt) = xi_start + (i + 0.5) * dxi
      coords_cube(2,ipt) = eta_start + (i + 0.5) * deta
      coords_cube(3,ipt) = -(depth_start + j * ddepth)
      ipt = ipt + 1
    enddo
  enddo
  do j = 0, nver-1
    do i = 0, nhor
      coords_cube(1,ipt) = xi_start + i * dxi
      coords_cube(2,ipt) = eta_start + i * deta
      coords_cube(3,ipt) = -(depth_start + (j + 0.5) * ddepth)
      ipt = ipt + 1
    enddo
  enddo
  if (ipt /= npts) print *, 'incorrect indexing'
  allocate(coords(NDIM,0:npts-1))
  call cube2sph_trans(coords_cube, coords, npts, R_EARTH, &
                      center_lat,center_lon,rotation_azi)
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
  write(77,'(a,2i12)') 'CELLS', nelem, nelem * (nv+1)
  !write(77,'(a,2i12)') 'CELLS', nelem, nelem * 5
  do ielem = 0, nelem-1
    write(77,'(9i12)') 8, ipt_map(0:7,ielem)
    !write(77,'(9i12)') 4, ipt_map(0:3,ielem)
  enddo
  write(77,'(a)') ''
  write(77,'(a,i12)') 'CELL_TYPES', nelem
  do ielem = 0, nelem-1
    write(77,'(i12)') 23
  enddo
  close(55)
  deallocate(ipt_map,coords,coords_cube)
end program write_slice_vtk
