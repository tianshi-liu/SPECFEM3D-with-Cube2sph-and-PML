subroutine read_netcdf_model(FILE_NAME, NDPTS, NLATS, NLONS, depth_grid, lat_grid, lon_grid, vs_rec)
  use netcdf
  implicit none
  !character (len = *), parameter :: FILE_NAME = "Alaska.ANT+RF.Ward.2018_kmps.nc"
  character (len = *) :: FILE_NAME
  integer :: ncid, dpt_varid, lat_varid, lon_varid, vs_varid
  integer, parameter :: NDIMS = 3 
  integer :: NDPTS, NLATS, NLONS
  character (len = *), parameter :: DPT_NAME = "depth"
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: VS_NAME = "vs"
  !real, dimension(NDPTS, NLATS, NLONS) :: vs_rec
  real, dimension(NLONS, NLATS, NDPTS) :: vs_rec
  real, dimension(NLATS) :: lat_grid
  real, dimension(NLONS) :: lon_grid
  real, dimension(NDPTS) :: depth_grid
  !integer :: start(NDIMS), ncount(NDIMS)
  !allocate(vs_rec(NDPTS, NLATS, NLONS))
  call check(nf90_open(FILE_NAME, nf90_nowrite, ncid)) ! open file
  call check(nf90_inq_varid(ncid, DPT_NAME, dpt_varid))
  call check(nf90_inq_varid(ncid, LAT_NAME, lat_varid))
  call check(nf90_inq_varid(ncid, LON_NAME, lon_varid)) ! query variable id
  call check(nf90_get_var(ncid, dpt_varid, depth_grid))
  call check(nf90_get_var(ncid, lat_varid, lat_grid))
  call check(nf90_get_var(ncid, lon_varid, lon_grid)) ! read grid data
  call check(nf90_inq_varid(ncid, VS_NAME, vs_varid))
  !start = (/1, 1, 1/)
  !ncount = (/NDPTS, NLATS, NLONS/)
  !ncount = (/NLONS, NLATS, NDPTS/)
  !call check(nf90_get_var(ncid, vs_varid, vs_rec, start=start, count=ncount))
  call check(nf90_get_var(ncid, vs_varid, vs_rec))
  call check(nf90_close(ncid))
  !deallocate(vs_rec)

contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check
end subroutine read_netcdf_model
