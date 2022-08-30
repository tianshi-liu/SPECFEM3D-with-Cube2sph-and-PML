subroutine update_parameters_from_netcdf(model_fn, &
                x, y, z, rho, vp, vs, rho_new, vp_new, vs_new, &
                nx, ny, nz, nspec)
  !! subroutine to assign a netcdf model to a sphere-like mesh in SPECFEM3D
  !! Tianshi Liu, 2021.1
  implicit none
  integer :: nx, ny, nz, nspec, i, j, k, ispec
  integer, parameter :: SIZE_REAL = 4, SIZE_DOUBLE = 8
  integer, parameter :: CUSTOM_REAL = SIZE_REAL
  real(kind=CUSTOM_REAL),dimension(nx,ny,nz,nspec) :: rho, vp, vs, &
          rho_new, vp_new, vs_new
  real, parameter :: vs_min = 1000.0
  double precision, dimension(nx,ny,nz,nspec) :: x, y, z
  double precision :: r, theta, phi, lat, lon, depth, xmesh, ymesh, zmesh
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: PI_OVER_TWO = PI / 2.0d0
  double precision, parameter :: RADIANS_TO_DEGREES = 180.d0 / PI
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0
  double precision, parameter :: R_EARTH = 6371000.d0
  double precision, parameter :: R_EARTH_KM = R_EARTH / 1000.d0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real, dimension(:,:,:), allocatable :: vs_rec
  real, dimension(:), allocatable :: depth_grid
  real, dimension(:), allocatable :: lat_grid
  real, dimension(:), allocatable :: lon_grid
  real, dimension(2, 2, 2) :: vs_square
  integer :: idepth, ilat, ilon
  double precision :: ratio_depth, ratio_lat, ratio_lon

  !!!!netcdf variables!!
  !character(len = *), parameter :: model_fn = 'AK_EMB_2019JGRSE/MCMC_EMBerg_highres_2019_AK.nc' ! MODIFY this line to set path to your .nc file
  !character(len = *), parameter :: model_fn = 'model_Ward18/Alaska.ANT+RF.Ward.2018_kmps.nc' ! MODIFY this line to set path to your .nc file
  character(len = *) :: model_fn ! MODIFY this line to set path to your .nc file
  real :: lat_start, lon_start, ddepth, dlat, dlon
  integer :: ndepth, nlat, nlon
  integer :: ncid,dimid,dimid_lonlatd(3),var_num,&
             dimids_temp(1),dimids_3d(3),map_3d(3),map_3d_temp(3),shape_dim(3)
  real :: vs_fill
  integer :: len_char
  character (len=10) :: vs_unit
  !!!!!!!!!!!!!!!!!!!!!!
  
  !!!! load netcdf model
  call open_netcdf_model(trim(model_fn),ncid)
  
  ! get longitude dimension
  call inquire_dimension_netcdf_model(ncid,'longitude',dimid,nlon)
  call inquire_variable_netcdf_model(ncid,'longitude',1,var_num,dimids_temp)
  allocate(lon_grid(nlon))
  dimid_lonlatd(dimid) = 1
  shape_dim(dimid) = nlon
  call get_variable1d_netcdf_model(ncid, var_num, lon_grid, shape_dim(dimid))
  !print *, 'longitude dimension: ', nlon 
  
  ! get latitude dimension
  call inquire_dimension_netcdf_model(ncid,'latitude',dimid,nlat)
  call inquire_variable_netcdf_model(ncid,'latitude',1,var_num,dimids_temp)
  allocate(lat_grid(nlat))
  dimid_lonlatd(dimid) = 2
  shape_dim(dimid) = nlat
  call get_variable1d_netcdf_model(ncid, var_num, lat_grid, shape_dim(dimid))
  !print *, 'latitude dimension: ', nlat
  
  ! get depth dimension
  call inquire_dimension_netcdf_model(ncid,'depth',dimid,ndepth)
  call inquire_variable_netcdf_model(ncid,'depth',1,var_num,dimids_temp)
  allocate(depth_grid(ndepth))
  dimid_lonlatd(dimid) = 3
  shape_dim(dimid) = ndepth
  call get_variable1d_netcdf_model(ncid, var_num, depth_grid, shape_dim(dimid))
  !print *, 'depth dimension: ', ndepth

  ! get vs dataset
  call inquire_variable_netcdf_model(ncid, 'vs', 3, var_num, dimids_3d)
  allocate(vs_rec(nlon,nlat,ndepth))
  map_3d_temp(1) = 1
  map_3d_temp(2) = map_3d_temp(1) * shape_dim(dimids_3d(1))
  map_3d_temp(3) = map_3d_temp(2) * shape_dim(dimids_3d(2))
  map_3d(dimid_lonlatd(dimids_3d(1))) = map_3d_temp(1)
  map_3d(dimid_lonlatd(dimids_3d(2))) = map_3d_temp(2)
  map_3d(dimid_lonlatd(dimids_3d(3))) = map_3d_temp(3)
  !print *, shape_dim
  !print *, map_3d
  call get_variable3d_netcdf_model(ncid, var_num, vs_rec,(/nlon,nlat,ndepth/),&
                                   map_3d, vs_fill)
  !print *, vs_fill
  
  dlat = lat_grid(2) - lat_grid(1)
  dlon = lon_grid(2) - lon_grid(1)
  lat_start = lat_grid(1)
  lon_start = lon_grid(1)
  ! check if the unit for vs is km/s, convert to m/s
  !! TL: assume all EMC models are in km/s
  !call inquire_char_parameter_netcdf_model(ncid,'vs','units',vs_unit,len_char)
  !if (vs_unit(1:len_char)=='kmps') then
    vs_rec(:,:,:) = vs_rec(:,:,:) * 1000.d0 ! m/s
    vs_fill = vs_fill*1000.d0 - 1.0
  !endif

  call close_netcdf_model(ncid)
  !!!!!!!!!!!!!!!!!!!!!

  do ispec = 1,nspec
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          xmesh = x(i,j,k,ispec)
          ymesh = y(i,j,k,ispec)
          zmesh = z(i,j,k,ispec)
          !!! MODIFY the following transform !!!!!!
          !!! (xmesh,ymesh,zmesh) -> (lat,lon,depth)
          call xyz_2_rlatlon_dble(xmesh,ymesh,zmesh,r,lat,lon)
          if (lon > 180.0d0 ) lon = lon - 360.0d0
          depth = (ONE - r) * R_EARTH_KM
          !depth = R_EARTH_KM - r / 1000.0d0
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!! interpolation !!!!
          !do idepth = 1,ndepth
          !  if(depth_grid(idepth) > depth) exit
          !enddo
          idepth = 1
          do while (idepth .le. ndepth)
            if (depth_grid(idepth) > depth) exit
            idepth = idepth + 1
          enddo
          if(depth_grid(idepth) > depth) idepth = idepth - 1
          ! depth_grid(idepth) <= depth, depth_grid(idepth+1) >= depth 
          ! depth grid not necessarily regular
          ilat = floor((lat - lat_start) / dlat) + 1
          ilon = floor((lon - lon_start) / dlon) + 1
          if ((idepth .ge. ndepth) .or. &
              (ilat .le. 0) .or. (ilat .ge. nlat) .or. &
              (ilon .le. 0) .or. (ilon .ge. nlon)) then  !! out of grid
            vs_new(i,j,k,ispec) = vs(i,j,k,ispec)
          else
            if (idepth .le. 0) then
              vs_square(:,:,1) = vs_rec(ilon:ilon + 1, ilat:ilat + 1, 1)
              vs_square(:,:,2) = vs_rec(ilon:ilon + 1, ilat:ilat + 1, 1)
            else
              vs_square(:,:,:) = vs_rec(ilon:ilon + 1, ilat:ilat + 1, idepth:idepth + 1)
            endif
            if ((vs_square(1,1,1) .ge. vs_fill) .or. &
                (vs_square(2,1,1) .ge. vs_fill) .or. &
                (vs_square(1,2,1) .ge. vs_fill) .or. &
                (vs_square(2,2,1) .ge. vs_fill) .or. &
                (vs_square(1,1,2) .ge. vs_fill) .or. &
                (vs_square(2,1,2) .ge. vs_fill) .or. &
                (vs_square(1,2,2) .ge. vs_fill) .or. &
                (vs_square(2,2,2) .ge. vs_fill)) then !! undefined value
              vs_new(i,j,k,ispec) = vs(i,j,k,ispec)
            else
              if (idepth .le. 0) then
                ratio_depth = 0.0
              else
                ddepth = depth_grid(idepth+1) - depth_grid(idepth)
                ratio_depth = (depth - depth_grid(idepth)) / ddepth
              endif
              ratio_lat = (lat - lat_grid(ilat)) / dlat
              ratio_lon = (lon - lon_grid(ilon)) / dlon
              vs_new(i,j,k,ispec) = vs_square(1,1,1) * (ONE - ratio_depth) * (ONE - ratio_lat) * (ONE - ratio_lon) + &
                      vs_square(1,1,2) * ratio_depth * (ONE - ratio_lat) * (ONE - ratio_lon) + &
                      vs_square(1,2,1) * (ONE - ratio_depth) * ratio_lat * (ONE - ratio_lon) + &
                      vs_square(1,2,2) * ratio_depth * ratio_lat * (ONE - ratio_lon) + &
                      vs_square(2,1,1) * (ONE - ratio_depth) * (ONE - ratio_lat) * ratio_lon + &
                      vs_square(2,1,2) * ratio_depth * (ONE - ratio_lat) * ratio_lon + &
                      vs_square(2,2,1) * (ONE - ratio_depth) * ratio_lat * ratio_lon + &
                      vs_square(2,2,2) * ratio_depth * ratio_lat * ratio_lon
            end if
          end if
          ! check for bad vs values
          !if (vs_new(i,j,k,ispec) > vp(i,j,k,ispec)) then
          if (vs_new(i,j,k,ispec) < 1.0) then
            print *, 'bad vs value.'
            print *, 'new vs: ', vs_new(i,j,k,ispec), ', old vs: ', &
            vs(i,j,k,ispec), ', vp: ', vp(i,j,k,ispec)
            print *, 'lat: ', lat, ', lon: ', lon, ', depth: ', depth
            print *, 'ratio_depth:', ratio_depth, ' ratio_lat:', ratio_lat, &
                     ' ratio_lon:', ratio_lon
          endif
          if ((depth .ge. 150.0) .and. (vs_new(i,j,k,ispec) < 2000.0)) then
            print *, 'bad vs = ', vs_new(i,j,k,ispec), ' at ', lat, lon, depth
            print *, 'vs_square = ', vs_square(1,1,1), vs_square(1,1,2), &
                                     vs_square(1,2,1), vs_square(1,2,2), &
                                     vs_square(2,1,1), vs_square(2,1,2), &
                                     vs_square(2,2,1), vs_square(2,2,2)
          endif
          if (vs_new(i,j,k,ispec) < vs_min) then
            vs_new(i,j,k,ispec) = vs_min
          endif
          ! set bulk velocity constant
          !vp_new(i,j,k,ispec) = sqrt(vp(i,j,k,ispec) * vp(i,j,k,ispec) + 4.d0 * &
          !(vs_new(i,j,k,ispec) * vs_new(i,j,k,ispec) - vs(i,j,k,ispec) * &
          ! vs(i,j,k,ispec)) / 3.d0)
          ! set the Poisson ratio constant
          vp_new(i,j,k,ispec) = vp(i,j,k,ispec) * (vs_new(i,j,k,ispec) &
                                / vs(i,j,k,ispec))
          ! scale dlnrho=0.33dlnvs
          rho_new(i,j,k,ispec) = rho(i,j,k,ispec)*(1.0+0.33* &
                    (vs_new(i,j,k,ispec)-vs(i,j,k,ispec))/vs(i,j,k,ispec))
        enddo
      enddo
    enddo
  enddo
  !rho_new(:,:,:,:) = rho(:,:,:,:)
  !vp_new(:,:,:,:) = vp(:,:,:,:)
  !vs_new = vs
  deallocate(lon_grid, lat_grid, depth_grid, vs_rec)
end subroutine update_parameters_from_netcdf
