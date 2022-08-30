subroutine open_netcdf_model(file_name,ncid)
  use netcdf
  implicit none
  integer :: ncid,stat
  character (len=*) :: file_name
  stat = nf90_open(trim(file_name), nf90_nowrite, ncid)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
end subroutine open_netcdf_model

subroutine close_netcdf_model(ncid)
  use netcdf
  implicit none
  integer :: ncid,stat
  stat = nf90_close(ncid)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
end subroutine close_netcdf_model

subroutine inquire_dimension_netcdf_model(ncid, var_name, dimid, dim_len)
  use netcdf
  implicit none
  integer, intent(in) :: ncid
  character (len = *), intent(in) :: var_name
  integer, intent(out) :: dimid, dim_len
  ! local parameter
  integer :: stat
  character (len = 100) :: var_name_dummy
  stat = nf90_inq_dimid(ncid, trim(var_name), dimid)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
  stat = nf90_inquire_dimension(ncid, dimid, var_name_dummy, dim_len)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
end subroutine inquire_dimension_netcdf_model

subroutine inquire_variable_netcdf_model(ncid, var_name, ndims, var_num, dimids)
  use netcdf
  implicit none
  integer, intent(in) :: ncid, ndims
  character (len = *), intent(in) :: var_name
  integer, intent(out) :: var_num, dimids(ndims)
  ! local parameter
  integer :: stat, ndims_dummy, var_type
  character (len = 100) :: var_name_dummy
  stat = nf90_inq_varid(ncid, trim(var_name), var_num)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
  stat = nf90_inquire_variable(ncid, var_num, var_name_dummy, var_type, &
                               ndims_dummy, dimids)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
  if (var_type /= nf90_float) print *, trim(var_name)//' is not float type'
end subroutine inquire_variable_netcdf_model

subroutine get_variable1d_netcdf_model(ncid, var_num, var_val, var_shape)
  use netcdf
  implicit none
  integer, intent(in) :: ncid, var_num
  integer, intent(in) :: var_shape(1)
  real, intent(out) :: var_val(var_shape(1))
  ! local parameter
  integer :: stat
  stat = nf90_get_var(ncid, var_num, var_val)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
end subroutine get_variable1d_netcdf_model

subroutine get_variable3d_netcdf_model(ncid, var_num, var_val, var_shape, &
                                       map_3d, val_fill)
  use netcdf
  implicit none
  integer, intent(in) :: ncid, var_num
  integer, intent(in) :: var_shape(3),map_3d(3)
  real, intent(out) :: var_val(var_shape(1),var_shape(2),var_shape(3)),val_fill
  ! local parameter
  integer :: stat, no_fill
  stat = nf90_get_var(ncid, var_num, var_val, map=map_3d)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
  stat = nf90_inq_var_fill(ncid, var_num, no_fill, val_fill)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))

end subroutine get_variable3d_netcdf_model

subroutine inquire_int_parameter_netcdf_model(ncid,var_name,att_name,att_val,att_len)
  use netcdf
  implicit none
  integer :: ncid,att_len,stat,att_num,var_num,xtype,att_val
  character (len=100) :: att_name,var_name
  if (trim(var_name)=='global') then
    var_num = nf90_global
  else
    stat = nf90_inq_varid(ncid, var_name, var_num)
    if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
  endif
  stat = nf90_inquire_attribute(ncid, var_num, trim(att_name), xtype, att_len, att_num)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
  stat = nf90_get_att(ncid, var_num, trim(att_name), att_val)

end subroutine inquire_int_parameter_netcdf_model

subroutine inquire_char_parameter_netcdf_model(ncid,var_name,att_name,att_val,att_len)
  use netcdf
  implicit none
  integer :: ncid,att_len,stat,att_num,var_num,xtype
  character (len=100) :: att_name,var_name,att_val
  if (trim(var_name)=='global') then
    var_num = nf90_global
  else
    stat = nf90_inq_varid(ncid, var_name, var_num)
    if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
  endif
  stat = nf90_inquire_attribute(ncid, var_num, trim(att_name), xtype, att_len, att_num)
  if (stat /= nf90_noerr) print *, trim(nf90_strerror(stat))
  stat = nf90_get_att(ncid, var_num, trim(att_name), att_val)

end subroutine inquire_char_parameter_netcdf_model
