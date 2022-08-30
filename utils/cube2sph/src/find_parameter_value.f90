subroutine find_parameter_value(key, val, f_unit)
  character(len=*) key, val
  integer :: f_unit
  integer :: ier, ind
  character(len=256) :: string
  rewind(unit=f_unit)
  ier = 0
  do while (ier == 0)
    read(f_unit,"(a)",iostat=ier) string
    if (ier /= 0) exit
    string = adjustl(string)
    !ind = index(string, trim(key))
    if ((index(string, trim(key)//' ') /= 1) .and. &
        (index(string, trim(key)//'=') /= 1)) cycle !line does not start with key
    ind = index(string, '#')
    if (ind > 1) then
      ! has a comment string, discard that
      string = string(1:ind-1)
    endif
    ind = index(string, '=')
    if (ind == 0) stop 'incorrect format'
    string = string(ind+1: len_trim(string))
    val = adjustl(string)
    return
  enddo 
  print *, 'key not found:' // trim(key)
  stop 1
end subroutine find_parameter_value
