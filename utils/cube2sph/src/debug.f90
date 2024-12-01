subroutine wait_for_attach(myrank)
  implicit none
  integer :: myrank, ierr, pid
  logical, volatile :: escape = .false.
  pid = 0 !getpid()
  print *, 'PID ', pid, ' ready to attach for rank ', myrank
  do 
    call sleep(1)
    if (escape) exit
  enddo
  print *, myrank, ' has attached.'
end subroutine wait_for_attach
