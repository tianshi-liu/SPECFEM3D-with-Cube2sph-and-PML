program create_ckbd
  use mpi
  use meshfem3D_models_par, only: myrank
  use constants_solver
  use specfem_par, only: &
      ONE_CRUST,rspl,espl,espl2,nspl,ibathy_topo,TOPOGRAPHY
  implicit none
  integer :: ier
  integer, dimension(:,:,:,:), allocatable :: ibool
  real, dimension(:,:,:,:),allocatable :: v_read, x,y,z, dv, v_new
  real, dimension(:), allocatable :: xstore, ystore, zstore
  !integer, parameter :: MAX_STRING_LEN = 512  ! constants.h
  character(len=MAX_STRING_LEN) :: file_ckbd,dirin,dirout,prname_lp,filename, &
                                   par_name
  character(len=MAX_STRING_LEN), parameter :: LOCAL_PATH = 'DATABASES_MPI'
  !integer, parameter :: NGLLX = 5
  !integer, parameter :: NGLLY = NGLLX
  !integer, parameter :: NGLLZ = NGLLX
  !integer, parameter :: IIN = 40,IOUT = 41
  integer :: nspec, i, j, k, ispec, nspec_irregular, nglob, iglob
  integer :: np,ip
  real, dimension(:), allocatable :: scalep,xp,yp,zp,&
                                                 sigmahp,sigmavp
  real :: dvmin,dvmax,dvmin_all,dvmax_all,exp_val
  double precision :: latp,lonp,dp,theta,phi,elevation,r0,ell,&
                      dcost,p20,radius
  logical :: ELLIPTICITY = .true.
  ONE_CRUST = .true.
  TOPOGRAPHY = .false.
  allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
  call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  ibathy_topo(:,:) = 0
  call read_topo_bathy_file(ibathy_topo)
 
  call get_command_argument(1,file_ckbd)
  call get_command_argument(2,dirin)
  call get_command_argument(3,dirout)
  call get_command_argument(4,par_name)

  open(unit=IIN,file=trim(file_ckbd),action='read',form='formatted',iostat=ier)
  read(IIN, *) np
  allocate(xp(np),yp(np),zp(np),sigmahp(np),sigmavp(np),scalep(np))
  do ip = 1, np
    read(IIN, *) latp,lonp,dp
    read(IIN, *) sigmahp(ip),sigmavp(ip),scalep(ip)
    ! limits longitude to [0.0,360.0]
    if (lonp < 0.d0 ) lonp = lonp + 360.d0
    if (lonp > 360.d0 ) lonp = lonp - 360.d0
    call lat_2_geocentric_colat_dble(latp,theta)   
    phi = lonp * DEGREES_TO_RADIANS
    call reduce(theta,phi)
    r0 = R_UNIT_SPHERE
    if (TOPOGRAPHY) then
      call get_topo_bathy(latp,lonp,elevation,ibathy_topo)
      r0 = r0 + elevation/R_EARTH
    endif
    dp = dp / 1000.0 ! convert to km for ellipticity calculation
    if (ELLIPTICITY) then
      dcost = dcos(theta)
      ! this is the Legendre polynomial of degree two, P2(cos(theta)), see the
      ! discussion above eq (14.4) in Dahlen and Tromp (1998)
      p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
      radius = r0 - dp*1000.0d0/R_EARTH
      ! get ellipticity using spline evaluation
      call spline_evaluation(rspl,espl,espl2,nspl,radius,ell)
      ! this is eq (14.4) in Dahlen and Tromp (1998)
      r0 = r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
    endif
    r0 = r0 - dp*1000.0d0/R_EARTH
    xp(ip) = real(r0*dsin(theta)*dcos(phi)*R_EARTH)
    yp(ip) = real(r0*dsin(theta)*dsin(phi)*R_EARTH)
    zp(ip) = real(r0*dcos(theta)*R_EARTH)
    !print *, r0*R_EARTH, theta*RADIANS_TO_DEGREES, phi*RADIANS_TO_DEGREES
  enddo

  call MPI_Init(ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
 
  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)// '/' //'proc',myrank,'_'
  ! reading in number of elements
  if (myrank == 0) write(*,*) '     reading in: external_mesh.bin'
  filename = prname_lp(1:len_trim(prname_lp))//'external_mesh.bin'
  open(unit=IIN,file=trim(filename),status='unknown',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening database proc######_external_mesh.bin'
  read(IIN) nspec
  allocate(ibool(NGLLX, NGLLY, NGLLZ, nspec), &
           x(NGLLX, NGLLY, NGLLZ, nspec), &
           y(NGLLX, NGLLY, NGLLZ, nspec), &
           z(NGLLX, NGLLY, NGLLZ, nspec), stat=ier)
  if (ier/=0) stop 'error allocating xyz array'
  read(IIN) nglob
  allocate(xstore(nglob), ystore(nglob), zstore(nglob), stat=ier)
  if (ier/=0) stop 'error allocating xyz store array'
  read(IIN) nspec_irregular
  read(IIN) ibool
  read(IIN) xstore
  read(IIN) ystore
  read(IIN) zstore
  close(IIN)
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          x(i,j,k,ispec) = xstore(iglob)
          y(i,j,k,ispec) = ystore(iglob)
          z(i,j,k,ispec) = zstore(iglob)
        enddo
      enddo
    enddo
  enddo
  deallocate(xstore, ystore, zstore, ibool)

  allocate(dv(NGLLX,NGLLY,NGLLZ,nspec))
  dv(:,:,:,:) = 0.0
  do ip = 1, np
    if (myrank == 0) print *, xp(ip), yp(ip), zp(ip), sigmahp(ip), &
                              sigmavp(ip), scalep(ip)
    do ispec = 1, nspec
      do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX
        call get_gaussian(x(i,j,k,ispec),y(i,j,k,ispec),z(i,j,k,ispec),&
                          xp(ip),yp(ip),zp(ip),sigmahp(ip),sigmavp(ip),&
                          exp_val)
        dv(i,j,k,ispec) = dv(i,j,k,ispec) + exp_val * scalep(ip)
      enddo;enddo;enddo
    enddo  
  enddo

  dvmin = minval(dv(:,:,:,:))
  dvmax = maxval(dv(:,:,:,:))

  dvmin_all = HUGEVAL
  dvmax_all = - HUGEVAL

  call min_all_all_cr(dvmin,dvmin_all)
  call max_all_all_cr(dvmax,dvmax_all)
  if (myrank == 0) print *, dvmin_all, dvmax_all

  write(prname_lp,'(a,i6.6,a)') trim(dirin)// '/' //'proc',myrank,'_'
  allocate(v_read(NGLLX,NGLLY,NGLLZ,nspec), v_new(NGLLX,NGLLY,NGLLZ,nspec))
  filename = prname_lp(1:len_trim(prname_lp))//trim(par_name)//'.bin'
  if (myrank == 0) print *, 'reading ', filename
  open(unit=IIN,file=trim(filename),status='old',action='read',&
       form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading input file'
  endif

  read(IIN) v_read
  close(IIN)
  
  do ispec = 1, nspec  
    do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX
      v_new(i,j,k,ispec) = v_read(i,j,k,ispec) * exp(dv(i,j,k,ispec))
      if (.not.((v_new(i,j,k,ispec)>1.0).and.&
                      (v_new(i,j,k,ispec)<10000.0))) then
        print *, myrank, 'bad dv', dv(i,j,k,ispec)
      endif
    enddo;enddo;enddo
  enddo

  !! write model
  write(prname_lp,'(a,i6.6,a)') trim(dirout)// '/' //'proc',myrank,'_'
  if (myrank == 0) print *, 'writing to ', prname_lp(1:len_trim(prname_lp))//trim(par_name)//'.bin'
  open(unit=IOUT,file=prname_lp(1:len_trim(prname_lp))//trim(par_name)//'.bin',&
       status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening output file'
  write(IOUT) v_new
  close(IOUT)
  !! write perturbation
  write(prname_lp,'(a,i6.6,a)') trim(dirout)// '/' //'proc',myrank,'_'
  if (myrank == 0) &
    print *, 'writing to ', prname_lp(1:len_trim(prname_lp))//'d'//trim(par_name)//'.bin'
  open(unit=IOUT,&
       file=prname_lp(1:len_trim(prname_lp))//'d'//trim(par_name)//'.bin',&
       status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'error opening perturbation file'
  write(IOUT) dv
  close(IOUT)
 
  deallocate(xp,yp,zp,sigmahp,sigmavp,scalep)
  deallocate(v_read,v_new,dv)
  call MPI_Finalize(ier)

end program create_ckbd

subroutine get_gaussian(x1,y1,z1,x0,y0,z0,sigmah,sigmav,exp_val)
  implicit none
  real,intent(in) :: x1,y1,z1,x0,y0,z0,sigmah,sigmav
  real,intent(out) :: exp_val
  
  ! local parameters
  real :: dist_h_sq,dist_v_sq
  real :: val
  real, parameter :: R_EARTH = 6371000.0
  real :: r0, r1, theta, ratio
  r1 = sqrt( x1*x1 + y1*y1 + z1*z1 )
  r0 = sqrt( x0*x0 + y0*y0 + z0*z0 )
  dist_v_sq = (r0 - r1) * (r0 - r1)
  ratio = (x0*x1 + y0*y1 + z0*z1) / (r0 * r1)

  ! checks boundaries of ratio (due to numerical inaccuracies)
  if (ratio > 1.0) then
    ratio = 1.0
  else if (ratio < -1.0) then
    ratio = -1.0
  endif
  theta = acos( ratio )
  dist_h_sq = (R_EARTH * theta) * (R_EARTH * theta)
  val = - dist_h_sq / sigmah / sigmah - dist_v_sq / sigmav / sigmav
  if (val < -20.0) then
    exp_val = 0.0
  else
    exp_val = exp(val) 
  endif
end subroutine get_gaussian
