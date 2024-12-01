module interp3D
  use constants
  use interp3D_par
  use locate_points
  double precision                           :: x_target, y_target, z_target
  double precision                           :: x_found,  y_found,  z_found
  double precision                           :: xi_found, eta_found, gamma_found
  integer                                    :: ispec_selected, domain
  double precision                           :: final_distance_squared
contains
!use to locate the interp points xi,eta,gamma in mesh
subroutine interp_initialize()
  implicit none
  integer  ::  ipts
  double precision  :: points_location(2)
  !debug one point with point_min_allprocs
  !x_target=1002.d0
  !y_target=-1001.d0
  !z_target=-999.d0
  allocate(hxi_store(npts,NGLLX),stat=ier)
  allocate(heta_store(npts,NGLLX),stat=ier)
  allocate(hgamma_store(npts,NGLLX),stat=ier)
  allocate(points_location_store(npts,NGLLX),stat=ier)
  allocate(ispec_selected_store(npts),stat=ier)
  allocate(xin_found(npts),stat=ier)
  allocate(yin_found(npts),stat=ier)
  allocate(zin_found(npts),stat=ier)
  call setup_GLL_points()
  call get_elem_minmaxsize_glob()
    
  do ipts=1,npts
    x_target=xin(ipts)
    y_target=yin(ipts)
    z_target=zin(ipts)
    call locate_points_in_mesh(x_target, y_target, z_target, elemsize_max_glob, & 
                                ispec_selected, xi_found, eta_found, gamma_found, &
                                x_found, y_found, z_found, domain, final_distance_squared)
    x_xi(ipts)=xi_found
    y_eta(ipts)=eta_found
    z_gamma(ipts)=gamma_found
    xin_found(ipts)=x_found
    yin_found(ipts)=y_found
    zin_found(ipts)=z_found
    dist_xyz_newxyz(ipts)=final_distance_squared
    call lagrange_any(xi_found,NGLLX,xigll,hxi,hpxi)
    call lagrange_any(eta_found,NGLLY,yigll,heta,hpeta)
    call lagrange_any(gamma_found,NGLLZ,zigll,hgamma,hpgamma)
    hxi_store(ipts,:) = hxi(:)
    heta_store(ipts,:) = heta(:)
    hgamma_store(ipts,:) = hgamma(:)
    ispec_selected_store(ipts)=ispec_selected
    call point_min_allprocs(final_distance_squared,points_location)
    points_location_store(ipts,:)=real(points_location(:))
    !debug make sure the dist for each is close to zero
    ! if(rank == int(points_location_store(ipts,2))) then
    !   print*,'ipts,dist,rank,ispec_selected',ipts,dist_xyz_newxyz(ipts),rank,ispec_selected
    ! endif
  enddo
  
    
end subroutine interp_initialize


subroutine write_profile(COOR_PATH,SUPRESS_UTM)
  use mpi_f08
  use ieee_arithmetic,only : ieee_value,ieee_quiet_nan
  implicit none
  character(len=MAX_STR_LEN) :: OUT_NAME,SUPRESS_UTM,COOR_PATH
  real(kind=8)  :: rlon,rlat,rx,ry
  double precision :: dist_dummy,dx_dummy,dy_dummy,dz_dummy
  integer :: ipts,i
  real, allocatable :: myprofile(:,:),myprofile_temp(:,:)
  integer,allocatable :: flag(:),flag_temp(:),nspec_all(:)
  allocate(myprofile(npts,5),myprofile_temp(npts,5),flag(npts),flag_temp(npts),nspec_all(nproc))
  myprofile(:,:) = 0.
  flag(:) = 0

  if(rank == 0) then 
    print*,'writing profiles ...'
  endif

  do ipts=1,npts 
      !debug make sure the foubd points is very close the true locations
      if(rank == int(points_location_store(ipts,2))) then
          dist_dummy= dist_xyz_newxyz(ipts)
          dx_dummy=(xin(ipts)-xin_found(ipts))*(xin(ipts)-xin_found(ipts))
          dy_dummy=(yin(ipts)-yin_found(ipts))*(yin(ipts)-yin_found(ipts))
          dz_dummy=(zin(ipts)-zin_found(ipts))*(zin(ipts)-zin_found(ipts))
          if(dist_dummy>100 .or. dx_dummy>100 .or. dy_dummy>100 .or. dz_dummy>100 )then
            flag(ipts) = 0
          else
            flag(ipts) = 1
            myprofile(ipts,1) = x_xi(ipts)
            myprofile(ipts,2) = y_eta(ipts)
            myprofile(ipts,3) = z_gamma(ipts)
            myprofile(ipts,4) = rank 
            myprofile(ipts,5) = ispec_selected_store(ipts)
          endif
          ! the points are find very accurately
          ! print*,'ipts,x,x_f,dist',ipts,xin(ipts),xin_found(ipts),dist_xyz_newxyz(ipts)
          ! print*,'ipts,y,y_f,dist',ipts,yin(ipts),yin_found(ipts),dist_xyz_newxyz(ipts)
          ! print*,'ipts,z,z_f,dist',ipts,zin(ipts),zin_found(ipts),dist_xyz_newxyz(ipts)
      endif
  enddo

  call MPI_Reduce(myprofile,myprofile_temp,npts*5,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call MPI_Reduce(flag,flag_temp,npts,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(NSPEC_AB,1,MPI_INTEGER,nspec_all,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  if(rank == 0) then 
    myprofile = myprofile_temp
    flag = flag_temp
    OUT_NAME = trim(COOR_PATH) // '.loc'
    open(IOUT,file=trim(OUT_NAME))
    do ipts=1,npts 
      rx=xin(ipts)
      ry=yin(ipts)
      if(SUPRESS_UTM=='false')then
        call utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,1) ! geo2utm
        rx = real(rlon,kind=CUSTOM_REAL)
        ry = real(rlat,kind=CUSTOM_REAL)
        !print*,rlon,rlat
      endif
      if(flag(ipts) == 0) myprofile(ipts,4) = -1
      i = myprofile(ipts,4)
      if(i == -1) i = 1
      write(IOUT,'(6(G0,1x),I0,1x,I0,1x,I0)')rx,ry,zin(ipts) * 0.001,myprofile(ipts,1:3), &
                                      int(myprofile(ipts,4:5)),nspec_all(i)
    enddo

    close(IOUT)
  endif

  deallocate(myprofile,myprofile_temp,flag,flag_temp)

  call MPI_Barrier(MPI_COMM_WORLD)

end subroutine write_profile


subroutine interp()
  use mpi_f08
  implicit none
  integer  ::  ipts,i,j,k,ispec
  double precision :: hlagrange,dist_dummy,dx_dummy,dy_dummy,dz_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: model_out_dummy
  allocate(model_out_dummy(npts),stat=ier)
  model_out_dummy(:)=0.d0
  model_out(:)=0.d0
  do ipts=1,npts
    !debug make sure the foubd points is very close the true locations
    if(rank == int(points_location_store(ipts,2))) then
        dist_dummy= dist_xyz_newxyz(ipts)
        dx_dummy=(xin(ipts)-xin_found(ipts))*(xin(ipts)-xin_found(ipts))
        dy_dummy=(yin(ipts)-yin_found(ipts))*(yin(ipts)-yin_found(ipts))
        dz_dummy=(zin(ipts)-zin_found(ipts))*(zin(ipts)-zin_found(ipts))
        if(dist_dummy>100 .or. dx_dummy>100 .or. dy_dummy>100 .or. dz_dummy>100 )then
          !print*,'min dist,dx,dy,dz in rank',dist_dummy,dx_dummy,dy_dummy,dz_dummy
          !print*,'x,x_found,y,y_found,z,z_found in rank',xin(ipts),xin_found(ipts),yin(ipts),&
          !        yin_found(ipts),zin(ipts),zin_found(ipts)
          
        else
          do k = 1,NGLLZ
          do j = 1,NGLLY
          do i = 1,NGLLX
            hlagrange = hxi_store(ipts,i) * heta_store(ipts,j) * hgamma_store(ipts,k)
            ispec=ispec_selected_store(ipts)
            model_out_dummy(ipts) = model_out_dummy(ipts) + dble(model_in(i,j,k,ispec))*hlagrange
            enddo
            enddo
          enddo
        endif
        ! the points are find very accurately
        ! print*,'ipts,x,x_f,dist',ipts,xin(ipts),xin_found(ipts),dist_xyz_newxyz(ipts)
        ! print*,'ipts,y,y_f,dist',ipts,yin(ipts),yin_found(ipts),dist_xyz_newxyz(ipts)
        ! print*,'ipts,z,z_f,dist',ipts,zin(ipts),zin_found(ipts),dist_xyz_newxyz(ipts)
    endif
  enddo
  call MPI_Reduce(model_out_dummy,model_out,npts,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ier)
end subroutine interp

subroutine find_neareast()

  use mpi_f08
  implicit none
  integer  ::  ipts,i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dist_local,model_local,dist_global,model_global
  allocate(dist_local(npts),stat=ier)
  allocate(dist_global(npts),stat=ier)
  allocate(model_local(npts),stat=ier)
  allocate(model_global(npts),stat=ier)
  dist_local(:) = huge(xin(1))
  do ispec = 1,NSPEC_AB
    do k=1,NGLLZ 
      do j=1,NGLLY 
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dist_global(:) = sqrt((xin(1:npts)-xstore(iglob))**2 &
                    +(yin(1:npts)-ystore(iglob))**2 &
                    +(zin(1:npts)-zstore(iglob))**2)
          do ipts=1,npts 
            if (dist_global(ipts) < dist_local(ipts)) then 
              model_local(ipts) = model_in(i,j,k,ispec)
              dist_local(ipts) = dist_global(ipts)
            endif 
          enddo
        enddo
      enddo
    enddo
  enddo
  if(rank == 0)  then 
    print*,'finished finding local distmin'
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! find the parameters for the nearest point in global
  call min_allprocs(model_local,dist_local,model_global,dist_global,npts)
  model_out(:)=model_global(:)
  
  DEALLOCATE(model_local,dist_local,dist_global,model_global)

end subroutine find_neareast

subroutine point_min_allprocs(dist_slice,point_location)
  use mpi_f08
  implicit none
  double precision,intent(in) :: dist_slice

  ! local
  double precision  :: distin(2),point_location(2)

  ! get current rank 
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank,ier)
  
  ! find min distance in global
  distin(1)=dble(dist_slice)
  distin(2)=dble(rank)
  ! find minloc
  call MPI_ALLREDUCE(distin,point_location,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ier)

end subroutine point_min_allprocs

subroutine min_allprocs(vpmin_slice,distmin_slice,vpout,distout,npts)
  use mpi_f08
  implicit none
  integer,intent(in)   :: npts
  real(kind=CUSTOM_REAL),intent(in) :: vpmin_slice(npts),distmin_slice(npts)
  real(kind=CUSTOM_REAL),intent(inout) :: vpout(npts),distout(npts)

  ! local
  integer                 :: myrank,ierr,ipt
  real(kind=CUSTOM_REAL)  :: distin(2,npts),distall(2,npts),vmin(npts),dmin(npts)

  ! get current rank 
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank,ierr)

  ! find min distance in global
  distin(1,:) = distmin_slice(:)
  distin(2,:) = real(myrank,kind=CUSTOM_REAL)

  ! find minloc
  call MPI_ALLREDUCE(distin,distall,npts,MPI_2REAL,MPI_MINLOC,MPI_COMM_WORLD,ierr)
  vmin(:) = 0.0_CUSTOM_REAL
  dmin(:) = 0.0_CUSTOM_REAL

  do ipt=1,npts 
    if(myrank == int(distall(2,ipt))) then 
      vmin(ipt) = vpmin_slice(ipt)
      dmin(ipt) = distmin_slice(ipt)
    endif
  enddo
  
  call MPI_Reduce(vmin,vpout,npts,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(dmin,distout,npts,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

end subroutine min_allprocs
subroutine interp_finalize()
  !free arrays
  if(allocated(ibool))deallocate(ibool)
  if(allocated(xstore))deallocate(xstore)
  if(allocated(ystore))deallocate(ystore)
  if(allocated(zstore))deallocate(zstore)
  if(allocated(irregular_element_number))deallocate(irregular_element_number)
  if(allocated(xix))deallocate(xix)
  if(allocated(xiy))deallocate(xiy)
  if(allocated(xiz))deallocate(xiz)
  if(allocated(etax))deallocate(etax)
  if(allocated(etay))deallocate(etay)
  if(allocated(etaz))deallocate(etaz)
  if(allocated(gammax))deallocate(gammax)
  if(allocated(gammay))deallocate(gammay)
  if(allocated(gammaz))deallocate(gammaz)
  if(allocated(jacobian))deallocate(jacobian)
  if(allocated(kappastore))deallocate(kappastore)
  if(allocated(rhostore))deallocate(rhostore)
  if(allocated(mustore))deallocate(mustore)
  if(allocated(ispec_is_poroelastic))deallocate(ispec_is_poroelastic)
  if(allocated(ispec_is_elastic))deallocate(ispec_is_elastic)
  if(allocated(ispec_is_acoustic))deallocate(ispec_is_acoustic)
  if(allocated(hxi))deallocate(hxi)
  if(allocated(heta))deallocate(heta)
  if(allocated(hgamma))deallocate(hgamma)
  if(allocated(hpxi))deallocate(hpxi)
  if(allocated(hpeta))deallocate(hpeta)
  if(allocated(hpgamma))deallocate(hpgamma)
  if(allocated(hxi_store))deallocate(hxi_store)
  if(allocated(heta_store))deallocate(heta_store)
  if(allocated(hgamma_store))deallocate(hgamma_store)
  if(allocated(points_location_store))deallocate(points_location_store)
  if(allocated(ispec_selected_store))deallocate(ispec_selected_store)
  if(allocated(model_in))deallocate(model_in)
  if(allocated(model_out))deallocate(model_out)
  if(allocated(xin))deallocate(xin)
  if(allocated(yin))deallocate(yin)
  if(allocated(zin))deallocate(zin)
  if(allocated(xin_found))deallocate(xin_found)
  if(allocated(yin_found))deallocate(yin_found)
  if(allocated(zin_found))deallocate(zin_found)
  if(allocated(x_xi))deallocate(x_xi)
  if(allocated(y_eta))deallocate(y_eta)
  if(allocated(z_gamma))deallocate(z_gamma)
  if(allocated(dist_xyz_newxyz))deallocate(dist_xyz_newxyz)
  



end subroutine interp_finalize
end module interp3D
