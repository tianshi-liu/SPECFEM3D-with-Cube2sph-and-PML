subroutine adepml_set_index_glob(myrank,xstore,ystore,zstore)
  use generate_databases_par, only: ibool,NSPEC_AB,NGLOB_AB,NDIM,CUSTOM_REAL,&
    NGLLX,NGLLY,NGLLZ,CPML_to_spec,PML_INSTEAD_OF_FREE_SURFACE,&
    ADEPML_DEFORMED,MAX_STRING_LEN,LOCAL_PATH,&
    CPML_width_x,CPML_width_y,CPML_width_z,nspec_cpml,&
    CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,&
    CPML_YZ_ONLY,CPML_XYZ,ZERO,ONE,TWO,HUGEVAL,&
    IMAIN,SMALLVAL,MIDX,MIDY,MIDZ,NGLLSQUARE,CPML_regions,&
    ibool_CPML,nglob_CPML,CPML_to_glob,&
    num_interfaces_PML,max_nibool_interfaces_PML,&
    num_interfaces_ext_mesh,nibool_interfaces_ext_mesh,&
    ibool_interfaces_ext_mesh,ibool_interfaces_PML,my_neighbors_PML,&
    nibool_interfaces_PML,my_neighbors_ext_mesh, &
    nglob_pml_in, pml_in_iglob, nglob_dirichlet, iglob_dirichlet
  implicit none
  integer, intent(in) :: myrank
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: xstore,ystore,zstore
  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: &
    xstore_undeformed,ystore_undeformed,zstore_undeformed
  double precision,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: x,y,z
  integer :: i,j,k,ispec,iglob,iglob_CPML,ispec_CPML,ier,iface,igll
  real(kind=CUSTOM_REAL) :: xoriginleft,xoriginright,yoriginfront,yoriginback,zoriginbottom,zorigintop
  real(kind=CUSTOM_REAL) :: x_min,x_min_all,y_min,y_min_all,z_min,z_min_all, &
                            x_max,x_max_all,y_max,y_max_all,z_max,z_max_all, &
                            x_origin,y_origin,z_origin,xtol,ytol,ztol,&
                            xc,yc,zc
  real(kind=CUSTOM_REAL) :: CPML_width_x_left, CPML_width_x_right, &
                            CPML_width_y_front,CPML_width_y_back, &
                            CPML_width_z_top,CPML_width_z_bottom, &
                            CPML_x_left, CPML_x_right, &
                            CPML_y_front,CPML_y_back, &
                            CPML_z_top,CPML_z_bottom, &
                       CPML_width_x_left_max_all, CPML_width_x_right_max_all, &
                       CPML_width_y_front_max_all,CPML_width_y_back_max_all,&
                       CPML_width_z_top_max_all,CPML_width_z_bottom_max_all
  integer, allocatable :: iglob_CPML_temp(:),nibool_interfaces_PML_rough(:)
  logical :: POINT_EXIST
  character(len=MAX_STRING_LEN) :: infn, procname

  write(procname,"('/proc',i6.6,'_')") myrank

  infn = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(procname)//'undeformed_xyz.bin'
  
  if (ADEPML_DEFORMED) then
  !! read in the undeformed coordinates
    open(unit=99,file=trim(infn),status='unknown',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) print *, 'file' ,trim(infn), 'not found'
    read(99) x
    read(99) y
    read(99) z
    do ispec = 1, NSPEC_AB
      do i=1,NGLLX;do j=1,NGLLY;do k=1,NGLLZ
        iglob = ibool(i,j,k,ispec)
        xstore_undeformed(iglob) = real(x(i,j,k,ispec),kind=CUSTOM_REAL)
        ystore_undeformed(iglob) = real(y(i,j,k,ispec),kind=CUSTOM_REAL)
        zstore_undeformed(iglob) = real(z(i,j,k,ispec),kind=CUSTOM_REAL)
      enddo;enddo;enddo
    enddo
  else
    xstore_undeformed(:) = xstore(:)
    ystore_undeformed(:) = ystore(:)
    zstore_undeformed(:) = zstore(:)
  endif
! Assuming the computational domain is convex and can be approximatly seen as a
! box
! Calculation of origin of whole computational domain
  x_min = minval(xstore_undeformed(:))
  x_max = maxval(xstore_undeformed(:))
  y_min = minval(ystore_undeformed(:))
  y_max = maxval(ystore_undeformed(:))
  z_min = minval(zstore_undeformed(:))
  z_max = maxval(zstore_undeformed(:))

  x_min_all = HUGEVAL
  y_min_all = HUGEVAL
  z_min_all = HUGEVAL

  x_max_all = - HUGEVAL
  y_max_all = - HUGEVAL
  z_max_all = - HUGEVAL

  call min_all_all_cr(x_min,x_min_all)
  call min_all_all_cr(y_min,y_min_all)
  call min_all_all_cr(z_min,z_min_all)

  call max_all_all_cr(x_max,x_max_all)
  call max_all_all_cr(y_max,y_max_all)
  call max_all_all_cr(z_max,z_max_all)

  x_origin = (x_min_all + x_max_all) / TWO
  y_origin = (y_min_all + y_max_all) / TWO
  z_origin = (z_max_all + z_min_all) / TWO

! Assuming CPML_width_x,CPML_width_y,CPML_width_Z are constants inside PML layer
! Calculation of width of PML along x, y and z direction, such as
! CPML_width_x,CPML_width_y,CPML_width_Z
  CPML_width_x_left   = ZERO
  CPML_width_x_right  = ZERO
  CPML_width_y_front  = ZERO
  CPML_width_y_back   = ZERO
  CPML_width_z_top    = ZERO
  CPML_width_z_bottom = ZERO

  CPML_x_right  = x_max_all
  CPML_x_left   = x_min_all
  CPML_y_front  = y_max_all
  CPML_y_back   = y_min_all
  CPML_z_top    = z_max_all
  CPML_z_bottom = z_min_all

  do ispec_CPML = 1,nspec_cpml
    ispec = CPML_to_spec(ispec_CPML)

    do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX

      iglob = ibool(i,j,k,ispec)

      if (CPML_regions(ispec_CPML) == CPML_X_ONLY .or. CPML_regions(ispec_CPML) == CPML_XY_ONLY .or. &
          CPML_regions(ispec_CPML) == CPML_XZ_ONLY .or. CPML_regions(ispec_CPML) == CPML_XYZ) then
        if (xstore_undeformed(iglob) - x_origin > ZERO) then
          if (xstore_undeformed(iglob) - x_origin <= CPML_x_right - x_origin) then
            CPML_x_right = xstore_undeformed(iglob)
          endif
        else
          if (abs(xstore_undeformed(iglob) - x_origin) <= abs(CPML_x_left-x_origin)) then
            CPML_x_left = xstore_undeformed(iglob)
          endif
        endif
      endif

      if (CPML_regions(ispec_CPML) == CPML_Y_ONLY .or. CPML_regions(ispec_CPML) == CPML_XY_ONLY .or. &
          CPML_regions(ispec_CPML) == CPML_YZ_ONLY .or. CPML_regions(ispec_CPML) == CPML_XYZ) then
        if (ystore_undeformed(iglob) - y_origin > ZERO) then
          if (ystore_undeformed(iglob) - y_origin <= CPML_y_front - y_origin) then
            CPML_y_front = ystore_undeformed(iglob)
          endif
        else
          if (abs(ystore_undeformed(iglob) - y_origin) <= abs(CPML_y_back-y_origin)) then
            CPML_y_back = ystore_undeformed(iglob)
          endif
        endif
      endif

      if (CPML_regions(ispec_CPML) == CPML_Z_ONLY .or. CPML_regions(ispec_CPML) == CPML_YZ_ONLY .or. &
        CPML_regions(ispec_CPML) == CPML_XZ_ONLY .or. CPML_regions(ispec_CPML) == CPML_XYZ) then
        if (zstore_undeformed(iglob) - z_origin > ZERO) then
          if (zstore_undeformed(iglob) - z_origin <= CPML_z_top - z_origin) then
            CPML_z_top = zstore_undeformed(iglob)
          endif
        else
          if (abs(zstore_undeformed(iglob) - z_origin) <= abs(CPML_z_bottom-z_origin)) then
            CPML_z_bottom = zstore_undeformed(iglob)
          endif
        endif
      endif

    enddo; enddo; enddo
  enddo ! ispec_CPML = 1,nspec_cpml
    
  CPML_width_x_right  = x_max_all - CPML_x_right
  CPML_width_x_left   = CPML_x_left - x_min_all
  CPML_width_y_front  = y_max_all - CPML_y_front
  CPML_width_y_back   = CPML_y_back - y_min_all
  CPML_width_z_top    = z_max_all - CPML_z_top
  CPML_width_z_bottom = CPML_z_bottom - z_min_all

  call max_all_all_cr(CPML_width_x_left,CPML_width_x_left_max_all)
  call max_all_all_cr(CPML_width_x_right,CPML_width_x_right_max_all)
  call max_all_all_cr(CPML_width_y_front,CPML_width_y_front_max_all)
  call max_all_all_cr(CPML_width_y_back,CPML_width_y_back_max_all)
  call max_all_all_cr(CPML_width_z_top,CPML_width_z_top_max_all)
  call max_all_all_cr(CPML_width_z_bottom,CPML_width_z_bottom_max_all)

  xoriginleft   = x_min_all + CPML_width_x_left_max_all
  xoriginright  = x_max_all - CPML_width_x_right_max_all
  yoriginback   = y_min_all + CPML_width_y_back_max_all
  yoriginfront  = y_max_all - CPML_width_y_front_max_all
  zoriginbottom = z_min_all + CPML_width_z_bottom_max_all

  CPML_width_x = max(CPML_width_x_left_max_all,CPML_width_x_right_max_all)
  CPML_width_y = max(CPML_width_y_front_max_all,CPML_width_y_back_max_all)
  CPML_width_z = max(CPML_width_z_bottom_max_all,CPML_width_z_top_max_all)
  ! tolerations for geometrical detection
  xtol = CPML_width_x*SMALLVAL
  ytol = CPML_width_y*SMALLVAL
  ztol = CPML_width_z*SMALLVAL

  if (myrank == 0) then
    print *, 'xmin: ', x_min_all, 'xmax:', x_max_all
    print *, 'ymin: ', y_min_all, 'ymax:', y_max_all
    print *, 'zmin: ', z_min_all, 'zmax:', z_max_all
    print *, 'CPML_width_x: ', CPML_width_x
    print *, 'CPML_width_y: ', CPML_width_y
    print *, 'CPML_width_z: ', CPML_width_z
  endif

  if (PML_INSTEAD_OF_FREE_SURFACE) then
    zorigintop = z_max_all - CPML_width_z_top_max_all
  else
    zorigintop = z_max_all
  endif

  call synchronize_all()

  nglob_dirichlet = 0 
  do iglob = 1, NGLOB_AB
    xc = xstore_undeformed(iglob)
    yc = ystore_undeformed(iglob)
    zc = zstore_undeformed(iglob)
    if (is_dirichlet(xc,yc,zc)) nglob_dirichlet =nglob_dirichlet + 1
  enddo 
  if (nglob_dirichlet > 0) then
    allocate(iglob_dirichlet(nglob_dirichlet))
    nglob_dirichlet = 0
    do iglob = 1, NGLOB_AB
      xc = xstore_undeformed(iglob)
      yc = ystore_undeformed(iglob)
      zc = zstore_undeformed(iglob)
      if (is_dirichlet(xc,yc,zc)) then
        nglob_dirichlet =nglob_dirichlet + 1
        iglob_dirichlet(nglob_dirichlet) = iglob
      endif
    enddo
  endif

  !! indexing in the CPML domain
  allocate(ibool_CPML(NGLLX,NGLLY,NGLLZ,nspec_cpml))
  allocate(iglob_CPML_temp(NGLLX*NGLLY*NGLLZ*nspec_cpml))
  nglob_CPML = 0
  ibool_CPML(:,:,:,:) = 0
  do ispec_CPML = 1, nspec_cpml
    do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
      ispec = CPML_to_spec(ispec_CPML)
      iglob = ibool(i,j,k,ispec)
      POINT_EXIST = .false.
      do iglob_CPML = 1, nglob_CPML
        if (iglob == iglob_CPML_temp(iglob_CPML)) then
          POINT_EXIST = .true.
          exit
        endif
      enddo
      if (.not. POINT_EXIST) then
        nglob_CPML = nglob_CPML + 1
        iglob_CPML_temp(nglob_CPML) = iglob
        ibool_CPML(i,j,k,ispec_CPML) = nglob_CPML
      else
        ibool_CPML(i,j,k,ispec_CPML) = iglob_CPML  
      endif
    enddo;enddo;enddo
  enddo
  allocate(CPML_to_glob(nglob_CPML))
  CPML_to_glob(1:nglob_CPML) = iglob_CPML_temp(1:nglob_CPML)
  !deallocate(iglob_CPML_temp)  

  !! indexing at the interface between PML and physical domain
  ! a count first
  nglob_pml_in = 0
  do iglob_CPML = 1, nglob_CPML
    iglob = CPML_to_glob(iglob_CPML)
    xc = xstore_undeformed(iglob)
    yc = ystore_undeformed(iglob)
    zc = zstore_undeformed(iglob)
    if (.not. is_point_inside_PML(xc,yc,zc)) nglob_pml_in = nglob_pml_in + 1
  enddo
  if (nglob_pml_in > 0) then
    allocate(pml_in_iglob(nglob_pml_in))
    nglob_pml_in = 0
    do iglob_CPML = 1, nglob_CPML
      iglob = CPML_to_glob(iglob_CPML)
      xc = xstore_undeformed(iglob)
      yc = ystore_undeformed(iglob)
      zc = zstore_undeformed(iglob)
      if (.not. is_point_inside_PML(xc,yc,zc)) then
        nglob_pml_in = nglob_pml_in + 1
        pml_in_iglob(nglob_pml_in) = iglob_CPML
      endif
    enddo
  endif

  !! MPI interface
  num_interfaces_PML = 0
  max_nibool_interfaces_PML = 0
  if (num_interfaces_ext_mesh > 0) then
    !! preliminary counting
    allocate(nibool_interfaces_PML_rough(num_interfaces_ext_mesh))
    nibool_interfaces_PML_rough(:) = 0
    do iface = 1, num_interfaces_ext_mesh
      do igll = 1, nibool_interfaces_ext_mesh(iface)
        iglob = ibool_interfaces_ext_mesh(igll,iface)
        xc = xstore_undeformed(iglob)
        yc = ystore_undeformed(iglob)
        zc = zstore_undeformed(iglob)
        if (is_point_inside_PML(xc,yc,zc)) &
          nibool_interfaces_PML_rough(iface) = &
              nibool_interfaces_PML_rough(iface)+1
      enddo
      if (nibool_interfaces_PML_rough(iface) > 0) &
        num_interfaces_PML = num_interfaces_PML + 1
    enddo
    max_nibool_interfaces_PML = maxval(nibool_interfaces_PML_rough(:))
    if (max_nibool_interfaces_PML > 0) then
      !! we actually have PML interfaces, start actual indexing
      allocate(nibool_interfaces_PML(num_interfaces_PML))
      allocate(ibool_interfaces_PML(max_nibool_interfaces_PML,&
                                               num_interfaces_PML))
      allocate(my_neighbors_PML(num_interfaces_PML))
      nibool_interfaces_PML(:) = 0
      num_interfaces_PML = 0
      do iface = 1, num_interfaces_ext_mesh
        if (nibool_interfaces_PML_rough(iface) .eq. 0) cycle
        num_interfaces_PML = num_interfaces_PML + 1
        do igll = 1, nibool_interfaces_ext_mesh(iface)
          iglob = ibool_interfaces_ext_mesh(igll,iface)
          xc = xstore_undeformed(iglob)
          yc = ystore_undeformed(iglob)
          zc = zstore_undeformed(iglob)
          if (is_point_inside_PML(xc,yc,zc)) then
            nibool_interfaces_PML(num_interfaces_PML) = &
              nibool_interfaces_PML(num_interfaces_PML)+1
            POINT_EXIST = .false.
            do iglob_CPML = 1, nglob_CPML
              if (iglob == iglob_CPML_temp(iglob_CPML)) then
                ibool_interfaces_PML(nibool_interfaces_PML(num_interfaces_PML)&
                      ,num_interfaces_PML)=iglob_CPML
                POINT_EXIST = .true.
                exit
              endif
            enddo
            if (.not. POINT_EXIST) print *, 'PML GLL wrong indexing'
          endif
        enddo
        my_neighbors_PML(num_interfaces_PML) = my_neighbors_ext_mesh(iface)
        if (.not. (nibool_interfaces_PML(num_interfaces_PML) .eq. &
                nibool_interfaces_PML_rough(iface))) &
            print *, 'indexing not match for PML interface'
      enddo
    endif
    deallocate(nibool_interfaces_PML_rough)
  endif
  deallocate(iglob_CPML_temp)  

  contains
    logical function is_point_inside_PML(xp,yp,zp)
      real,intent(in) :: xp,yp,zp
      is_point_inside_PML = .false.
      if ((xp<(xoriginleft-xtol)).or.(xp>(xoriginright+xtol))&
        .or.(yp<(yoriginback-ytol)).or.(yp>(yoriginfront+ytol))&
        .or.(zp<(zoriginbottom-ztol)).or. &
            ((zp>(zorigintop+ztol)).or.PML_INSTEAD_OF_FREE_SURFACE)) then
        !if ((xp>(x_min_all+xtol)).and.(xp<(x_max_all-xtol)) &
        !    .and.(yp>(y_min_all+ytol)).and.(yp<(y_max_all-ytol)) &
        !    .and.(zp>(z_min_all+ztol)).and.(zp<(z_max_all-ztol))) &
        is_point_inside_PML = .true.
      endif
    end function is_point_inside_PML

    logical function is_dirichlet(xp,yp,zp)
      real,intent(in) :: xp,yp,zp
      is_dirichlet = .false.
      if((xp<(x_min_all+xtol)).or.(xp>(x_max_all-xtol))&
        .or.(yp<(y_min_all+ytol)).or.(yp>(y_max_all-ytol))&
        .or.(zp<(z_min_all+ztol)).or. &
            ((zp>(z_max_all-ztol)).and.PML_INSTEAD_OF_FREE_SURFACE)) then
        is_dirichlet = .true.
      endif
    end function is_dirichlet

end subroutine adepml_set_index_glob
