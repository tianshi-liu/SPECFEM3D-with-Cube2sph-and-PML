subroutine adepml_set_local_dampingcoeff(myrank,xstore,ystore,zstore)
  !! Note that this subroutine must be called before cube2sph transform
  use generate_databases_par, only: ibool,NGLOB_AB,NDIM,CUSTOM_REAL,ZERO,&
          NGLLX,NGLLY,NGLLZ,pml_d,pml_kappa,pml_beta,NPOWER,PML_P,&
          PML_KAPPA0,MPML_P,ONE,PI,f0_FOR_PML,CPML_to_spec,&
          CPML_width_x,CPML_width_y,CPML_width_z,nspec_cpml,&
          PML_INSTEAD_OF_FREE_SURFACE,&
          CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,&
          CPML_YZ_ONLY,CPML_XYZ,is_CPML,ONE,TWO,HUGEVAL,&
          IMAIN,SMALLVAL,MIDX,MIDY,MIDZ,NGLLSQUARE,&
          CPML_regions,num_pml_out,pml_out_ispec,pml_out_ijk
          !ibool_CPML,&
          !nglob_CPML,CPML_to_glob,CPML_regions
          !num_interfaces_PML,max_nibool_interfaces_PML,&
          !num_interfaces_ext_mesh,nibool_interfaces_ext_mesh,&
          !ibool_interfaces_ext_mesh,ibool_interfaces_PML,my_neighbors_PML,&
          !nibool_interfaces_PML,my_neighbors_ext_mesh
  use create_regions_mesh_ext_par, only: rhostore,rho_vp
  implicit none
  integer, intent(in) :: myrank
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: xstore,ystore,zstore
  ! local parameters
  integer :: i,j,k,ispec,iglob,iglobc,ispec_CPML,ier
  real(kind=CUSTOM_REAL) :: pml_damping_profile_l,dist,vp,d_temp
  real(kind=CUSTOM_REAL) :: xoriginleft,xoriginright,yoriginfront,yoriginback,zoriginbottom,zorigintop
  real(kind=CUSTOM_REAL) :: abscissa_in_PML_x,abscissa_in_PML_y,abscissa_in_PML_z
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
                            CPML_width_y_front_max_all,CPML_width_y_back_max_all, &
                            CPML_width_z_top_max_all,CPML_width_z_bottom_max_all, &
                            vp_max,vp_max_all
  !integer, allocatable :: iglob_CPML_temp(:)!,nibool_interfaces_PML_rough(:)
  !logical :: POINT_EXIST
! for robust parameter separation of PML damping parameter

  ! checks number of PML elements
  if (count(is_CPML(:)) /= NSPEC_CPML) then
    print *,'Error in slice ',myrank,': number of PML elements ',NSPEC_CPML,' but only ',count(is_CPML(:)),' flags set'
    stop 'Error C-PML array has invalid number of PML flags set'
  endif

  ! checks if C-PML flags assigned correctly
  do ispec_CPML = 1,NSPEC_CPML
    ispec = CPML_to_spec(ispec_CPML)
    if (.not. is_CPML(ispec)) stop 'Error found C-PML element with invalid PML flag'
  enddo
  ! stores damping profiles & auxiliary coefficients
  allocate(pml_d(NDIM,NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 841')
  if (ier /= 0) stop 'error allocating array pml_d'
  allocate(pml_kappa(NDIM,NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 844')
  if (ier /= 0) stop 'error allocating array pml_kappa'
  allocate(pml_beta(NDIM,NGLLX,NGLLY,NGLLZ,nspec_cpml),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 847')
  if (ier /= 0) stop 'error allocating array pml_beta'

  pml_d(:,:,:,:,:) = ZERO
  pml_kappa(:,:,:,:,:) = ONE
  pml_beta(:,:,:,:,:) = PI*f0_FOR_PML

! Assuming the computational domain is convex and can be approximatly seen as a box
! Calculation of origin of whole computational domain
  x_min = minval(xstore(:))
  x_max = maxval(xstore(:))
  y_min = minval(ystore(:))
  y_max = maxval(ystore(:))
  z_min = minval(zstore(:))
  z_max = maxval(zstore(:))

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
! Calculation of width of PML along x, y and z direction, such as CPML_width_x,CPML_width_y,CPML_width_Z
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
        if (xstore(iglob) - x_origin > ZERO) then
          if (xstore(iglob) - x_origin <= CPML_x_right - x_origin) then
            CPML_x_right = xstore(iglob)
          endif
        else
          if (abs(xstore(iglob) - x_origin) <= abs(CPML_x_left-x_origin)) then
            CPML_x_left = xstore(iglob)
          endif
        endif
      endif

      if (CPML_regions(ispec_CPML) == CPML_Y_ONLY .or. CPML_regions(ispec_CPML) == CPML_XY_ONLY .or. &
          CPML_regions(ispec_CPML) == CPML_YZ_ONLY .or. CPML_regions(ispec_CPML) == CPML_XYZ) then
        if (ystore(iglob) - y_origin > ZERO) then
          if (ystore(iglob) - y_origin <= CPML_y_front - y_origin) then
            CPML_y_front = ystore(iglob)
          endif
        else
          if (abs(ystore(iglob) - y_origin) <= abs(CPML_y_back-y_origin)) then
            CPML_y_back = ystore(iglob)
          endif
        endif
      endif

      if (CPML_regions(ispec_CPML) == CPML_Z_ONLY .or. CPML_regions(ispec_CPML) == CPML_YZ_ONLY .or. &
        CPML_regions(ispec_CPML) == CPML_XZ_ONLY .or. CPML_regions(ispec_CPML) == CPML_XYZ) then
        if (zstore(iglob) - z_origin > ZERO) then
          if (zstore(iglob) - z_origin <= CPML_z_top - z_origin) then
            CPML_z_top = zstore(iglob)
          endif
        else
          if (abs(zstore(iglob) - z_origin) <= abs(CPML_z_bottom-z_origin)) then
            CPML_z_bottom = zstore(iglob)
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

  if (PML_INSTEAD_OF_FREE_SURFACE) then
    zorigintop = z_max_all - CPML_width_z_top_max_all
  else
    zorigintop = z_max_all
  endif

! Calculation of maximum p velocity inside PML
  vp_max = ZERO
  do ispec_CPML = 1,nspec_cpml
    ispec = CPML_to_spec(ispec_CPML)
    do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
      vp = rho_vp(i,j,k,ispec) / rhostore(i,j,k,ispec)
      if (vp >= vp_max) then
        vp_max = vp
      endif
    enddo; enddo; enddo
  enddo

  call max_all_all_cr(vp_max,vp_max_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Boundary values of X-/Y-/Z-regions'
    write(IMAIN,*) minval(xstore(:)), maxval(xstore(:))
    write(IMAIN,*) minval(ystore(:)), maxval(ystore(:))
    write(IMAIN,*) minval(zstore(:)), maxval(zstore(:))
    write(IMAIN,*)
    write(IMAIN,*) 'Origins of right/left X-surface C-PML',xoriginright,xoriginleft
    write(IMAIN,*) 'Origins of front/back Y-surface C-PML',yoriginfront,yoriginback
    write(IMAIN,*) 'Origin of bottom Z-surface C-PML',zoriginbottom
    if (PML_INSTEAD_OF_FREE_SURFACE) then
      write(IMAIN,*) 'Origin of top Z-surface C-PML',zorigintop
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'CPML_width_x: ',CPML_width_x
    write(IMAIN,*) 'CPML_width_y: ',CPML_width_y
    write(IMAIN,*) 'CPML_width_z: ',CPML_width_z
    write(IMAIN,*)
  endif

  call synchronize_all()

  !call wait_for_attach(myrank)

  ! loops over all C-PML elements
  do ispec_CPML = 1,nspec_cpml

    ispec = CPML_to_spec(ispec_CPML)
    vp = vp_max_all

    iglobc = ibool(MIDX,MIDY,MIDZ,ispec)
    !! right PML
    if (xstore(iglobc) - xoriginright > xtol) then

      do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)
        ! gets abscissa of current grid point along the damping profile
        abscissa_in_PML_x = xstore(iglob) - xoriginright

        ! determines distance to C-PML/mesh interface
        dist = abscissa_in_PML_x / CPML_width_x
          
        d_temp = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
        pml_d(1,i,j,k,ispec_CPML)=pml_d(1,i,j,k,ispec_CPML)+d_temp
        pml_d(2,i,j,k,ispec_CPML)=pml_d(2,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_d(3,i,j,k,ispec_CPML)=pml_d(3,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_beta(1,i,j,k,ispec_CPML)=pml_beta(1,i,j,k,ispec_CPML)*&
                (ONE-(dist**PML_P))
        pml_kappa(1,i,j,k,ispec_CPML)=ONE+(PML_KAPPA0-ONE)*(dist**NPOWER)
      enddo;enddo;enddo
    endif

    !! left PML
    if (xstore(iglobc) - xoriginleft < -xtol) then

      do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)
        ! gets abscissa of current grid point along the damping profile
        abscissa_in_PML_x = xoriginleft - xstore(iglob)

        ! determines distance to C-PML/mesh interface
        dist = abscissa_in_PML_x / CPML_width_x

        d_temp = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_x)
        pml_d(1,i,j,k,ispec_CPML)=pml_d(1,i,j,k,ispec_CPML)+d_temp
        pml_d(2,i,j,k,ispec_CPML)=pml_d(2,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_d(3,i,j,k,ispec_CPML)=pml_d(3,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_beta(1,i,j,k,ispec_CPML)=pml_beta(1,i,j,k,ispec_CPML)*&
                (ONE-(dist**PML_P))
        pml_kappa(1,i,j,k,ispec_CPML)=ONE+(PML_KAPPA0-ONE)*(dist**NPOWER)
      enddo;enddo;enddo
    endif

    !! front PML
    if (ystore(iglobc) - yoriginfront > ytol) then

      do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)
        ! gets abscissa of current grid point along the damping profile
        abscissa_in_PML_y = ystore(iglob) - yoriginfront

        ! determines distance to C-PML/mesh interface
        dist = abscissa_in_PML_y / CPML_width_y

        d_temp = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
        pml_d(2,i,j,k,ispec_CPML)=pml_d(2,i,j,k,ispec_CPML)+d_temp
        pml_d(1,i,j,k,ispec_CPML)=pml_d(1,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_d(3,i,j,k,ispec_CPML)=pml_d(3,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_beta(2,i,j,k,ispec_CPML)=pml_beta(2,i,j,k,ispec_CPML)*&
                (ONE-(dist**PML_P))
        pml_kappa(2,i,j,k,ispec_CPML)=ONE+(PML_KAPPA0-ONE)*(dist**NPOWER)
      enddo;enddo;enddo
    endif

    !! back PML
    if (ystore(iglobc) - yoriginback < -ytol) then

      do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)
        ! gets abscissa of current grid point along the damping profile
        abscissa_in_PML_y = yoriginback - ystore(iglob)

        ! determines distance to C-PML/mesh interface
        dist = abscissa_in_PML_y / CPML_width_y

        d_temp = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_y)
        pml_d(2,i,j,k,ispec_CPML)=pml_d(2,i,j,k,ispec_CPML)+d_temp
        pml_d(1,i,j,k,ispec_CPML)=pml_d(1,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_d(3,i,j,k,ispec_CPML)=pml_d(3,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_beta(2,i,j,k,ispec_CPML)=pml_beta(2,i,j,k,ispec_CPML)*&
                (ONE-(dist**PML_P))
        pml_kappa(2,i,j,k,ispec_CPML)=ONE+(PML_KAPPA0-ONE)*(dist**NPOWER)
      enddo;enddo;enddo
    endif

    !! top PML
    if ((zstore(iglobc)-zorigintop>ztol).and.PML_INSTEAD_OF_FREE_SURFACE) then

      do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)
        ! gets abscissa of current grid point along the damping profile
        abscissa_in_PML_z = zstore(iglob) - zorigintop

        ! determines distance to C-PML/mesh interface
        dist = abscissa_in_PML_z / CPML_width_z

        d_temp = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
        pml_d(3,i,j,k,ispec_CPML)=pml_d(3,i,j,k,ispec_CPML)+d_temp
        pml_d(1,i,j,k,ispec_CPML)=pml_d(1,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_d(2,i,j,k,ispec_CPML)=pml_d(2,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_beta(3,i,j,k,ispec_CPML)=pml_beta(3,i,j,k,ispec_CPML)*&
                (ONE-(dist**PML_P))
        pml_kappa(3,i,j,k,ispec_CPML)=ONE+(PML_KAPPA0-ONE)*(dist**NPOWER)
      enddo;enddo;enddo
    endif

    !! bottom PML
    if (zstore(iglobc) - zoriginbottom < -ztol) then

      do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)
        ! gets abscissa of current grid point along the damping profile
        abscissa_in_PML_z = zoriginbottom - zstore(iglob)

        ! determines distance to C-PML/mesh interface
        dist = abscissa_in_PML_z / CPML_width_z

        d_temp = pml_damping_profile_l(myrank,iglob,dist,vp,CPML_width_z)
        pml_d(3,i,j,k,ispec_CPML)=pml_d(3,i,j,k,ispec_CPML)+d_temp
        pml_d(1,i,j,k,ispec_CPML)=pml_d(1,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_d(2,i,j,k,ispec_CPML)=pml_d(2,i,j,k,ispec_CPML)+&
                MPML_P*(dist**NPOWER)*d_temp
        pml_beta(3,i,j,k,ispec_CPML)=pml_beta(3,i,j,k,ispec_CPML)*&
                (ONE-(dist**PML_P))
        pml_kappa(3,i,j,k,ispec_CPML)=ONE+(PML_KAPPA0-ONE)*(dist**NPOWER)
      enddo;enddo;enddo
    endif



  enddo  !ispec_CPML = 1,nspec_cpml
  pml_d(:,:,:,:,:) = pml_d(:,:,:,:,:)&
          /(pml_kappa(:,:,:,:,:)**2)
  pml_beta(:,:,:,:,:) = pml_beta(:,:,:,:,:)&
          +pml_d(:,:,:,:,:)*pml_kappa(:,:,:,:,:)

  !!! identifying PML outer interface and indexing
  ! do a counting first
  num_pml_out = 0
  do ispec_CPML = 1,nspec_cpml
    ispec = CPML_to_spec(ispec_CPML)
    iglobc = ibool(MIDX,MIDY,1,ispec)
    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
    if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
    iglobc = ibool(MIDX,MIDY,NGLLZ,ispec)
    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
    if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
    iglobc = ibool(MIDX,1,MIDZ,ispec)
    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
    if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
    iglobc = ibool(MIDX,NGLLY,MIDZ,ispec)
    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
    if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
    iglobc = ibool(1,MIDY,MIDZ,ispec)
    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
    if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
    iglobc = ibool(NGLLX,MIDY,MIDZ,ispec)
    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
    if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
  enddo !ispec_CPML = 1,nspec_cpml
  if (num_pml_out > 0) then
    allocate(pml_out_ispec(num_pml_out))
    allocate(pml_out_ijk(NDIM, NGLLSQUARE, num_pml_out))
    ! actual indexing
    num_pml_out = 0
    do ispec_CPML = 1,nspec_cpml
      ispec = CPML_to_spec(ispec_CPML)
      iglobc = ibool(MIDX,MIDY,1,ispec)
      xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
      if (is_interface_out(xc,yc,zc)) call include_this_face_out(3,1)
      iglobc = ibool(MIDX,MIDY,NGLLZ,ispec)
      xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
      if (is_interface_out(xc,yc,zc)) call include_this_face_out(3,NGLLZ)
      iglobc = ibool(MIDX,1,MIDZ,ispec)
      xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
      if (is_interface_out(xc,yc,zc)) call include_this_face_out(2,1)
      iglobc = ibool(MIDX,NGLLY,MIDZ,ispec)
      xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
      if (is_interface_out(xc,yc,zc)) call include_this_face_out(2,NGLLY)
      iglobc = ibool(1,MIDY,MIDZ,ispec)
      xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
      if (is_interface_out(xc,yc,zc)) call include_this_face_out(1,1)
      iglobc = ibool(NGLLX,MIDY,MIDZ,ispec)
      xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
      if (is_interface_out(xc,yc,zc)) call include_this_face_out(1,NGLLX)
    enddo !ispec_CPML = 1,nspec_cpml
  endif
  
  !!! identifying PML inner interface and indexing
  ! do a counting first
  !num_pml_in = 0
  !do ispec_CPML = 1,nspec_cpml
  !  ispec = CPML_to_spec(ispec_CPML)
  !  iglobc = ibool(MIDX,MIDY,1,ispec)
  !  xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !  if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
  !  iglobc = ibool(MIDX,MIDY,NGLLZ,ispec)
  !  xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !  if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
  !  iglobc = ibool(MIDX,1,MIDZ,ispec)
  !  xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !  if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
  !  iglobc = ibool(MIDX,NGLLY,MIDZ,ispec)
  !  xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !  if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
  !  iglobc = ibool(1,MIDY,MIDZ,ispec)
  !  xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !  if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
  !  iglobc = ibool(NGLLX,MIDY,MIDZ,ispec)
  !  xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !  if (is_interface_out(xc,yc,zc)) num_pml_out = num_pml_out + 1
  !enddo !ispec_CPML = 1,nspec_cpml
  !if (num_pml_out > 0) then
  !  allocate(pml_out_ispec(num_pml_out))
  !  allocate(pml_out_ijk(NDIM, NGLLSQUARE, num_pml_out))
  !  ! actual indexing
  !  num_pml_out = 0
  !  do ispec_CPML = 1,nspec_cpml
  !    ispec = CPML_to_spec(ispec_CPML)
  !    iglobc = ibool(MIDX,MIDY,1,ispec)
  !    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !    if (is_interface_in(xc,yc,zc)) call include_this_face_in(3,1)
  !    iglobc = ibool(MIDX,MIDY,NGLLZ,ispec)
  !    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !    if (is_interface_in(xc,yc,zc)) call include_this_face_in(3,NGLLZ)
  !    iglobc = ibool(MIDX,1,MIDZ,ispec)
  !    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !    if (is_interface_in(xc,yc,zc)) call include_this_face_in(2,1)
  !    iglobc = ibool(MIDX,NGLLY,MIDZ,ispec)
  !    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !    if (is_interface_in(xc,yc,zc)) call include_this_face_in(2,NGLLY)
  !    iglobc = ibool(1,MIDY,MIDZ,ispec)
  !    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !    if (is_interface_in(xc,yc,zc)) call include_this_face_in(1,1)
  !    iglobc = ibool(NGLLX,MIDY,MIDZ,ispec)
  !    xc = xstore(iglobc); yc = ystore(iglobc); zc = zstore(iglobc)
  !    if (is_interface_in(xc,yc,zc)) call include_this_face_in(1,NGLLX)
  !  enddo !ispec_CPML = 1,nspec_cpml
  !endif
  
  !!! indexing in the CPML domain
  !allocate(ibool_CPML(NGLLX,NGLLY,NGLLZ,nspec_cpml))
  !allocate(iglob_CPML_temp(NGLLX*NGLLY*NGLLZ*nspec_cpml))
  !nglob_CPML = 0
  !ibool_CPML(:,:,:,:) = 0
  !do ispec_CPML = 1, nspec_cpml
  !  do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
  !    ispec = CPML_to_spec(ispec_CPML)
  !    iglob = ibool(i,j,k,ispec)
  !    POINT_EXIST = .false.
  !    do iglob_CPML = 1, nglob_CPML
  !      if (iglob == iglob_CPML_temp(iglob_CPML)) then
  !        POINT_EXIST = .true.
  !        exit
  !      endif
  !    enddo
  !    if (.not. POINT_EXIST) then
  !      nglob_CPML = nglob_CPML + 1
  !      iglob_CPML_temp(nglob_CPML) = iglob
  !      ibool_CPML(i,j,k,ispec_CPML) = nglob_CPML
  !    else
  !      ibool_CPML(i,j,k,ispec_CPML) = iglob_CPML  
  !    endif
  !  enddo;enddo;enddo
  !enddo
  !allocate(CPML_to_glob(nglob_CPML))
  !CPML_to_glob(1:nglob_CPML) = iglob_CPML_temp(1:nglob_CPML)
  !deallocate(iglob_CPML_temp)

  !!! MPI interface
  !num_interfaces_PML = 0
  !max_nibool_interfaces_PML = 0
  !if (num_interfaces_ext_mesh > 0) then
  !  !! preliminary counting
  !  allocate(nibool_interfaces_PML_rough(num_interfaces_ext_mesh))
  !  nibool_interfaces_PML_rough(:) = 0
  !  do iface = 1, num_interfaces_ext_mesh
  !    do igll = 1, nibool_interfaces_ext_mesh(iface)
  !      iglob = ibool_interfaces_ext_mesh(igll,iface)
  !      xc = xstore(iglob);yc = ystore(iglob);zc = zstore(iglob)
  !      if (is_point_inside_PML(xc,yc,zc)) &
  !        nibool_interfaces_PML_rough(iface) = &
  !            nibool_interfaces_PML_rough(iface)+1
  !    enddo
  !    if (nibool_interfaces_PML_rough(iface) > 0) &
  !      num_interfaces_PML = num_interfaces_PML + 1
  !  enddo
  !  max_nibool_interfaces_PML = maxval(nibool_interfaces_PML_rough(:))
  !  if (max_nibool_interfaces_PML > 0) then
  !    !! we actually have PML interfaces, start actual indexing
  !    allocate(nibool_interfaces_PML(num_interfaces_PML))
  !    allocate(ibool_interfaces_PML(max_nibool_interfaces_PML,&
  !                                             num_interfaces_PML))
  !    allocate(my_neighbors_PML(num_interfaces_PML))
  !    nibool_interfaces_PML(:) = 0
  !    num_interfaces_PML = 0
  !    do iface = 1, num_interfaces_ext_mesh
  !      if (nibool_interfaces_PML_rough(iface) .eq. 0) cycle
  !      num_interfaces_PML = num_interfaces_PML + 1
  !      do igll = 1, nibool_interfaces_ext_mesh(iface)
  !        iglob = ibool_interfaces_ext_mesh(igll,iface)
  !        xc = xstore(iglob);yc = ystore(iglob);zc = zstore(iglob)
  !        if (is_point_inside_PML(xc,yc,zc)) then
  !          nibool_interfaces_PML(num_interfaces_PML) = &
  !            nibool_interfaces_PML(num_interfaces_PML)+1
  !          POINT_EXIST = .false.
  !          do iglob_CPML = 1, nglob_CPML
  !            if (iglob == iglob_CPML_temp(iglob_CPML)) then
  !              ibool_interfaces_PML(nibool_interfaces_PML(num_interfaces_PML)&
  !                    ,num_interfaces_PML)=iglob_CPML
  !              POINT_EXIST = .true.
  !              exit
  !            endif
  !          enddo
  !          if (.not. POINT_EXIST) print *, 'PML GLL wrong indexing'
  !        endif
  !      enddo
  !      my_neighbors_PML(num_interfaces_PML) = my_neighbors_ext_mesh(iface)
  !      if (.not. (nibool_interfaces_PML(num_interfaces_PML) .eq. &
  !              nibool_interfaces_PML_rough(iface))) &
  !          print *, 'indexing not match for PML interface'
  !    enddo
  !  endif
  !  deallocate(nibool_interfaces_PML_rough)
  !endif
  !deallocate(iglob_CPML_temp)

  !!!!!!!!!!!!!!!!!
  contains
    logical function is_point_inside_PML(xp,yp,zp)
      real,intent(in) :: xp,yp,zp
      is_point_inside_PML = .false.
      if ((xp<(xoriginleft-xtol)).or.(xp>(xoriginright+xtol))&
          .or.(yp<(yoriginback-ytol)).or.(yp>(yoriginfront+ytol))&
          .or.(zp<(zoriginbottom-ztol)).or.(zp>(zorigintop+ztol))) then
        if ((xp>(x_min_all+xtol)).and.(xp<(x_max_all-xtol)) &
            .and.(yp>(y_min_all+ytol)).and.(yp<(y_max_all-ytol)) &
            .and.(zp>(z_min_all+ztol)).and.(zp<(z_max_all-ztol))) &
          is_point_inside_PML = .true.
      endif
    end function is_point_inside_PML
    logical function is_interface_in(xp,yp,zp)
      real, intent(in) :: xp,yp,zp
      is_interface_in = .false.
      !! internal interface
      if (((abs(xp-xoriginleft)<xtol).or.(abs(xp-xoriginright)<xtol)).and.&
         (yp.gt.yoriginback).and.(yp.lt.yoriginfront).and.&
         (zp.gt.zoriginbottom).and.(zp.lt.zorigintop)) &
        is_interface_in = .true.
      if (((abs(yp-yoriginback)<ytol).or.(abs(yp-yoriginfront)<ytol)).and.&
         (xp.gt.xoriginleft).and.(xp.lt.xoriginright).and.&
         (zp.gt.zoriginbottom).and.(zp.lt.zorigintop)) &
        is_interface_in = .true.
      if (((abs(zp-zoriginbottom)<ztol).or.(abs(zp-zorigintop)<ztol)).and.&
         (xp.gt.xoriginleft).and.(xp.lt.xoriginright).and.&
         (yp.gt.yoriginback).and.(yp.lt.yoriginfront)) &
        is_interface_in = .true.
    end function is_interface_in
    logical function is_interface_out(xp,yp,zp)
      real, intent(in) :: xp,yp,zp
      is_interface_out = .false.
      !! external interface
      if ((abs(xp-x_min_all)<xtol).or.(abs(xp-x_max_all)<xtol)&
            .or.(abs(yp-y_min_all)<ytol).or.(abs(yp-y_max_all)<ytol)&
            .or.(abs(zp-z_min_all)<ztol).or.(abs(zp-z_max_all)<ztol)) &
        is_interface_out = .true.
    end function is_interface_out
    subroutine include_this_face_out(icor,igll_face)
      integer, intent(in) :: icor,igll_face
      integer :: igll_temp
      num_pml_out = num_pml_out + 1
      pml_out_ispec(num_pml_out) = ispec_CPML
      igll_temp = 0
      if (icor == 3) then
        do i = 1,NGLLX; do j = 1,NGLLY
          igll_temp = igll_temp + 1
          pml_out_ijk(1,igll_temp,num_pml_out) = i
          pml_out_ijk(2,igll_temp,num_pml_out) = j
          pml_out_ijk(3,igll_temp,num_pml_out) = igll_face
        enddo; enddo
      elseif (icor == 2) then
        do i = 1,NGLLX; do k = 1,NGLLZ
          igll_temp = igll_temp + 1
          pml_out_ijk(1,igll_temp,num_pml_out) = i
          pml_out_ijk(2,igll_temp,num_pml_out) = igll_face
          pml_out_ijk(3,igll_temp,num_pml_out) = k
        enddo; enddo
      elseif (icor == 1) then
        do j = 1,NGLLY; do k = 1,NGLLZ
          igll_temp = igll_temp + 1
          pml_out_ijk(1,igll_temp,num_pml_out) = igll_face
          pml_out_ijk(2,igll_temp,num_pml_out) = j
          pml_out_ijk(3,igll_temp,num_pml_out) = k
        enddo; enddo
      endif
    end subroutine include_this_face_out
    !subroutine include_this_face_in(icor,igll_face)
    !  integer, intent(in) :: icor,igll_face
    !  integer :: igll_temp
    !  num_pml_in = num_pml_in + 1
    !  pml_in_ispec(num_pml_in) = ispec_CPML
    !  igll_temp = 0
    !  if (icor == 3) then
    !    do i = 1,NGLLX; do j = 1,NGLLY
    !      igll_temp = igll_temp + 1
    !      pml_in_ijk(1,igll_temp,num_pml_in) = i
    !      pml_in_ijk(2,igll_temp,num_pml_in) = j
    !      pml_in_ijk(3,igll_temp,num_pml_in) = igll_face
    !    enddo; enddo
    !  elseif (icor == 2) then
    !    do i = 1,NGLLX; do k = 1,NGLLZ
    !      igll_temp = igll_temp + 1
    !      pml_in_ijk(1,igll_temp,num_pml_in) = i
    !      pml_in_ijk(2,igll_temp,num_pml_in) = igll_face
    !      pml_in_ijk(3,igll_temp,num_pml_in) = k
    !    enddo; enddo
    !  elseif (icor == 1) then
    !    do j = 1,NGLLY; do k = 1,NGLLZ
    !      igll_temp = igll_temp + 1
    !      pml_in_ijk(1,igll_temp,num_pml_in) = igll_face
    !      pml_in_ijk(2,igll_temp,num_pml_in) = j
    !      pml_in_ijk(3,igll_temp,num_pml_in) = k
    !    enddo; enddo
    !  endif
    !end subroutine include_this_face_in
end subroutine adepml_set_local_dampingcoeff


