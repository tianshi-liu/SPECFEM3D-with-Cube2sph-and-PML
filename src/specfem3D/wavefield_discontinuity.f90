!! Solving the wavefield discontinuity problem with a non-split-node
!! scheme
!! Tianshi Liu, 2023.5

subroutine read_mesh_databases_wavefield_discontinuity()
  use specfem_par, only: CUSTOM_REAL,&
                         NSPEC_AB, NGLLX, NGLLY, NGLLZ, NDIM, NGLLSQUARE
  use wavefield_discontinuity_solver
  implicit none
  integer :: IIN = 27 ! file number for proc*_external_mesh.bin in this version
  allocate(ispec_to_elem_wd(NSPEC_AB))
  read(IIN) ispec_to_elem_wd
  read(IIN) nglob_wd
  read(IIN) nspec_wd
  allocate(ibool_wd(NGLLX, NGLLY, NGLLZ, nspec_wd), &
           boundary_to_iglob_wd(nglob_wd), &
           mass_in_wd(nglob_wd))
  read(IIN) ibool_wd
  read(IIN) boundary_to_iglob_wd
  read(IIN) mass_in_wd
  read(IIN) nfaces_wd
  allocate(face_ijk_wd(NDIM, NGLLSQUARE, nfaces_wd), &
           face_ispec_wd(nfaces_wd), &
           face_normal_wd(NDIM, NGLLSQUARE, nfaces_wd), &
           face_jacobian2Dw_wd(NGLLSQUARE, nfaces_wd))
  read(IIN) face_ijk_wd
  read(IIN) face_ispec_wd
  read(IIN) face_normal_wd
  read(IIN) face_jacobian2Dw_wd
  allocate(displ_wd(NDIM, nglob_wd), & 
           accel_wd(NDIM, nglob_wd), &
           traction_wd(NDIM, NGLLSQUARE, nfaces_wd))
end subroutine read_mesh_databases_wavefield_discontinuity

subroutine open_wavefield_discontinuity_file()
  use specfem_par, only: prname,IMAIN,myrank,LOCAL_PATH
  use specfem_par,only: SIMULATION_TYPE,NDIM,NSTEP,DT,NGLLSQUARE,CUSTOM_REAL

  ! nqdu added 
  use wavefield_discontinuity_solver
  implicit none

  ! local
  logical :: iexist  
  integer :: ier,it 
  integer(kind=8) :: file_size, block_size 

  ! check if dt/nstep file exists
  inquire(file=trim(LOCAL_PATH)//'wavefield_discontinuity_info.txt', &
          exist=iexist)
  if(.not.iexist) then
    DT_wd = real(DT,kind=CUSTOM_REAL) 
    NSTEP_wd = NSTEP
    SAVE_DOWNSAMPLED_WD = .false.

    ! open it in the old way
    if(SIMULATION_TYPE == 1) then 
      open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
          file=trim(prname)//'wavefield_discontinuity.bin', &
          status='old',action='read',form='unformatted')
    else 
      open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
          file=trim(prname)//'wavefield_discontinuity.bin', &
          status='old',action='read',form='unformatted',&
          access='stream')
    endif

    ! allocate dummpy arrays 
    allocate(field_a_wd(1,1,1),field_d_wd(1,1,1),field_t_wd(1,1,1,1))

  else 
    open(unit=IFILE_WAVEFIELD_DISCONTINUITY+100, &
        file=trim(LOCAL_PATH)//'wavefield_discontinuity_info.txt', &
        status='old',action='read',form='formatted')
    read(IFILE_WAVEFIELD_DISCONTINUITY+100,*) DT_wd
    read(IFILE_WAVEFIELD_DISCONTINUITY+100,*) NSTEP_wd
    close(IFILE_WAVEFIELD_DISCONTINUITY+100)
    SAVE_DOWNSAMPLED_WD = .true.

    ! open file to get size 
    inquire(file=trim(prname)//'wavefield_discontinuity.bin', &
            exist=iexist, size=file_size)
    if(.not.iexist) then
      if(myrank == 0) &
        print*, 'wavefield_discontinuity.bin does not exist, please check!'
      stop
    endif

    ! compute block size
    block_size = (nglob_wd*NDIM*2 + nfaces_wd*NGLLSQUARE*NDIM)*CUSTOM_REAL + 2 * 3 * CUSTOM_REAL
    block_size = block_size * NSTEP_wd 

    ! warnings if memory usage is too large
    if(block_size > file_size / 5) then 
      if ( myrank == 0 ) then
        write(IMAIN,*) '****************************************************************'
        write(IMAIN,*) 'Warning: the size of downsampled wavefield_discontinuity.bin is much larger than expected!'
        write(IMAIN,*) '         Please check if DT and NSTEP in wavefield_discontinuity_info.txt are correct!'
        write(IMAIN,*) '         Current DT and DT_wd = ', DT,DT_wd
        write(IMAIN,*) '         Current NSTEP = ', NSTEP,NSTEP_wd
        write(IMAIN,*) '         total/downsample file size (bytes) = ', file_size, block_size
        write(IMAIN,*) '****************************************************************'
      end if
    endif

    ! allocate space 
    allocate(field_d_wd(NDIM,nglob_wd,NSTEP_wd), &
             field_a_wd(NDIM,nglob_wd,NSTEP_wd), &
             field_t_wd(NDIM,NGLLSQUARE,nfaces_wd,NSTEP_wd),&
             stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 89001')

    ! read everything into memory
    open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
        file=trim(prname)//'wavefield_discontinuity.bin', &
        status='old',action='read',form='unformatted')
    do it = 1, NSTEP_wd
      read(IFILE_WAVEFIELD_DISCONTINUITY) field_d_wd(:,:,it)
      read(IFILE_WAVEFIELD_DISCONTINUITY) field_a_wd(:,:,it)
      read(IFILE_WAVEFIELD_DISCONTINUITY) field_t_wd(:,:,:,it)
    enddo

    ! close in finalize
    !close(IFILE_WAVEFIELD_DISCONTINUITY)
  endif

end subroutine open_wavefield_discontinuity_file

subroutine read_wavefield_discontinuity_file()
  use wavefield_discontinuity_solver, only: IFILE_WAVEFIELD_DISCONTINUITY, &
                 displ_wd, accel_wd, traction_wd
  
  !nqdu added
  use wavefield_discontinuity_solver,only: nglob_wd,nfaces_wd,DT_wd,NSTEP_wd
  use wavefield_discontinuity_solver,only: field_a_wd,field_d_wd,field_t_wd,SAVE_DOWNSAMPLED_WD
  use specfem_par,only: NDIM,NGLLSQUARE,CUSTOM_REAL
  use specfem_par,only: it ,NSTEP,SIMULATION_TYPE,DT 

  ! use specfem_par, only: NDIM
  ! use wavefield_discontinuity_solver, only: nglob_wd, nfaces_wd
  implicit none
  integer(kind=8) :: offset(3),block_bytes
  real(kind=CUSTOM_REAL) :: coef,t_now
  integer :: it1,it2,it_sem 

  if(SAVE_DOWNSAMPLED_WD) then

    ! get current time 
    if (SIMULATION_TYPE == 1) then
      it_sem = it 
    else 
      it_sem = NSTEP - it + 1
    endif
    t_now = real((it_sem - 1) * DT,kind=CUSTOM_REAL)

    ! Interpolate to get the correct index in downsampled arrays
    it1 = int(t_now / DT_wd) + 1 
    it2 = it1 + 1
    coef = (t_now - real((it1 - 1) * DT_wd,kind=CUSTOM_REAL)) / DT_wd
    if(it2 > NSTEP_wd) then
      it2 = NSTEP_wd
      it1 = NSTEP_wd
      coef = 0.0_CUSTOM_REAL
    endif

    ! interpolate 
    displ_wd(:,:) = field_d_wd(:,:,it1) * (field_d_wd(:,:,it2)-field_d_wd(:,:,it1)) * coef
    accel_wd(:,:) = field_a_wd(:,:,it1) * (field_a_wd(:,:,it2)-field_a_wd(:,:,it1)) * coef
    traction_wd(:,:,:) = field_t_wd(:,:,:,it1) * (field_t_wd(:,:,:,it2)-field_t_wd(:,:,:,it1)) * coef
  
  else 
    if(SIMULATION_TYPE == 1) then
      read(IFILE_WAVEFIELD_DISCONTINUITY) displ_wd
      read(IFILE_WAVEFIELD_DISCONTINUITY) accel_wd
      read(IFILE_WAVEFIELD_DISCONTINUITY) traction_wd
    else 
      block_bytes = CUSTOM_REAL *(nglob_wd*NDIM*2 + nfaces_wd*NGLLSQUARE*NDIM) + 8 * 3
      offset(1) = (NSTEP - it) * block_bytes + 5 
      offset(2) = offset(1) + CUSTOM_REAL * nglob_wd*NDIM + 8
      offset(3) = offset(2) + CUSTOM_REAL * nglob_wd*NDIM + 8

      ! read
      read(IFILE_WAVEFIELD_DISCONTINUITY,pos=offset(1)) displ_wd
      read(IFILE_WAVEFIELD_DISCONTINUITY,pos=offset(2)) accel_wd
      read(IFILE_WAVEFIELD_DISCONTINUITY,pos=offset(3)) traction_wd
    endif
  endif
  ! if(size(displ_wd)> 0) then
  !   print*,maxval(displ_wd),maxval(traction_wd),maxval(accel_wd),trim(prname),maxloc(traction_wd)
  ! endif
end subroutine read_wavefield_discontinuity_file

subroutine finalize_wavefield_discontinuity()
  use wavefield_discontinuity_solver
  implicit none
  close(IFILE_WAVEFIELD_DISCONTINUITY)
  ! close(IFILE_WAVEFIELD_DISCONTINUITY+1)
  deallocate(ispec_to_elem_wd, ibool_wd, boundary_to_iglob_wd, mass_in_wd, &
             face_ijk_wd, face_ispec_wd, face_normal_wd, face_jacobian2Dw_wd, &
             displ_wd, accel_wd, traction_wd)

  ! free arrays
  deallocate(field_d_wd,field_a_wd,field_t_wd)
end subroutine finalize_wavefield_discontinuity

subroutine add_displacement_discontinuity_element(ispec, dummyx_loc, &
                                                  dummyy_loc, dummyz_loc)
  use specfem_par, only: CUSTOM_REAL,NGLLX, NGLLY, NGLLZ
  use wavefield_discontinuity_solver, only: ispec_to_elem_wd, ibool_wd, &
                                            displ_wd
  implicit none
  integer, intent(in) :: ispec
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: dummyx_loc, &
                        dummyy_loc, dummyz_loc
  integer :: ispec_wd, i, j, k, iglob_wd
  ispec_wd = ispec_to_elem_wd(ispec)
  if (ispec_wd /= 0) then
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob_wd = ibool_wd(i,j,k,ispec_wd)
          if (iglob_wd /= 0) then
            dummyx_loc(i,j,k) = dummyx_loc(i,j,k) + displ_wd(1, iglob_wd)
            dummyy_loc(i,j,k) = dummyy_loc(i,j,k) + displ_wd(2, iglob_wd)
            dummyz_loc(i,j,k) = dummyz_loc(i,j,k) + displ_wd(3, iglob_wd)
          endif
        enddo
      enddo
    enddo
  endif 
end subroutine add_displacement_discontinuity_element

subroutine add_traction_discontinuity()
  use wavefield_discontinuity_solver, only: accel_wd, mass_in_wd, nglob_wd, &
                                    boundary_to_iglob_wd, traction_wd,&
                                    nfaces_wd, face_ijk_wd, face_ispec_wd, &
                                    face_jacobian2Dw_wd
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NGLLSQUARE, ibool
  use specfem_par_elastic, only: accel,b_accel

  !nqdu added
  use specfem_par,only: SIMULATION_TYPE

  implicit none
  integer :: iglob_wd, iglob, ispec, i, j, k, iface_wd, igll
  real(kind=CUSTOM_REAL) :: jacobianw

  if(SIMULATION_TYPE == 1) then
    do iglob_wd = 1, nglob_wd
      iglob = boundary_to_iglob_wd(iglob_wd)
      accel(:,iglob) = accel(:,iglob) - &
                      accel_wd(:,iglob_wd) * mass_in_wd(iglob_wd)
    enddo
    do iface_wd = 1, nfaces_wd
      do igll = 1, NGLLSQUARE
        i = face_ijk_wd(1, igll, iface_wd)
        j = face_ijk_wd(2, igll, iface_wd)
        k = face_ijk_wd(3, igll, iface_wd)
        ispec = face_ispec_wd(iface_wd)
        iglob = ibool(i,j,k,ispec)
        jacobianw = face_jacobian2Dw_wd(igll, iface_wd)
        accel(:,iglob) = accel(:,iglob) + &
                      traction_wd(:,igll,iface_wd) * jacobianw
      enddo
    enddo
  else 
    do iglob_wd = 1, nglob_wd
      iglob = boundary_to_iglob_wd(iglob_wd)
      b_accel(:,iglob) = b_accel(:,iglob) - &
                      accel_wd(:,iglob_wd) * mass_in_wd(iglob_wd)
    enddo
    do iface_wd = 1, nfaces_wd
      do igll = 1, NGLLSQUARE
        i = face_ijk_wd(1, igll, iface_wd)
        j = face_ijk_wd(2, igll, iface_wd)
        k = face_ijk_wd(3, igll, iface_wd)
        ispec = face_ispec_wd(iface_wd)
        iglob = ibool(i,j,k,ispec)
        jacobianw = face_jacobian2Dw_wd(igll, iface_wd)
        b_accel(:,iglob) = b_accel(:,iglob) + &
                      traction_wd(:,igll,iface_wd) * jacobianw
      enddo
    enddo
  endif
end subroutine add_traction_discontinuity

subroutine transfer_wavefield_discontinuity_to_GPU()
  use specfem_par, only: Mesh_pointer, NDIM, NGLLSQUARE
  use wavefield_discontinuity_solver
  implicit none
  call transfer_wavefield_discontinuity_to_device(nglob_wd*NDIM, &
                                                  nfaces_wd*NDIM*NGLLSQUARE, &
                                                  displ_wd, accel_wd, &
                                                  traction_wd, Mesh_pointer)
end subroutine transfer_wavefield_discontinuity_to_GPU

subroutine prepare_wavefield_discontinuity_GPU()
  use specfem_par, only: Mesh_pointer
  use wavefield_discontinuity_solver
  implicit none
  call prepare_wavefield_discontinuity_device(Mesh_pointer, ispec_to_elem_wd, &
                                                nglob_wd, nspec_wd, ibool_wd, &
                                                boundary_to_iglob_wd, &
                                                mass_in_wd, &
                                                nfaces_wd, face_ijk_wd, &
                                                face_ispec_wd, face_normal_wd, &
                                                face_jacobian2dw_wd)
end subroutine prepare_wavefield_discontinuity_GPU



subroutine add_traction_discontinuity_GPU()
  use wavefield_discontinuity_solver
  use specfem_par, only: Mesh_pointer
  implicit none
  call wavefield_discontinuity_add_traction_cuda(nglob_wd, nfaces_wd, &
                                                 Mesh_pointer)
end subroutine add_traction_discontinuity_GPU