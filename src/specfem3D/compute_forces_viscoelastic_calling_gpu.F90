
! elastic solver


subroutine compute_forces_viscoelastic_ADE_GPU_calling()
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par

  !nqdu added
  use wavefield_discontinuity_par,only: IS_WAVEFIELD_DISCONTINUITY

  implicit none

  !local 
  integer:: iphase
  logical :: ADE_CONTRIB
  
  ADE_CONTRIB = PML_CONDITIONS .and. USE_ADE_PML.and.(nglob_CPML>0)

  !! Tianshi Liu: for solving wavefield discontinuity problem with
  !! non-split-node scheme
  !nqdu if (IS_WAVEFIELD_DISCONTINUITY) then
  if (IS_WAVEFIELD_DISCONTINUITY .and. COUPLE_WITH_INJECTION_TECHNIQUE) then
    call read_wavefield_discontinuity_file()
    call transfer_wavefield_discontinuity_to_GPU()
  endif

  do iphase = 1,2 
    ! contains both forward SIM_TYPE==1 and backward SIM_TYPE==3 simulations
    call compute_forces_viscoelastic_cuda_ade(Mesh_pointer, iphase, deltat, &
                                          nspec_outer_elastic, &
                                          nspec_inner_elastic, &
                                          COMPUTE_AND_STORE_STRAIN,ATTENUATION,ANISOTROPY)

    ! computes additional contributions
    if(iphase == 1)then

      if (IS_WAVEFIELD_DISCONTINUITY .and. COUPLE_WITH_INJECTION_TECHNIQUE) then
        call add_traction_discontinuity_GPU()
      endif
    
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                                            SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                                            b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
                                            Mesh_pointer)
      endif

      ! acoustic coupling
      if (ACOUSTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          call compute_coupling_el_ac_cuda(Mesh_pointer,iphase, &
                                          num_coupling_ac_el_faces)
        endif
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION) then
        ! note:
        ! these routines are not implemented as CUDA kernels, we just transfer the fields
        ! from the GPU to the CPU and vice versa

        ! transfers displacement & acceleration to the CPU
        call transfer_displ_from_device(NDIM*NGLOB_AB,displ, Mesh_pointer)
        call transfer_accel_from_device(NDIM*NGLOB_AB,accel, Mesh_pointer)

        call compute_coupling_viscoelastic_po(iphase)

        ! transfers acceleration back to GPU
        call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
      endif

      ! adds source term (single-force/moment-tensor solution)
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !call transfer_accel_from_device(NGLOB_AB*NDIM,accel,Mesh_pointer)
      call compute_add_sources_viscoelastic_GPU()
    endif 

    ! communication
#ifndef USE_CUDA_AWARE_MPI
    if(iphase == 1) then  
      ! copy accel boundary to accel->d_send_buffer -> (async) h_send_buf
      call sync_accel_bdry_buffers(Mesh_pointer,iphase,buffer_send_vector_ext_mesh)
      call assemble_MPI_vector_send_cuda(NPROC, &
                  buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                  num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                  nibool_interfaces_ext_mesh, &
                  my_neighbors_ext_mesh, &
                  request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
        
      if(ADE_CONTRIB) then 
        call sync_ade_bdry_buffers(Mesh_pointer,iphase,buffer_send_matrix_PML)
        call assemble_MPI_matrix_send_cuda(NPROC,buffer_send_matrix_PML,&
                                          buffer_recv_matrix_PML,num_interfaces_PML,&
                                          max_nibool_interfaces_PML,nibool_interfaces_PML,&
                                          my_neighbors_PML,&
                                          request_send_matrix_PML,request_recv_matrix_PML)
      endif
    
    else 
      call assemble_MPI_accel_update(NPROC,Mesh_pointer,num_interfaces_ext_mesh,&
                                     max_nibool_interfaces_ext_mesh,&
                                     buffer_recv_vector_ext_mesh,&
                                     request_send_vector_ext_mesh,&
                                     request_recv_vector_ext_mesh)

      !ADE Qt_t
      if(ADE_CONTRIB)  then
        call assemble_MPI_ADE_update(NPROC,Mesh_pointer,num_interfaces_PML,&
                                  max_nibool_interfaces_PML,&
                                  buffer_recv_matrix_PML,&
                                  request_send_matrix_PML,&
                                  request_recv_matrix_PML)
      endif
    endif
#else  
    ! copy accel boundary to accel->d_send_buffer -> (async) h_send_buf
    call sync_accel_bdry_buffers(Mesh_pointer,iphase,my_neighbors_ext_mesh,&
                                 nibool_interfaces_ext_mesh)
    if(ADE_CONTRIB) then 
      call sync_ade_bdry_buffers(Mesh_pointer,iphase,my_neighbors_PML,&
                                nibool_interfaces_PML)
    endif
#endif
  enddo      

  ! corrector to update PML auxiliary variables
  if(ADE_CONTRIB)  then
    call include_adepml_accel_aux_GPU() 
  endif 

  ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  ! dirichlet boundary condition has been added to rmass
  call apply_massmat_device(Mesh_pointer)

  ! update velocity
  call update_velocity_device(Mesh_pointer,deltatover2)

  if(ADE_CONTRIB)  then
    call update_Qt_conv_GPU()
    call update_Qu_conv_GPU()
  endif
end subroutine compute_forces_viscoelastic_ADE_GPU_calling

! subroutine compute_forces_viscoelastic_ADE_GPU_calling()

!   use specfem_par
!   use specfem_par_acoustic
!   use specfem_par_elastic
!   use specfem_par_poroelastic
!   use pml_par
!   use fault_solver_dynamic, only: bc_dynflt_set3d_all,SIMULATION_TYPE_DYN,synchronize_GPU
!   use fault_solver_kinematic, only: bc_kinflt_set_all,SIMULATION_TYPE_KIN
  
!   !nqdu added
!   use wavefield_discontinuity_par,only: IS_WAVEFIELD_DISCONTINUITY

!   implicit none

!   integer:: iphase,i

!   ! check
!   !if (PML_CONDITIONS) call exit_MPI(myrank,'PML conditions not yet implemented on GPUs')

!   !! Tianshi Liu: for solving wavefield discontinuity problem with
!   !! non-split-node scheme
!   !nqdu if (IS_WAVEFIELD_DISCONTINUITY) then
!   if (IS_WAVEFIELD_DISCONTINUITY .and. COUPLE_WITH_INJECTION_TECHNIQUE) then
!     call read_wavefield_discontinuity_file()
!     call transfer_wavefield_discontinuity_to_GPU()
!   endif

!   ! distinguishes two runs: for elements in contact with MPI interfaces, and elements within the partitions
!   do iphase = 1,2

!     ! contains both forward SIM_TYPE==1 and backward SIM_TYPE==3 simulations
!     call compute_forces_viscoelastic_cuda_ade(Mesh_pointer, iphase, deltat, &
!                                           nspec_outer_elastic, &
!                                           nspec_inner_elastic, &
!                                           COMPUTE_AND_STORE_STRAIN,ATTENUATION,ANISOTROPY)
!     ! computes additional contributions
!     if (iphase == 1) then

!       if (IS_WAVEFIELD_DISCONTINUITY .and. COUPLE_WITH_INJECTION_TECHNIQUE) then
!           call add_traction_discontinuity_GPU()
!       endif
      
!       ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
!       if (STACEY_ABSORBING_CONDITIONS) then
!         call compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
!                                              SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
!                                              b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
!                                              Mesh_pointer)
!       endif

!       ! acoustic coupling
!       if (ACOUSTIC_SIMULATION) then
!         if (num_coupling_ac_el_faces > 0) then
!           call compute_coupling_el_ac_cuda(Mesh_pointer,iphase, &
!                                            num_coupling_ac_el_faces)
!         endif
!       endif

!       ! poroelastic coupling
!       if (POROELASTIC_SIMULATION) then
!         ! note:
!         ! these routines are not implemented as CUDA kernels, we just transfer the fields
!         ! from the GPU to the CPU and vice versa

!         ! transfers displacement & acceleration to the CPU
!         call transfer_displ_from_device(NDIM*NGLOB_AB,displ, Mesh_pointer)
!         call transfer_accel_from_device(NDIM*NGLOB_AB,accel, Mesh_pointer)

!         call compute_coupling_viscoelastic_po(iphase)

!         ! transfers acceleration back to GPU
!         call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
!       endif

!       ! adds source term (single-force/moment-tensor solution)
!       ! note: we will add all source contributions in the first pass, when iphase == 1
!       !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
!       !call transfer_accel_from_device(NGLOB_AB*NDIM,accel,Mesh_pointer)
!       call compute_add_sources_viscoelastic_GPU()
!       ! call compute_add_sources_viscoelastic()
!       ! call transfer_accel_to_device(NGLOB_AB*NDIM,accel,Mesh_pointer)
    
!     endif
  
!     if(iphase == 2) then 
!       ! while inner elements compute "Kernel_2", we wait for MPI to
!       ! finish and transfer the boundary terms to the device asynchronousl
!       !daniel: todo - this avoids calling the Fortran vector send from CUDA routine
!       ! wait for asynchronous copy to finish
!       call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_vector_ext_mesh)

!       ! sends MPI buffers
!       call assemble_MPI_vector_send_cuda(NPROC, &
!                   buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
!                   num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
!                   nibool_interfaces_ext_mesh, &
!                   my_neighbors_ext_mesh, &
!                   request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
  
!       if(PML_CONDITIONS .and. USE_ADE_PML.and.(nglob_CPML>0)) then 
!         call sync_copy_Qt_t_from_device(Mesh_pointer,iphase,buffer_send_matrix_PML)
    
!         call assemble_MPI_matrix_send_cuda(NPROC,buffer_send_matrix_PML,&
!                                           buffer_recv_matrix_PML,num_interfaces_PML,&
!                                           max_nibool_interfaces_PML,nibool_interfaces_PML,&
!                                           my_neighbors_PML,&
!                                           request_send_matrix_PML,request_recv_matrix_PML)
!       endif

!       ! transfers MPI buffers onto GPU
!       call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_ext_mesh, &
!                   num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
!                   request_recv_vector_ext_mesh)
!       if(PML_CONDITIONS .and. USE_ADE_PML.and.(nglob_CPML>0)) then 
!         call transfer_ade_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_matrix_PML,&
!                                               num_interfaces_PML,max_nibool_interfaces_PML,&
!                                               request_recv_matrix_PML,request_send_matrix_PML)
!       endif

!     endif ! iphase


!     ! assemble all the contributions between slices using MPI
!     if (iphase == 1) then
!       ! sends accel values to corresponding MPI interface neighbors

!       ! transfers boundary region to host asynchronously. The
!       ! MPI-send is done from within compute_forces_viscoelastic_cuda,
!       ! once the inner element kernels are launched, and the
!       ! memcpy has finished. see compute_forces_viscoelastic_cuda: ~ line 1655
!       call transfer_boundary_from_device_a(Mesh_pointer,nspec_outer_elastic)
!       if(PML_CONDITIONS .and. USE_ADE_PML.and.(nglob_CPML>0))  then
!         call transfer_ade_boundary_from_device(Mesh_pointer,nspec_outer_elastic)
!       endif

!     else
!       ! waits for send/receive requests to be completed and assembles values
!       call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB,accel, Mesh_pointer, &
!                       buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
!                       max_nibool_interfaces_ext_mesh, &
!                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
!                       request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
!                       1)
!       if(PML_CONDITIONS .and. USE_ADE_PML.and.(nglob_CPML>0)) &
!         call assemble_MPI_matrix_write_cuda(NPROC,Mesh_pointer,&
!                                             num_interfaces_PML,request_send_matrix_PML,&
!                                             request_recv_matrix_PML)
        
!     endif

!   enddo

!   !Percy , Fault boundary term B*tau is added to the assembled forces
!   !        which at this point are stored in the array 'accel'
!   if (SIMULATION_TYPE_DYN .or. SIMULATION_TYPE_KIN) then
!     ! transfers wavefields to the CPU
!     ! call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc,accel, Mesh_pointer)
!     ! will remove later if GPU fault solver is fully tested

!     ! adds dynamic source
!     ! if (SIMULATION_TYPE_DYN) call bc_dynflt_set3d_all(accel,veloc,displ)
!     ! if (SIMULATION_TYPE_KIN) call bc_kinflt_set_all(accel,veloc,displ)
!     call fault_solver_gpu(Mesh_pointer,Fault_pointer,deltat,myrank,it)  ! GPU fault solver
!     !call transfer_boundary_from_device_a(Mesh_pointer,nspec_outer_elastic)
!     ! transfer data from mp->d_boundary to mp->h_boundary
!     !call sync_copy_from_device(Mesh_pointer,2,buffer_send_vector_ext_mesh)
!     ! transfer data from mp->h_boundary to send_buffer
!     !call assemble_MPI_vector_send_cuda(NPROC, &
!     !              buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
!     !              num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
!     !              nibool_interfaces_ext_mesh, &
!     !              my_neighbors_ext_mesh, &
!     !              request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

!     ! transfers MPI buffers onto GPU
!     !call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_ext_mesh, &
!     !              num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
!     !              request_recv_vector_ext_mesh)
!     ! waits for send/receive requests to be completed and assembles values
!     !call synchronize_MPI_vector_write_cuda(NPROC,NGLOB_AB,accel, Mesh_pointer, &
!     !                  buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
!     !                  max_nibool_interfaces_ext_mesh, &
!     !                  nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
!     !                  request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
!     !                  1)

!     if (mod(it,500) == 0 .and. it /= 0) call synchronize_GPU(it)  ! output results every 500 steps

!     ! transfers acceleration back to GPU
!     ! call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
!     ! will remove later if GPU fault solver is fully tested
!   endif

!   ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
!   call apply_massmat_device(Mesh_pointer)

!   ! corrector to update PML auxiliary variables
!   if(PML_CONDITIONS .and. USE_ADE_PML.and.(nglob_CPML>0))  then
!     call include_adepml_accel_aux_GPU() 
!     call apply_dirichlet_condition(Mesh_pointer)
!   endif 

!   ! update velocity
!   call update_velocity_device(Mesh_pointer,deltatover2)

!   if(PML_CONDITIONS .and. USE_ADE_PML.and.(nglob_CPML>0))  then
!     call update_Qt_conv_GPU()
!     call update_Qu_conv_GPU()
!   endif
  
!   ! updates acceleration with ocean load term
!   if (APPROXIMATE_OCEAN_LOAD) then
!     call compute_coupling_ocean_cuda(Mesh_pointer)

!     ! updates velocities
!     ! Newmark finite-difference time scheme with elastic domains:
!     ! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!     !
!     ! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
!     ! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
!     ! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t))
!     !
!     ! where
!     !   u, v, a are displacement,velocity & acceleration
!     !   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!     !   f denotes a source term (acoustic/elastic)
!     !   chi_dot_dot is acoustic (fluid) potential ( dotted twice with respect to time)
!     !
!     ! corrector:
!     ! updates the velocity term which requires a(t+delta)
!     ! GPU_MODE: this is handled in 'kernel_3' at the same time as accel*rmass
!     call kernel_3_b_cuda(Mesh_pointer,deltatover2,b_deltatover2)
!   endif

! end subroutine compute_forces_viscoelastic_ADE_GPU_calling