program node_stretching
   use mpi
   use constants
   use meshfem3D_par, only: iregion_code, idoubling, MIN_ATTENUATION_PERIOD, &
           MAX_ATTENUATION_PERIOD, ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R120,R220,R400, &
                      R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
                      R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
                      RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER, &
                      ISOTROPIC_3D_MANTLE, USE_FULL_TISO_MANTLE, LOCAL_PATH

   use meshfem3D_models_par, only: myrank,CASE_3D, &
           REFERENCE_1D_MODEL, THREE_D_MODEL, ONE_CRUST, TRANSVERSE_ISOTROPY,&
           TOPOGRAPHY,ELLIPTICITY,CRUSTAL,ibathy_topo,nspl,rspl,espl,espl2
   implicit none
   integer :: npts, ipts, nspec, ispec, ipass, ipts_dummy, inode, ier
   !integer, parameter :: NGNOD = 27
   double precision, dimension(NGNOD) :: xelm, yelm, zelm
   double precision, dimension(:), allocatable :: xnode, ynode, znode, &
           xnode_new, ynode_new, znode_new
   character(len=MAX_STRING_LEN) :: infn, outfn, meshfn
   integer, dimension(:,:), allocatable :: element_node
   logical :: elem_in_crust, elem_in_mantle, elem_is_tiso
   !integer, dimension(:), allocatable :: idoubling
   logical, dimension(:), allocatable :: ispec_is_tiso
   double precision :: rmid
   call MPI_Init(ier)
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
   ! input xyz of nodes from file
   infn = 'nodes_coords_file_sph'
   outfn = 'nodes_coords_file_topo'
   open(unit=IIN, file=infn(1:len_trim(infn)), &
            status='unknown', form='formatted', action='read', iostat=ier)
   read(IIN, *) npts
   allocate(xnode(npts), ynode(npts), znode(npts))
   allocate(xnode_new(npts), ynode_new(npts), znode_new(npts))
   do ipts = 1, npts
     read(IIN, *) ipts_dummy, xnode(ipts), ynode(ipts), znode(ipts)
   enddo
   close(IIN)
   meshfn = 'mesh_file'
   open(unit=IIN, file=meshfn(1:len_trim(meshfn)), &
            status='unknown', form='formatted', action='read', iostat=ier)
   read(IIN, *) nspec
   allocate(element_node(NGNOD, nspec))
   do ispec = 1, nspec
     read(IIN, *) ipts_dummy, (element_node(inode, ispec), inode=1, NGNOD)
   enddo
   close(IIN)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !set parameters
   myrank = 0
   LOCAL_PATH = 'DATABASES_MPI'
   ELLIPTICITY = .true.
   TOPOGRAPHY = .true.
   CASE_3D = .true.
   CRUSTAL = .true.
   ISOTROPIC_3D_MANTLE = .true.
   ONE_CRUST = .true.
   !REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
   !THREE_D_MODEL = THREE_D_MODEL_GLL
   TRANSVERSE_ISOTROPY = .true.
   !CRUSTAL = .true.
   REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
   THREE_D_MODEL = THREE_D_MODEL_S29EA
   RMIDDLE_CRUST = 6356000.d0
   RMOHO = 6346600.d0 ! at 24.4km depth
   R80  = 6291000.d0
   R220 = 6151000.d0
   R400 = 5961000.d0 ! 410km discontinuity
   R600 = 5771000.d0
   R670 = 5721000.d0 ! 650km discontinuity
   R771 = 5600000.d0
   ! 3D models do not honor PREM moho but a fictitious moho at 40km depth:
   ! either to make simulation cheaper or to have a 3D crustal structure
   RMOHO_FICTITIOUS_IN_MESHER = (R80 + R_EARTH) / 2.0d0
   iregion_code = IREGION_CRUST_MANTLE
   R80_FICTITIOUS_IN_MESHER = R80
   if (CRUSTAL .and. CASE_3D) then
     !> Hejun
     ! mesh will honor 3D crustal moho topography
     ! moves MOHO up 5km to honor moho topography deeper than 35 km
     ! moves R80 down to 120km depth in order to have less squeezing for elements below moho
     RMOHO_FICTITIOUS_IN_MESHER = RMOHO_FICTITIOUS_IN_MESHER + RMOHO_STRETCH_ADJUSTMENT
     R80_FICTITIOUS_IN_MESHER = R80_FICTITIOUS_IN_MESHER + R80_STRETCH_ADJUSTMENT
   endif
   allocate(idoubling(nspec), ispec_is_tiso(nspec))
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! sets up spline coefficients for ellipticity
   if (ELLIPTICITY) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
   if (TOPOGRAPHY) then
     ! arrays for elevations
     allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
     if (ier /= 0) stop 'Error allocating ibathy_topo array'

     ! initializes
     ibathy_topo(:,:) = 0

     ! sets up topo/bathy
     call model_topo_bathy_broadcast(ibathy_topo,LOCAL_PATH)
   endif
   !call model_1dref_broadcast(CRUSTAL)
   !call model_s362ani_broadcast(THREE_D_MODEL)
   ! crustal model
   if (CRUSTAL) &
     call meshfem3D_crust_broadcast()
   do ispec = 1, nspec
     do inode = 1, NGNOD
       xelm(inode) = xnode(element_node(inode, ispec)) / R_EARTH
       yelm(inode) = ynode(element_node(inode, ispec)) / R_EARTH
       zelm(inode) = znode(element_node(inode, ispec)) / R_EARTH
     enddo
     rmid = sqrt(xelm(NGNOD) * xelm(NGNOD) + yelm(NGNOD) * yelm(NGNOD) &
                  + zelm(NGNOD) * zelm(NGNOD))
     if (rmid > RMOHO / R_EARTH) then
       idoubling(ispec) = IFLAG_CRUST
     else if (rmid > R80 / R_EARTH .and. rmid <= RMOHO / R_EARTH) then
       idoubling(ispec) = IFLAG_80_MOHO
     else if (rmid > R220 / R_EARTH .and. rmid <= R80 / R_EARTH) then
       idoubling(ispec) = IFLAG_220_80
     else if (rmid > R670 / R_EARTH .and. rmid <= R220 / R_EARTH) then
       idoubling(ispec) = IFLAG_670_220
     else
       idoubling(ispec) = IFLAG_MANTLE_NORMAL
     end if
     elem_in_crust = .false.
     elem_in_mantle = .false.
     elem_is_tiso = .false.
     if (iregion_code == IREGION_CRUST_MANTLE) then
       if (CRUSTAL .and. CASE_3D) then
         ! 3D crustal models
         if (idoubling(ispec) == IFLAG_CRUST &
           .or. idoubling(ispec) == IFLAG_220_80 &
           .or. idoubling(ispec) == IFLAG_80_MOHO) then
           ! Stretch mesh to honor smoothed moho thickness from crust2.0
           ! mesh is stretched between surface and 220
           !
           ! differentiate between regional and global meshing
           if (REGIONAL_MOHO_MESH) then
             call moho_stretching_honor_crust_reg(myrank,xelm,yelm,zelm, &
                                                  elem_in_crust,elem_in_mantle)
           else
             call moho_stretching_honor_crust(myrank,xelm,yelm,zelm, &
                                              elem_in_crust,elem_in_mantle)
           endif
         else
           ! element below 220km
           ! sets element flag for mantle
           elem_in_mantle = .true.
         endif
       else
         ! 1D crust, no stretching
         ! sets element flags
         if (idoubling(ispec) == IFLAG_CRUST) then
           elem_in_crust = .true.
         else
           elem_in_mantle = .true.
         endif
       endif
       ! sets transverse isotropic flag for elements in mantle
       if (TRANSVERSE_ISOTROPY) then
         ! modifies tiso to have it for all mantle elements
         ! preferred for example, when using 1Dref (STW model)
         if (USE_FULL_TISO_MANTLE) then
           ! all elements below the actual moho will be used for transverse isotropy
           ! note: this will increase the computation time by ~ 45 %
           if (elem_in_mantle) then
             elem_is_tiso = .true.
           endif
         else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
           ! transverse isotropic mantle between fictitious moho to 670km depth
           ! preferred for Harvard (Kustowski's) models using STW 1D reference, i.e.
           ! THREE_D_MODEL_S362ANI
           ! THREE_D_MODEL_S362WMANI
           ! THREE_D_MODEL_S29EA
           ! THREE_D_MODEL_GLL
           ! which show significant transverse isotropy also below 220km depth
           if (USE_OLD_VERSION_5_1_5_FORMAT) then
             ! assigns TI only to elements below (2-layer) fictitious moho down to 670
             if (idoubling(ispec) == IFLAG_220_80 &
               .or. idoubling(ispec) == IFLAG_80_MOHO &
               .or. idoubling(ispec) == IFLAG_670_220) then
               elem_is_tiso = .true.
             endif
           else
             ! assigns TI to elements in mantle elements just below actual moho down to 670
             if (idoubling(ispec) == IFLAG_220_80 &
               .or. idoubling(ispec) == IFLAG_80_MOHO &
               .or. idoubling(ispec) == IFLAG_670_220 &
               .or. (idoubling(ispec) == IFLAG_CRUST .and. elem_in_mantle )) then
               elem_is_tiso = .true.
             endif
           endif
         else if (idoubling(ispec) == IFLAG_220_80 .or. idoubling(ispec) == IFLAG_80_MOHO) then
           ! default case for PREM reference models:
           ! models use only transverse isotropy between moho and 220 km depth
           elem_is_tiso = .true.
           ! checks mantle flag to be sure
           !if (elem_in_mantle .eqv. .false. ) stop 'Error mantle flag confused between moho and 220'
         endif
       endif
   
     endif ! IREGION_CRUST_MANTLE

     ! sets element tiso flag
     ispec_is_tiso(ispec) = elem_is_tiso
     ! adds surface topography
     if (TOPOGRAPHY) then
       if (idoubling(ispec) == IFLAG_CRUST .or. &
           idoubling(ispec) == IFLAG_220_80 .or. &
           idoubling(ispec) == IFLAG_80_MOHO) then
         ! stretches mesh between surface and R220 accordingly
         ! stretches anchor points only, interpolates GLL points later on
         call add_topography(xelm,yelm,zelm,ibathy_topo)
       endif
     endif
     ! make the Earth elliptical
     if (ELLIPTICITY) then
       ! note: after adding ellipticity, the mesh becomes elliptical and geocentric and geodetic/geographic colatitudes differ.
       ! make the Earth's ellipticity, use element anchor points
       call get_ellipticity(xelm,yelm,zelm,nspl,rspl,espl,espl2)
     endif


     do inode = 1, NGNOD
       xnode_new(element_node(inode, ispec)) = xelm(inode) * R_EARTH
       ynode_new(element_node(inode, ispec)) = yelm(inode) * R_EARTH
       znode_new(element_node(inode, ispec)) = zelm(inode) * R_EARTH
     enddo
   enddo
   ! output xyz to file
   open(unit=IOUT, file=outfn(1:len_trim(outfn)), &
            status='unknown', form='formatted', action='write', iostat=ier)
   write(IOUT, *) npts
   do ipts = 1, npts
     write(IOUT, "(i10,f21.6,f21.6,f21.6)") ipts, xnode_new(ipts), &
                   ynode_new(ipts), znode_new(ipts)
   enddo
   close(IOUT)
   !!!!!!!!!!!!!!!!!!!!
   deallocate(element_node, xnode, ynode, znode, &
              xnode_new, ynode_new, znode_new, idoubling, ispec_is_tiso)
   call MPI_Finalize(ier)
end program node_stretching
