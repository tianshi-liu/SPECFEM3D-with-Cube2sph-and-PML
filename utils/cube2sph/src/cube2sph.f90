program cube2sph
  use mpi
  use constants
  use meshfem3D_par, only: LOCAL_PATH
  use meshfem3D_models_par, only: myrank
  use generate_databases_par, only: &
          prname, npts, nodes_coords, nmat_ext_mesh, nundefMat_ext_mesh, &
          materials_ext_mesh, undef_mat_prop, nelmnts_ext_mesh, &
          elmnts_ext_mesh, mat_ext_mesh, NSPEC_AB, &
          nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, &
          nspec2D_bottom_ext, nspec2D_top_ext, NSPEC2D_BOTTOM, NSPEC2D_TOP, &
          ibelm_xmin, nodes_ibelm_xmin, ibelm_xmax, nodes_ibelm_xmax, &
          ibelm_ymin, nodes_ibelm_ymin, ibelm_ymax, nodes_ibelm_ymax, &
          ibelm_bottom, nodes_ibelm_bottom, ibelm_top, nodes_ibelm_top, &
          nspec_cpml, nspec_cpml_tot, CPML_to_spec, CPML_regions, &
          is_CPML, NPROC, num_interfaces_ext_mesh, &
          max_interface_size_ext_mesh, &
          my_neighbors_ext_mesh, my_nelmnts_neighbors_ext_mesh, &
          my_interfaces_ext_mesh, ibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh, num_interface, SAVE_MOHO_MESH, &
          boundary_number, nspec2D_moho_ext, ibelm_moho, nodes_ibelm_moho
  implicit none
  integer :: ipt, ier
  character(len=MAX_STRING_LEN) :: arg_string
  double precision, dimension(:,:), allocatable :: nodes_coords_new
  double precision :: center_lat, center_lon, rotation_azi, &
     alpha, beta, gamma, cosa, sina, cosb, sinb, cosg, sing, t1, t2, &
     M11, M12, M13, M21, M22, M23, M31, M32, M33, xi, eta, x, y, z
  call get_command_argument(1, arg_string)
  read(arg_string, *) center_lat
  call get_command_argument(2, arg_string)
  read(arg_string, *) center_lon
  call get_command_argument(3, arg_string)
  read(arg_string, *) rotation_azi
  call MPI_Init(ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  call MPI_Comm_size(MPI_COMM_WORLD, NPROC, ier)
  LOCAL_PATH = 'DATABASES_MPI'
  call read_partition_files()
  allocate(nodes_coords_new(NDIM, npts), stat=ier)
  alpha = center_lon / RADIANS_TO_DEGREES
  beta = PI / 2.0 - center_lat / RADIANS_TO_DEGREES
  gamma = rotation_azi / RADIANS_TO_DEGREES
  cosa = cos(alpha)
  sina = sin(alpha)
  cosb = cos(beta)
  sinb = sin(beta)
  cosg = cos(gamma)
  sing = sin(gamma)
  M11 = cosg * cosb * cosa - sing * sina
  M12 = - sing* cosb * cosa - cosg * sina
  M13 = sinb * cosa
  M21 = cosg * cosb * sina + sing * cosa
  M22 = - sing * cosb * sina + cosg * cosa
  M23 = sinb * sina
  M31 = - cosg * sinb
  M32 = sing * sinb
  M33 = cosb
  do ipt = 1, npts
    x = nodes_coords(1, ipt)
    y = nodes_coords(2, ipt)
    z = nodes_coords(3, ipt)
    xi = x / R_EARTH
    eta = y / R_EARTH
    t1 = tan(xi)
    t2 = tan(eta)
    z = (1.0+z/R_EARTH)/sqrt(1.0+t1*t1+t2*t2)
    x = -z*t2*R_EARTH
    y = z*t1*R_EARTH
    z = z*R_EARTH
    nodes_coords_new(1, ipt) = M11 * x + M12 * y + M13 * z
    nodes_coords_new(2, ipt) = M21 * x + M22 * y + M23 * z
    nodes_coords_new(3, ipt) = M31 * x + M32 * y + M33 * z
  enddo
  nodes_coords(:,:) = nodes_coords_new(:,:)
  call write_partition_files()
  deallocate(nodes_coords_new)
  call MPI_Finalize(ier)
end program cube2sph
