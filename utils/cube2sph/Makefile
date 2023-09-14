CC=icc
FC=ifort
MPIFC=mpif90

FCFLAGS=-O3 -traceback -module ./obj -I./obj -I. -I./include
CFLAGS=-O3 -traceback -I./include

NCDFCFLAGS=-I/usr/local/include
NCDFLDFLAGS=-L/usr/local/lib
NCDFLDLIBS=-lnetcdff

FC_MODEXT=mod
FC_MODDIR=./obj

O=./obj
E=./bin

get_moho_topo_obj= \
	$O/meshfem3D_models.check.o \
	$O/meshfem3D_par.check_module.o \
	$O/model_1dref.check.o \
	$O/model_1066a.check.o \
	$O/model_ak135.check.o \
	$O/model_sea1d.check.o \
	$O/model_aniso_inner_core.check.o \
	$O/model_aniso_mantle.check.o \
	$O/model_atten3D_QRFSI12.check.o \
	$O/model_attenuation.check.o \
	$O/model_crust_1_0.check.o \
	$O/model_crust_2_0.check.o \
	$O/model_crustmaps.check.o \
	$O/model_eucrust.check.o \
	$O/model_epcrust.check.o \
	$O/model_full_sh.check.o \
	$O/model_gapp2.check.o \
	$O/model_gll.check.o \
	$O/model_heterogen_mantle.check.o \
	$O/model_iasp91.check.o \
	$O/model_jp1d.check.o \
	$O/model_jp3d.check.o \
	$O/model_ppm.check.o \
	$O/model_s20rts.check.o \
	$O/model_s40rts.check.o \
	$O/model_s362ani.check.o \
	$O/model_sea99_s.check.o \
	$O/lgndr.check.o
setup_model_obj= \
	${get_moho_topo_obj} \
	$O/add_topography_410_650.check.o \
	$O/calc_jacobian.check.o \
	$O/compute_element_properties.check.o \
	$O/get_model.check.o \
	$O/get_shape3D.check.o \
	$O/get_ellipticity.check.o \
	$O/add_topography.check.o \
	$O/moho_stretching.check.o

node_stretching_obj= \
        ${get_moho_topo_obj} \
        $O/moho_stretching.check.o \
        $O/get_ellipticity.check.o \
        $O/add_topography.check.o

write_topo_obj= \
	${get_moho_topo_obj} \
	$O/specfem3D_par.solverstatic_module.o \
	$O/add_topography.check.o

write_boundary_obj= \
	${get_moho_topo_obj} \
	$O/netcdf_model_par.o \
	$O/get_boundary_from_netcdf.o \
	$O/get_boundary_from_crust_1_0.o \
	$O/netcdf_model_wrapper.o

write_model_tomo_obj= \
	${get_moho_topo_obj} \
	$O/specfem3D_par.solverstatic_module.o \
	$O/add_topography.check.o \
	$O/netcdf_model_par.o \
	$O/get_boundary_from_netcdf.o \
	$O/get_boundary_from_crust_1_0.o \
	$O/get_model_mantle.o \
	$O/netcdf_model_wrapper.o

write_force_solution_file_obj = \
	$O/get_force.solver.o \
	$O/get_cmt.solver.o \
	$O/specfem3D_par.solverstatic_module.o


write_cmt_solution_file_obj = \
	$O/get_force.solver.o \
	$O/get_cmt.solver.o \
	$O/specfem3D_par.solverstatic_module.o

write_stations_file_obj = \
	$O/specfem3D_par.solverstatic_module.o

shared_obj= \
	$O/calendar.shared.o \
	$O/shared_par.shared_module.o \
	$O/create_name_database.shared.o \
	$O/create_name_database_cartesian.shared.o \
	$O/binary_c_io.cc.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/gll_library.shared.o \
	$O/intgrl.shared.o \
	$O/model_prem.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/make_ellipticity.shared.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/smooth_weights_vec.shared.o \
	$O/spline_routines.shared.o \
	$O/heap_sort.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/lagrange_poly.shared.o \
	$O/get_all_eight_slices.shared.o \
	$O/hex_nodes.shared.o \
	$O/rthetaphi_xyz.shared.o

adios_shared_STUBS= \
	$O/adios_method_stubs.cc.o
generate_databases_obj=\
	$O/meshfem3D_par.check_module.o\
	$O/wavefield_discontinuity_cube2sph_mod.gen.o \
	$O/read_write_wavefield_discontinuity.gen.o \
	$O/read_partition_files.gen.o \
	$O/write_partition_files.gen.o \
	$O/generate_databases_par.gen_mod.o

setup_adepml_obj=\
	$O/get_gll_xyz.o\
	$O/get_shape3D.check.o \
	$O/calc_adepml_physical_jacobian.o\
	$O/create_mass_matrices_pml_elastic_loc.o

cube2sph_adepml_obj=\
        $O/get_gll_xyz.o\
	$O/cube2sph_trans.o\
	$O/get_shape3D.check.o \
        $O/calc_adepml_physical_jacobian.o\
        $O/create_mass_matrices_pml_elastic_loc.o\
	$O/calc_r_trans.o\
	$O/update_parameters_from_netcdf.o \
	$O/netcdf_model_wrapper.o \
	$O/find_parameter_value.o \
	$O/setup_model_cartesian_subroutine.o \
	$O/node_stretching_subroutine.o 

cube2sph_station_obj = \
	$O/cube2sph_trans.o

cube2sph_txt_obj = \
	$O/cube2sph_trans.o

create_ckbd_obj= \
	$O/meshfem3D_par.check_module.o \
	$O/specfem3D_par.solverstatic_module.o \
	$O/get_ellipticity.check.o \
	$O/add_topography.check.o


## rules
all: $E/setup_adepml $E/cube2sph $E/cube2sph_station_inv $E/cube2sph_boundary_gmt $E/cube2sph_force $E/cube2sph_station $E/write_horizontal_slice_vtk $E/write_vertical_slice_vtk $E/write_force_solution_file $E/write_cmt_solution_file $E/write_stations_file

$E/setup_adepml: ./src/setup_adepml.f90 $(shared_obj) $(adios_shared_STUBS) $(generate_databases_obj) $(setup_adepml_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph: ./src/cube2sph.f90 $(setup_model_obj) $(node_stretching_obj) $(shared_obj) $(adios_shared_STUBS) $(generate_databases_obj) $(cube2sph_adepml_obj)
	${MPIFC} ${FCFLAGS} ${NCDFLDFLAGS} -o $@ $^ ${NCDFLDLIBS}

$E/cube2sph_adepml: ./src/cube2sph_adepml.f90 $(setup_model_obj) $(node_stretching_obj) $(shared_obj) $(adios_shared_STUBS) $(generate_databases_obj) $(cube2sph_adepml_obj)
	${MPIFC} ${FCFLAGS} ${NCDFLDFLAGS} -o $@ $^ ${NCDFLDLIBS}

$E/cube2sph_nopml: ./src/cube2sph_nopml.f90 $(setup_model_obj) $(node_stretching_obj) $(shared_obj) $(adios_shared_STUBS) $(generate_databases_obj) $(cube2sph_adepml_obj)
	${MPIFC} ${FCFLAGS} ${NCDFLDFLAGS} -o $@ $^ ${NCDFLDLIBS}

$E/cube2sph_force: ./src/cube2sph_force.f90 $(cube2sph_station_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph_station: ./src/cube2sph_station.f90 $(cube2sph_station_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph_txt: ./src/cube2sph_txt.f90 $(cube2sph_txt_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/mask_depth: ./src/mask_depth.f90
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/add_topo_to_depth: ./src/add_topo_to_depth.f90
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/write_topo_txt: ./src/write_topo_txt.f90 $(write_topo_obj) $(shared_obj) $(adios_shared_STUBS)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph_boundary_gmt: ./src/cube2sph_boundary_gmt.f90 $(cube2sph_station_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph_station_inv: ./src/cube2sph_station_inv.f90 $(cube2sph_station_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/create_ckbd: ./src/create_ckbd.f90 $(create_ckbd_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$O/%.check_module.o: ./src/meshfem3D/%.f90 $O/shared_par.shared_module.o
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.check.o: ./src/meshfem3D/%.f90 $O/shared_par.shared_module.o $O/meshfem3D_par.check_module.o
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.check.o: ./src/meshfem3D/%.F90 $O/shared_par.shared_module.o $O/meshfem3D_par.check_module.o
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.checkmpi.o: ./src/meshfem3D/%.f90 $O/shared_par.shared_module.o $O/meshfem3D_par.check_module.o
	${MPIFC} ${FCFLAGS} -c -o $@ $<

$O/%.checkmpi.o: ./src/meshfem3D/%.F90 $O/shared_par.shared_module.o $O/meshfem3D_par.check_module.o
	${MPIFC} ${FCFLAGS} -c -o $@ $<

$O/%.shared_module.o: ./src/shared/%.f90 ./include/constants.h
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.shared.o: ./src/shared/%.f90 $O/shared_par.shared_module.o
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.shared.o: ./src/shared/%.F90 $O/shared_par.shared_module.o
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.sharedmpi.o: ./src/shared/%.f90 $O/shared_par.shared_module.o $O/read_parameter_file.shared.o $O/read_value_parameters.shared.o
	${MPIFC} ${FCFLAGS} -c -o $@ $<

$O/%.cc.o: ./src/shared/%.c ./include/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

$O/%.gen.o: ./src/generate_databases/%.f90 $O/generate_databases_par.gen_mod.o
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.gen.o: ./src/generate_databases/%.F90 $O/generate_databases_par.gen_mod.o
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.gen_mod.o: ./src/generate_databases/%.F90 ./include/constants.h
	${FC} ${FCFLAGS} -c -o $@ $<

$O/get_gll_xyz.o: ./src/get_gll_xyz.f90 
	${FC} ${FCFLAGS} -c -o $@ $<

$O/calc_adepml_physical_jacobian.o: ./src/calc_adepml_physical_jacobian.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/create_mass_matrices_pml_elastic_loc.o: ./src/create_mass_matrices_pml_elastic_loc.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/cube2sph_trans.o: ./src/cube2sph_trans.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/calc_r_trans.o: ./src/calc_r_trans.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/setup_model_cartesian_subroutine.o: ./src/setup_model_cartesian_subroutine.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/node_stretching_subroutine.o: ./src/node_stretching_subroutine.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/find_parameter_value.o: ./src/find_parameter_value.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/update_parameters_from_netcdf.o: ./src/update_parameters_from_netcdf.f90
	${FC} ${FCFLAGS} ${NCDFCFLAGS} -c -o $@ $<
$O/netcdf_model_wrapper.o: ./src/netcdf_model_wrapper.f90
	${FC} ${FCFLAGS} ${NCDFCFLAGS} -c -o $@ $<

$O/netcdf_model_par.o: ./src/netcdf_model_par.f90
	${FC} ${FCFLAGS} ${NCDFCFLAGS} -c -o $@ $<

$O/get_boundary_from_crust_1_0.o: ./src/get_boundary_from_crust_1_0.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/get_boundary_from_netcdf.o: ./src/get_boundary_from_netcdf.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/get_model_mantle.o: ./src/get_model_mantle.f90
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.solver.o: ./src/specfem3D/%.f90 $O/shared_par.shared_module.o
	${FC} ${FCFLAGS} -c -o $@ $<
$O/%.solver.o: ./src/specfem3D/%.F90 $O/shared_par.shared_module.o
	${FC} ${FCFLAGS} -c -o $@ $<
$O/%.solverstatic_module.o: ./src/specfem3D/%.f90 ./OUTPUT_FILES/values_from_mesher.h $O/shared_par.shared_module.o
	${FC} ${FCFLAGS} -c -o $@ $<

$O/%.solverstatic_module.o: ./src/specfem3D/%.F90 ./OUTPUT_FILES/values_from_mesher.h $O/shared_par.shared_module.o
	${FC} ${FCFLAGS} -c -o $@ $<
$E/write_force_solution_file: ./src/write_force_solution_file.f90 $(write_force_solution_file_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/write_cmt_solution_file: ./src/write_cmt_solution_file.f90 $(write_cmt_solution_file_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/write_stations_file: ./src/write_stations_file.f90 $(write_stations_file_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/write_horizontal_slice_vtk: ./src/write_horizontal_slice_vtk.f90 $O/cube2sph_trans.o $(write_stations_file_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^
$E/write_vertical_slice_vtk: ./src/write_vertical_slice_vtk.f90 $O/cube2sph_trans.o $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/write_boundary_depth: ./src/write_boundary_depth.f90 $(write_boundary_obj) $(shared_obj) $(adios_shared_STUBS)
	${MPIFC} ${FCFLAGS} ${NCDFLDFLAGS} -o $@ $^ ${NCDFLDLIBS}

$E/write_model_tomo: ./src/write_model_tomo.f90 $(write_model_tomo_obj) $(shared_obj) $(adios_shared_STUBS)
	${MPIFC} ${FCFLAGS} ${NCDFLDFLAGS} -o $@ $^ ${NCDFLDLIBS}

#$O/debug.o: ./src/debug.f90
#	${FC} ${FCFLAGS} -c -o $@ $<

clean: 
	-rm -rf $E/* $O/*

.PHONY: all clean
