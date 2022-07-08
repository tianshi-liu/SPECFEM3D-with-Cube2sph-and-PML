CC=icc
FC=ifort
MPIFC=mpif90

FCFLAGS=-g -traceback -module ./obj -I./obj -I. -I./include
CFLAGS=-g -traceback -I./include

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
	$O/setup_model_cartesian_subroutine.o \
	$O/node_stretching_subroutine.o

cube2sph_station_obj = \
	$O/cube2sph_trans.o


## rules
all: $E/setup_adepml $E/cube2sph_adepml $E/cube2sph_force $E/cube2sph_station $E/write_horizontal_slice_vtk $E/write_vertical_slice_vtk $E/cube2sph_nopml

$E/setup_adepml: ./src/setup_adepml.f90 $(shared_obj) $(adios_shared_STUBS) $(generate_databases_obj) $(setup_adepml_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph_adepml: ./src/cube2sph_adepml.f90 $(setup_model_obj) $(node_stretching_obj) $(shared_obj) $(adios_shared_STUBS) $(generate_databases_obj) $(cube2sph_adepml_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph_nopml: ./src/cube2sph_nopml.f90 $(setup_model_obj) $(node_stretching_obj) $(shared_obj) $(adios_shared_STUBS) $(generate_databases_obj) $(cube2sph_adepml_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph_force: ./src/cube2sph_force.f90 $(cube2sph_station_obj) $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^

$E/cube2sph_station: ./src/cube2sph_station.f90 $(cube2sph_station_obj) $(shared_obj)
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

$E/write_horizontal_slice_vtk: ./src/write_horizontal_slice_vtk.f90 $O/cube2sph_trans.o $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^
$E/write_vertical_slice_vtk: ./src/write_vertical_slice_vtk.f90 $O/cube2sph_trans.o $(shared_obj)
	${MPIFC} ${FCFLAGS} -o $@ $^
#$O/debug.o: ./src/debug.f90
#	${FC} ${FCFLAGS} -c -o $@ $<

clean: 
	-rm -rf $E/* $O/*

.PHONY: all clean
