#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  3 . 0
#               ---------------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                              CNRS, France
#                       and Princeton University, USA
#                 (there are currently many more authors!)
#                           (c) October 2017
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================
#
# United States Government Sponsorship Acknowledged.
#
# Makefile.  Generated from Makefile.in by configure.
#######################################

FC = ifort
FCFLAGS = -g ${DEBUG_COUPLED_FLAG}
FC_DEFINE = -D
MPIFC = mpif90
MPILIBS = 

FLAGS_CHECK = -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -std03 -diag-disable 6477 -implicitnone -gen-interfaces -warn all -O3 -check nobounds -DFORCE_VECTORIZATION

FCFLAGS_f90 = -module ./obj -I./obj -I.  -I${SETUP}

FC_MODEXT = mod
FC_MODDIR = ./obj

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK) $(COND_OPENMP_FFLAGS)

MPIFCCOMPILE_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_CHECK) $(COND_OPENMP_FFLAGS)

CC = icc
CFLAGS = -g -O2 $(CPPFLAGS)
CPPFLAGS = -I${SETUP}  -DFORCE_VECTORIZATION $(COND_MPI_CPPFLAGS)

# all linker flags
MPILIBS +=  

#######################################
####
#### MPI
####
#######################################

## serial or parallel
MPI = yes
#MPI = no

FCLINK = $(MPIFCCOMPILE_CHECK)
#FCLINK = $(FCCOMPILE_CHECK)

COND_MPI_CPPFLAGS = $(FC_DEFINE)WITH_MPI
#COND_MPI_CPPFLAGS =

# objects toggled between the parallel and serial version
COND_MPI_OBJECTS = $O/parallel.sharedmpi.o
#COND_MPI_OBJECTS = $O/serial.shared.o

MPI_INCLUDES =  -I/scinet/niagara/software/2019b/opt/intel-2019u4/openmpi/4.0.1-hpcx2.5/include

#######################################
####
#### SCOTCH
####
#######################################

SCOTCH = yes
#SCOTCH = no

USE_BUNDLED_SCOTCH = 1

SCOTCH_DIR = ./external_libs/scotch
SCOTCH_INCDIR = ./external_libs/scotch/include
SCOTCH_LIBDIR = ./external_libs/scotch/lib

SCOTCH_INC = -I${SCOTCH_INCDIR}
#SCOTCH_INC =

SCOTCH_LIBS = -L${SCOTCH_LIBDIR} -lscotch -lscotcherr
#SCOTCH_LIBS =

SCOTCH_FLAGS = $(SCOTCH_INC) $(FC_DEFINE)USE_SCOTCH
#SCOTCH_FLAGS =

### added support for METIS as well, thus uncomment the line below and compile METIS if you want to use it instead of SCOTCH
# add $(FC_DEFINE)USE_METIS_INSTEAD_OF_SCOTCH
#SCOTCH_LIBS = -L./external_libs/metis-4.0.3 -L./metis-4.0.3 -lmetis
#SCOTCH_FLAGS = $(FC_DEFINE)USE_METIS_INSTEAD_OF_SCOTCH

#######################################
####
#### CUDA
#### with configure: ./configure --with-cuda=cuda5 CUDA_FLAGS=.. CUDA_LIB=.. CUDA_INC=.. MPI_INC=.. ..
####
#######################################

#CUDA = yes
CUDA = no

#CUDA5 = yes
CUDA5 = no

#CUDA6 = yes
CUDA6 = no

#CUDA7 = yes
CUDA7 = no

#CUDA8 = yes
CUDA8 = no

#CUDA9 = yes
CUDA9 = no

# CUDA compilation with linking
#CUDA_PLUS = yes
CUDA_PLUS = no

# default cuda libraries
# runtime library -lcudart needed, others are optional -lcuda -lcublas

CUDA_FLAGS = 
CUDA_INC =  -I${SETUP}
CUDA_LINK =   -lstdc++

#NVCC = nvcc
NVCC = icc

##
## GPU architecture
##
# CUDA architecture / code version
# Fermi: -gencode=arch=compute_10,code=sm_10 not supported
# Tesla (default): -gencode=arch=compute_20,code=sm_20
# Geforce GT 650m: -gencode=arch=compute_30,code=sm_30
# Kepler (cuda5,K20) : -gencode=arch=compute_35,code=sm_35
# Kepler (cuda6.5,K80): -gencode=arch=compute_37,code=sm_37
# Maxwell (cuda7,K2200): -gencode=arch=compute_50,code=sm_50
# Pascal (cuda8,P100): -gencode=arch=compute_60,code=sm_60
# Volta (cuda9,V100): -gencode=arch=compute_70,code=sm_70
GENCODE_20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE_30 = -gencode=arch=compute_30,code=\"sm_30,compute_30\"
GENCODE_35 = -gencode=arch=compute_35,code=\"sm_35,compute_35\"
GENCODE_37 = -gencode=arch=compute_37,code=\"sm_37\"
GENCODE_50 = -gencode=arch=compute_50,code=\"sm_50,compute_50\"
GENCODE_60 = -gencode=arch=compute_60,code=\"sm_60,compute_60\"
GENCODE_70 = -gencode=arch=compute_70,code=\"sm_70,compute_70\"

# cuda preprocessor flag
# CUDA version 9.0
##GENCODE = $(GENCODE_70) $(FC_DEFINE)GPU_DEVICE_Volta
# CUDA version 8.0
##GENCODE = $(GENCODE_60) $(FC_DEFINE)GPU_DEVICE_Pascal
# CUDA version 7.x
##GENCODE = $(GENCODE_50) $(FC_DEFINE)GPU_DEVICE_Maxwell
# CUDA version 6.5
##GENCODE = $(GENCODE_37) $(FC_DEFINE)GPU_DEVICE_K80
# CUDA version 5.x
##GENCODE = $(GENCODE_35) $(FC_DEFINE)GPU_DEVICE_K20
# CUDA version 4.x
#GENCODE = $(GENCODE_20)

# CUDA flags and linking
#NVCC_FLAGS_BASE = $(CUDA_FLAGS) $(CUDA_INC) $(MPI_INCLUDES) $(COND_MPI_CPPFLAGS)
##NVCC_FLAGS = $(NVCC_FLAGS_BASE) -dc -DCUDA $(GENCODE)
#NVCC_FLAGS = $(NVCC_FLAGS_BASE) -DCUDA -DUSE_OLDER_CUDA4_GPU $(GENCODE)

##NVCCLINK_BASE = $(NVCC) $(CUDA_INC) $(MPI_INCLUDES) $(COND_MPI_CPPFLAGS) -DCUDA
##NVCCLINK = $(NVCCLINK_BASE) -dlink $(GENCODE)
#NVCCLINK = $(NVCCLINK_BASE) -DUSE_OLDER_CUDA4_GPU $(GENCODE)

NVCC_FLAGS = $(MPI_INCLUDES) $(COND_MPI_CPPFLAGS)
NVCCLINK = $(NVCC) $(NVCC_FLAGS)


#######################################
####
#### OpenMP
#### with configure: ./configure --enable-openmp OMP_FCFLAGS=".." OMP_LIB=..
####
#######################################

#OPENMP = yes
OPENMP = no

#OMP_FLAGS = -qopenmp $(FC_DEFINE)USE_OPENMP
OMP_FLAGS =

FCFLAGS += $(OMP_FLAGS)

#OMP_LIBS = $(OMP_LIB)
OMP_LIBS =

# objects toggled between openmp and non-openmp version
##COND_OMP_OBJECTS = $O/older_not_maintained_compute_forces_viscoelastic_Dev_openmp.openmp.o
#COND_OMP_OBJECTS =
COND_OMP_OBJECTS =

#######################################
####
#### VTK
#### with configure: ./configure --enable-vtk --with-vtk-version=5.8 ..
####
#######################################

#VTK = yes
VTK = no

VTK_MAJOR_VERSION = 

# additional libraries
ifeq ($(VTK),yes)
  ifeq ($(shell test $(VTK_MAJOR_VERSION) -gt 5; echo $$?),0)
    VTKLIBS = -lvtkRenderingOpenGL2-7.0
  endif
endif

#FCCOMPILE_CHECK += $(FC_DEFINE)VTK_VIS
#CPPFLAGS += 
#VTKLIBS +=  

#######################################
####
#### ADIOS
#### with configure: ./configure --with-adios ADIOS_CONFIG=..
####
#######################################

#ADIOS = yes
ADIOS = no

ADIOS_INC = 
ADIOS_LINK = 

FCFLAGS_f90 += $(ADIOS_INC)
MPILIBS += $(ADIOS_LINK)

#######################################
####
#### ASDF
#### with configure: ./configure --with-asdf ASDF_LIBS=..
####
#######################################

#ASDF = yes
ASDF = no

#FCFLAGS += @ASDF_FCFLAGS@
MPILIBS += 

#######################################
####
#### directories
####
#######################################

## compilation directories
# B : build directory
B = .
# E : executables directory
E = $B/bin
# O : objects directory
O = $B/obj
# S_TOP : source file root directory
S_TOP = .
# L : libraries directory
L = $B/lib
# setup file directory
SETUP = $B/setup
# output file directory
OUTPUT = $B/OUTPUT_FILES


#######################################
####
#### targets
####
#######################################

# code subdirectories
SUBDIRS = \
	auxiliaries \
	check_mesh_quality \
	cuda \
	decompose_mesh \
	generate_databases \
	meshfem3D \
	shared \
	specfem3D \
	tomography/postprocess_sensitivity_kernels \
	tomography \
	$(EMPTY_MACRO)

# default targets for the pure Fortran version
DEFAULT = \
	xdecompose_mesh \
	xmeshfem3D \
	xgenerate_databases \
	xspecfem3D \
	$(EMPTY_MACRO)

ifeq ($(MPI),yes)
DEFAULT += xdecompose_mesh_mpi
DEFAULT += xinverse_problem_for_model
SUBDIRS += inverse_problem_for_model
endif

default: $(DEFAULT)

all: default aux check_mesh postprocess tomography

backup:
	cp -rp src setup DATA/Par_file* Makefile bak

ifdef CLEAN
clean:
	@echo "cleaning by CLEAN"
	-rm -f $(foreach dir, $(CLEAN), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_SHARED_OBJECTS) $($(dir)_TARGETS))
	-rm -f ${E}/*__genmod.*
	-rm -f ${O}/*__genmod.*
else
clean:
	@echo "cleaning all"
	-rm -f $(foreach dir, $(SUBDIRS), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_TARGETS))
	-rm -f ${E}/*__genmod.*
	-rm -f ${O}/*__genmod.*
endif

realclean: clean
ifeq (${USE_BUNDLED_SCOTCH},1)
	@echo "cleaning bundled Scotch in directory: ${SCOTCH_DIR}/src"
	$(MAKE) -C ${SCOTCH_DIR}/src realclean
endif
	-rm -rf $E/* $O/*

# unit testing
# If the first argument is "test"...
ifeq (test,$(findstring test,firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  TEST_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # turn them into do-nothing targets
  $(eval $(TEST_ARGS):;@:)
endif

tests:
	@echo "testing in directory: ${S_TOP}/tests/"
	cd ${S_TOP}/tests; ./run_all_tests.sh $(TEST_ARGS)
	@echo ""

help:
	@echo "usage: make [executable]"
	@echo ""
	@echo "supported main executables:"
	@echo "    xdecompose_mesh"
	@echo "    xmeshfem3D"
	@echo "    xgenerate_databases"
	@echo "    xspecfem3D"
	@echo "    xinverse_problem_for_model"
	@echo ""
	@echo "defaults:"
	@echo "    xdecompose_mesh"
	@echo "    xmeshfem3D"
	@echo "    xgenerate_databases"
	@echo "    xspecfem3D"
	@echo "    xinverse_problem_for_model"
	@echo ""
	@echo "    xcombine_surf_data"
	@echo "    xcombine_vol_data"
	@echo "    xcombine_vol_data_vtk"
	@echo "    xconvolve_source_timefunction"
	@echo "    xcreate_movie_shakemap_AVS_DX_GMT"
	@echo ""
	@echo "    xcheck_mesh_quality"
	@echo "    xconvert_skewness_to_angle"
	@echo ""
	@echo "additional executables:"
	@echo "- auxiliary executables: [make aux]"
	@echo "    xcombine_surf_data"
	@echo "    xcombine_vol_data"
	@echo "    xcombine_vol_data_vtk"
	@echo "    xconvolve_source_timefunction"
	@echo "    xcreate_movie_shakemap_AVS_DX_GMT"
	@echo ""
	@echo "- check mesh executables: [make check_mesh]"
	@echo "    xcheck_mesh_quality"
	@echo "    xconvert_skewness_to_angle"
	@echo ""
	@echo "- sensitivity kernel postprocessing tools: [make postprocess]"
	@echo "    xclip_sem"
	@echo "    xcombine_sem"
	@echo "    xsmooth_sem"
	@echo ""
	@echo "- tomography tools: [make tomography]"
	@echo "    xmodel_update"
	@echo "    xsum_kernels"
	@echo "    xsum_preconditioned_kernels"
	@echo ""
	@echo "for unit testing:"
	@echo "    tests"
	@echo ""

.PHONY: all default backup clean realclean help tests

#######################################

${SETUP}/version.fh: ./.git/logs/HEAD
	@echo "GEN $@"
	@echo "! This file is generated by Make. Do not edit this file!" > $@
	@echo "character(len=*), parameter :: git_package_version = \"$$(cd ${S_TOP} && git describe --tags)\"" >> $@
	@echo "character(len=*), parameter :: git_commit_version = \"$$(cd ${S_TOP} && git log -n 1 | head -1 | /bin/grep commit)\"" >> $@
	@echo "character(len=*), parameter :: git_date_version = \"From $$(cd ${S_TOP} && git log -n 1 | head -5 | /bin/grep Date:)\"" >> $@

#######################################

# Get dependencies and rules for building stuff
include $(patsubst %, ${S_TOP}/src/%/rules.mk, $(SUBDIRS))

#######################################

##
## Shortcuts
##

# Shortcut for: <prog>/<xprog> -> bin/<xprog>
define target_shortcut
$(patsubst $E/%, %, $(1)): $(1)
.PHONY: $(patsubst $E/%, %, $(1))
$(patsubst $E/x%, %, $(1)): $(1)
.PHONY: $(patsubst $E/x%, %, $(1))
endef

# Shortcut for: dir -> src/dir/<targets in here>
define shortcut
$(1): $($(1)_TARGETS)
.PHONY: $(1)
$$(foreach target, $$(filter $E/%,$$($(1)_TARGETS)), $$(eval $$(call target_shortcut,$$(target))))
endef

$(foreach dir, $(SUBDIRS), $(eval $(call shortcut,$(dir))))

# testing
test : tests

# Other old shortcuts
bak: backup
mesh: $E/xmeshfem3D
gen: $E/xgenerate_databases
spec: $E/xspecfem3D
dec: $E/xdecompose_mesh

.PHONY: bak mesh gen spec dec

