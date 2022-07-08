#!/usr/bin/bash
module load intel openmpi python/3.6.8
echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo
#specfem_dir="${HOME}/specfem3d_pml"
specfem_dir="${HOME}/specfem3d_pml_1order"
# cleans output files
rm -rf DATABASES_MPI*
rm -rf MESH*
rm -rf OUTPUT_FILES_*
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*
mkdir -p MESH

# links executables
#cd ../../../
#make clean
#make all
#cd $currentdir
mkdir -p bin
cd bin/
rm -f *
#cp $specfem_dir/bin/xmeshfem3D ./
#cp $specfem_dir/bin/xgenerate_databases ./
#cp $specfem_dir/bin/xspecfem3D ./
#cp $specfem_dir/bin/xdecompose_mesh ./
cp $specfem_dir/bin/* ./
cp $specfem_dir/utils/cube2sph/bin/* ./
cd ../

# stores setup
cp DATA/Par_file_initmesh DATA/Par_file
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# runs in-house mesher
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  running mesher..."
  echo
  ./bin/xmeshfem3D
else
  # This is a MPI simulation
  echo
  echo "  running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem3D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
mv OUTPUT_FILES/ OUTPUT_FILES_initmesh
cp -R MESH/ MESH-default
rm -rf DATABASES_MPI

mkdir -p OUTPUT_FILES
# stores setup
cp DATA/Par_file_ref DATA/Par_file
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

cp MESH-default/* .
python3 hex8tohex27.py
bash change_name.sh

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# decomposes mesh using the pre-saved mesh files in MESH-default
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC ./MESH-default $BASEMPIDIR
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


#echo "start to run database generation: `date`"
## runs database generation
#if [ "$NPROC" -eq 1 ]; then
#  # This is a serial simulation
#  echo
#  echo "  running database generation..."
#  echo
#  ./bin/xgenerate_databases
#else
#  # This is a MPI simulation
#  echo
#  echo "  running database generation on $NPROC processors..."
#  echo
#  mpirun -np $NPROC ./bin/xgenerate_databases
#fi
## checks exit code
#if [[ $? -ne 0 ]]; then exit 1; fi
#echo "database generation done: `date`"
#
##mpirun --hostfile hostfile -np $NPROC ./bin/setup_adepml
##cp DATA/FORCESOLUTION_cart DATA/FORCESOLUTION
##cp DATA/STATIONS_cart DATA/STATIONS
#
#cp -R DATABASES_MPI DATABASES_MPI_REF
#mpirun -np $NPROC ./bin/cube2sph_adepml
#./bin/cube2sph_force
#./bin/cube2sph_station
#cp DATABASES_MPI_REF/proc*adepml_param.bin DATABASES_MPI
