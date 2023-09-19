#!/usr/bin/bash
module load intel openmpi python/3.6.8
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up a cube-shaped mesh..."
echo
### CHANGE THIS DIRECTORY
specfem_dir="${HOME}/specfem3d_pml_1order_repository"
##########################
# cleans output files
rm -rf DATABASES_MPI*
rm -rf MESH*
rm -rf OUTPUT_FILES_*
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*
mkdir -p MESH

mkdir -p bin
cd bin/
rm -f *
cp $specfem_dir/bin/* ./
cp $specfem_dir/utils/cube2sph/bin/* ./
cd ../
# setup parameter files
cp DATA/Par_file_initmesh DATA/Par_file
# stores setup
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/

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
  echo "  mesher must be run in serial for cube2sph mesh setup"
  echo
  exit 1
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
mv OUTPUT_FILES/ OUTPUT_FILES_initmesh
cp -R MESH/ MESH-default
rm -rf DATABASES_MPI

mkdir -p OUTPUT_FILES
# change parameter files for partitioning
sed -i "/^NPROC/c\NPROC                           = 400" DATA/Par_file
sed -i "/^NGNOD/c\NGNOD                           = 27" DATA/Par_file

cp DATA/Par_file OUTPUT_FILES/

cp MESH-default/* .
python3 hex8tohex27.py
bash change_name.sh
cp nummaterial_velocity_file_tomo MESH-default/nummaterial_velocity_file

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


