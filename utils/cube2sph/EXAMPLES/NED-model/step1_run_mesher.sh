#!/usr/bin/bash
#SBATCH --partition=debug 
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --output step1_mesher.o

# set param here 
myprocs=4
specfem_dir=${HOME}/src/specfem3d-cube2sph/
# specfem_dr=${HOME}/scratch/test/RR/SPECFEM3D-with-Cube2sph-and-PML/

## stop here ##############

cub2sph_dir=$specfem_dir/utils/cube2sph/
source module_load
echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo
#specfem_dir="${HOME}/specfem3d_pml"
#specfem_dir=${HOME}/SPECFEM3D-with-Cube2sph-and-PML/

# cleans output files
rm -rf DATABASES_MPI*
rm -rf MESH*
rm -rf OUTPUT_FILES_*
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*
mkdir -p MESH

# stores setup
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file_initmesh DATA/Par_file
cp DATA/Par_file OUTPUT_FILES/
#cp DATA/FORCESOLUTION OUTPUT_FILES/
#cp DATA/STATIONS OUTPUT_FILES/
if [ ! -e DATA/topo_bathy ]; then
  ln -sf $cub2sph_dir/DATA/topo_bathy DATA/
fi
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
  $specfem_dir/bin/xmeshfem3D
else
  # This is a MPI simulation
  echo
  echo "  running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC $specfem_dir//bin/xmeshfem3D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv OUTPUT_FILES/ OUTPUT_FILES_initmesh
cp -R MESH/ MESH-default
rm -rf DATABASES_MPI

mkdir -p OUTPUT_FILES
# stores setup
#cp DATA/Par_file_ref DATA/Par_file
sed -i "/^NPROC/c\NPROC                           = $myprocs" DATA/Par_file
sed -i "/^NGNOD/c\NGNOD                           = 27" DATA/Par_file

cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
#cp DATA/CMTSOLUTION OUTPUT_FILES/
#cp DATA/STATIONS OUTPUT_FILES/

cp MESH-default/* .
python hex8tohex27.py
bash change_name.sh
#python3 flip_num.py
cp nummaterial_velocity_file_tomo MESH-default/nummaterial_velocity_file

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# decomposes mesh using the pre-saved mesh files in MESH-default
echo
echo "  decomposing mesh..."
echo
$specfem_dir//bin/xdecompose_mesh $NPROC ./MESH-default $BASEMPIDIR
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
