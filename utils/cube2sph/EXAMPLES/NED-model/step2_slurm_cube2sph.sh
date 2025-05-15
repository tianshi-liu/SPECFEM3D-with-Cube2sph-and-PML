#!/usr/bin/bash

## job name and output file
#SBATCH --job-name specfem_mesher
#SBATCH --output step2_cube2sph.o

###########################################################
# USER PARAMETERS

## 40 CPUs ( 10*4 ), walltime 5 hour
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:18:00
#SBATCH --mem=0
#SBATCH --account=rrg-liuqy
#SBATCH --partition=debug
###########################################################
set -e

# set parameters
cube2sph_param="0. 0. 0.0"

#### stop here ###################

#cd $SLURM_SUBMIT_DIR
module load mpi hdf5/serial netcdf/serial

currentdir=`pwd`

#specfem_dir="${HOME}/SPECFEM3D-with-Cube2sph-and-PML"
specfem_dir=${HOME}/src/specfem3d-cube2sph/

# specfem_dr=${HOME}/scratch/test/RR/SPECFEM3D-with-Cube2sph-and-PML/
cub2sph_dir=$specfem_dir/utils/cube2sph/
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# generate database for the undeformed mesh, with a default model
echo -e ".true.\n.false." > adepml_stage
#cp $specfem_dir/bin/xgenerate_databases bin

echo "start to run database generation: `date`"
# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  running database generation..."
  echo
  $specfem_dir/bin/xgenerate_databases
else
  # This is a MPI simulation
  echo
  echo "  running database generation on $NPROC processors..."
  echo
  mpirun -np $NPROC $specfem_dir/bin/xgenerate_databases
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo "database generation done: `date`"
cp OUTPUT_FILES/output_generate_databases.txt .

# start cube2sph transformation
echo -e ".false.\n.true." > adepml_stage

#cp $specfem_dir/utils/cube2sph/bin/* bin/
cp -R DATABASES_MPI DATABASES_MPI_REF
#mpirun -np $NPROC ./bin/cube2sph_adepml

#./bin/cube2sph_force
#./bin/cube2sph_station
mpirun -np $NPROC $cub2sph_dir/bin/cube2sph_adepml DATABASES_MPI_REF DATABASES_MPI $cube2sph_param
cp DATABASES_MPI_REF/proc*adepml_param.bin DATABASES_MPI

#cp $specfem_dir/bin/xgenerate_databases bin/
#cp $specfem_dir/bin/xcombine_vol_data_vtk bin/
#cp DATA/Par_file_gll DATA/Par_file
sed -i "/^MODEL/c\MODEL                           = gll" DATA/Par_file
cp DATA/Par_file OUTPUT_FILES
echo "start to run database generation: `date`"
# runs database generation on isotropic model, so that smoothing can run properly.
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  running database generation..."
  echo
  $specfem_dir/bin/xgenerate_databases
else
  # This is a MPI simulation
  echo
  echo "  running database generation on $NPROC processors..."
  echo
  mpirun -np $NPROC $specfem_dir/bin/xgenerate_databases
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo "database generation done: `date`"

NPROC_MINUS_ONE=`echo "$NPROC-1" | bc`
$specfem_dir/bin/xcombine_vol_data_vtk 0 ${NPROC_MINUS_ONE} vs DATABASES_MPI/ . 0


