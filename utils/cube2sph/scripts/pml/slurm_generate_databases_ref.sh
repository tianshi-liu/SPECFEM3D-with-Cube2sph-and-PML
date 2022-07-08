#!/usr/bin/bash

## job name and output file
#SBATCH --job-name specfem_mesher
#SBATCH --output %j.o

###########################################################
# USER PARAMETERS

## 40 CPUs ( 10*4 ), walltime 5 hour
#SBATCH --nodes=10
#SBATCH --ntasks=400
#SBATCH --time=00:15:00

###########################################################

cd $SLURM_SUBMIT_DIR
module load intel openmpi python/3.6.8

currentdir=`pwd`

specfem_dir="${HOME}/specfem3d_pml_1order"
#specfem_dir="${HOME}/specfem3d"
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

cp $specfem_dir/bin/xgenerate_databases bin

echo "start to run database generation: `date`"
# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  running database generation..."
  echo
  ./bin/xgenerate_databases
else
  # This is a MPI simulation
  echo
  echo "  running database generation on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xgenerate_databases
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo "database generation done: `date`"

echo
echo "    generate vtk file for vs... "
echo
./bin/xcombine_vol_data_vtk 0 399 vs DATABASES_MPI/ . 0
