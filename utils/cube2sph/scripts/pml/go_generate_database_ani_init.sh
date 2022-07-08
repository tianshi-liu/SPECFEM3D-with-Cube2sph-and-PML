#!/usr/bin/bash

## job name and output file
#SBATCH --job-name specfem_mesher
#SBATCH --output %j.o

###########################################################
# USER PARAMETERS

## 40 CPUs ( 10*4 ), walltime 5 hour
#SBATCH --nodes=10
#SBATCH --ntasks=400
#SBATCH --time=00:20:00

###########################################################

cd $SLURM_SUBMIT_DIR

#module load gcc/7.3.0 openmpi/3.1.1

#module load intel/2018.3 intelmpi/2018.3

module load intel openmpi

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

#cp ../src_rec_sub/cartesian/FORCESOLUTION_AK.ANM_X DATA/FORCESOLUTION
cp DATA/Par_file_sub_ani DATA/Par_file
#cp -r ../model_M13 ../model_M13_ani
cp ~/specfem3d_subsample/bin/xwrite_zero_array bin/
mpirun -np $NPROC ./bin/xwrite_zero_array $SCRATCH/ANAT_ani/model_M13_ani/DATABASES_MPI Gc_nondim,Gs_nondim
cd ../model_M13_ani
ln -sf ~/specfem3d_subsample/bin/xgenerate_databases bin/
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
