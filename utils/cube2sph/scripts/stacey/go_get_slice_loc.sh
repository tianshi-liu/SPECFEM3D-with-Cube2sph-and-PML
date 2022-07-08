#!/usr/bin/bash

## job name and output file
#SBATCH --job-name specfem_mesher
#SBATCH --output %j.o

###########################################################
# USER PARAMETERS

## 40 CPUs ( 10*4 ), walltime 5 hour
#SBATCH --nodes=10
#SBATCH --ntasks=400
#SBATCH --time=00:30:00

###########################################################

cd $SLURM_SUBMIT_DIR

#module load gcc/7.3.0 openmpi/3.1.1

#module load intel/2018.3 intelmpi/2018.3

module load intel openmpi

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo "start: `date`"
# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  getting slice loc ..."
  echo
  ./get_slice_loc
else
  # This is a MPI simulation
  echo
  echo "  getting slice loc on $NPROC processors..."
  echo
  for depth in 10 15 20 25 30 35 40 45 50 60 70 80 90
  #for depth in 10 25 45 60
  do
  mpirun -np $NPROC ./bin/xcreate_slice_loc DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat
  #mpirun -np $NPROC ./bin/xcreate_slice_loc DATA/volume_0.1.xyz slices/volume_loc_0.1.dat
  done
  for lon in 144 152 160 1 2 3 4 5 E
  do
  mpirun -np $NPROC ./bin/xcreate_slice_loc DATA/slice_vert${lon}.xyz slices/slice_loc_vert${lon}.dat
  done
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo "done: `date`"
