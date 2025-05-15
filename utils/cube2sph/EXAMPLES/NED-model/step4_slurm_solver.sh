#!/usr/bin/bash

## job name and output file
#SBATCH --job-name specfem_mesher
#SBATCH --output step3_solver.o

###########################################################
# USER PARAMETERS

## 40 CPUs ( 10*4 ), walltime 5 hour
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:15:00
#SBATCH --account=rrg-liuqy
#SBATCH --mem=0
#SBATCH --partition=debug
###########################################################

cd $SLURM_SUBMIT_DIR
module load intel openmpi
#module load vtune 

currentdir=`pwd`

specfem_dir=${HOME}//specfem3d-cube2sph/
#specfem_dir="${HOME}/SPECFEM3D-with-Cube2sph-and-PML"
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

echo "start to run solver: `date`"
# runs solver
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  running solver..."
  echo
  ./bin/xspecfem3D
else
  # This is a MPI simulation
  echo
  echo "  running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC $specfem_dir/bin/xspecfem3D
  #mpirun -np $NPROC vtune -collect hotspots -trace-mpi -r out.dir  ./bin/xspecfem3D
fi

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo "solver done: `date`"

