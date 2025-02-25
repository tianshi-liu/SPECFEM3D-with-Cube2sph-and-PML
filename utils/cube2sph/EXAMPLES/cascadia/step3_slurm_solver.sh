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

# cd $SLURM_SUBMIT_DIR
# module load intel openmpi
#module load vtune 

change_parfile() {
    local param=$1
    local value=$2 
    local file="DATA/Par_file"

    local oldstr=`grep "^$param " $file`
    local newstr="$param     =       $value"

    sed  "s?$oldstr?$newstr?g" $file  > tmp
    mv tmp $file 
}

currentdir=`pwd`

specfem_dir=${HOME}/code/c_fortran/SPECFEM3D-with-Cube2sph-and-PML/
#specfem_dir="${HOME}/src/specfem3d-cube2sph/"
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

change_parfile GPU_MODE .true.
change_parfile COUPLE_WITH_INJECTION_TECHNIQUE .true.
change_parfile NSTEP 1200


echo "start to run solver: `date`"
rm -f OUTPUT_FILES/timestamp* OUTPUT_FILES/*.semd
# runs solver
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "  running solver..."
  echo
   $specfem_dir/bin/xspecfem3D
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

python compare.py
code test.jpg
