#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --time=00:25:59
#SBATCH --job-name=smooth
#SBATCH --output=step3_smooth.o
#SBATCH --partition=debug

set -e 

# sem path
specfem_dir=${HOME}//specfem3d-cube2sph/

# load modules 
module load intel openmpi hdf5 netcdf

LOCAL_PATH=./DATABASES_MPI
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

# smooth vp/vs/rho if required
mkdir -p savedir
rm -rf savedir/*
for param in vp vs rho; 
do
  mv $LOCAL_PATH/*_${param}.bin savedir/
  mpirun -np $NPROC $specfem_dir/bin/xsmooth_sem_sph_pde 8000 8000 $param savedir savedir .false.
  for i in `seq 1 $NPROC`;
  do
    ii=`echo $i |awk '{printf "%06d", $1-1}'`
    name=savedir/proc${ii}_$param
    mv ${name}_smooth.bin $name.bin 
  done 
done

# copy to local_path
\cp savedir/* $LOCAL_PATH/
rm -rf savedir

# regenerate database
echo "rerun database generation: `date`"
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

