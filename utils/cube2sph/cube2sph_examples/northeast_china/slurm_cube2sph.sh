#!/usr/bin/bash

## job name and output file
#SBATCH --job-name specfem_mesher
#SBATCH --output %j.o

###########################################################
# USER PARAMETERS

## 40 CPUs ( 10*4 ), walltime 5 hour
#SBATCH --nodes=5
#SBATCH --ntasks=200
#SBATCH --time=00:15:00

###########################################################

cd $SLURM_SUBMIT_DIR
module load intel openmpi python/3.6.8 hdf5/1.8.21 netcdf/4.6.3

currentdir=`pwd`

specfem_dir="${HOME}/specfem3d_pml_1order_repository"
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# generate database for the undeformed mesh, with a default model
echo -e ".true.\n.false." > adepml_stage
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

# start cube2sph transformation
echo -e ".false.\n.true." > adepml_stage

cp $specfem_dir/utils/cube2sph/bin/* bin/
cp -R DATABASES_MPI DATABASES_MPI_REF
#mpirun -np $NPROC ./bin/cube2sph_adepml

#./bin/cube2sph_force
#./bin/cube2sph_station
mpirun -np $NPROC ./bin/cube2sph_adepml DATABASES_MPI_REF DATABASES_MPI 43.5 121.8 0.0
cp DATABASES_MPI_REF/proc*adepml_param.bin DATABASES_MPI
cp $specfem_dir/bin/xgenerate_databases bin/
cp $specfem_dir/bin/xcombine_vol_data_vtk bin/
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

NPROC_MINUS_ONE=`echo "$NPROC-1" | bc`
./bin/xcombine_vol_data_vtk 0 ${NPROC_MINUS_ONE} vs DATABASES_MPI/ . 0


