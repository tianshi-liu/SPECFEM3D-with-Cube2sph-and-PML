#!/usr/bin/bash

## job name and output file
#SBATCH --job-name specfem_mesher
#SBATCH --output %j.o

###########################################################
# USER PARAMETERS

## 40 CPUs ( 10*4 ), walltime 5 hour
#SBATCH --nodes=10
#SBATCH --ntasks=400
#SBATCH --time=00:70:00

###########################################################

cd $SLURM_SUBMIT_DIR
module load intel openmpi python/3.6.8

currentdir=`pwd`

specfem_dir="${HOME}/specfem3d_pml_1order_repository"
#specfem_dir="${HOME}/specfem3d"
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

cp ${specfem_dir}/bin/xspecfem3D bin
cp ${specfem_dir}/utils/cube2sph/bin/write_stations_file bin
cp ${specfem_dir}/utils/cube2sph/bin/write_cmt_solution_file bin
cp ${specfem_dir}/utils/cube2sph/cube2sph_examples/rotate_seismogram.py .
#cp DATA/FORCESOLUTION_sph DATA/FORCESOLUTION_sph
#cp DATA/STATIONS_sph DATA/STATIONS_sph
ln -sf ${specfem_dir}/utils/cube2sph/DATA/topo_bathy DATA/
seis_dir=seis_rotate
./bin/write_stations_file DATA/STATIONS_sph DATA/STATIONS DATA/rotate_nu .true. .true.
./bin/write_cmt_solution_file DATA/CMTSOLUTION_sph DATA/CMTSOLUTION .true. .true.
echo "  running solver ..."
mpirun -np $NPROC ./bin/xspecfem3D
echo "  solver done"
mkdir -p ${seis_dir}
python3 rotate_seismogram.py --fn_matrix="DATA/rotate_nu" --rotate="XYZ->NEZ" --from_dir="OUTPUT_FILES" --to_dir="${seis_dir}" --from_template='${nt}.${sta}.BX${comp}.semd' --to_template='${nt}.${sta}.BX${comp}.sem.ascii'
cp DATA/CMTSOLUTION_sph DATA/STATIONS_sph DATA/CMTSOLUTION DATA/STATIONS ${seis_dir}
