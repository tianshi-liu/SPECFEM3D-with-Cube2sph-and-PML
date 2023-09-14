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

specfem_dir="${HOME}/specfem3d_pml_1order"
#specfem_dir="${HOME}/specfem3d"
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

cp ${specfem_dir}/bin/xspecfem3D bin
cp ${specfem_dir}/utils/cube2sph/bin/write_stations_file bin
cp ${specfem_dir}/utils/cube2sph/bin/write_force_solution_file bin
#cp DATA/FORCESOLUTION_sph DATA/FORCESOLUTION_sph
#cp DATA/STATIONS_sph DATA/STATIONS_sph
seis_dir=pml_prem
./bin/write_stations_file DATA/STATIONS_sph DATA/STATIONS DATA/rotate_nu .false. .false.
for comp in E; do
case $comp in
  Z)
    sed -i "/^component dir vect source E: /c\component dir vect source E:     0.d0" DATA/FORCESOLUTION_sph
    sed -i "/^component dir vect source N: /c\component dir vect source N:     0.d0" DATA/FORCESOLUTION_sph
    sed -i "/^component dir vect source Z_UP: /c\component dir vect source Z_UP:  1.d0" DATA/FORCESOLUTION_sph
  ;;
  N)
    sed -i "/^component dir vect source E: /c\component dir vect source E:     0.d0" DATA/FORCESOLUTION_sph
    sed -i "/^component dir vect source N: /c\component dir vect source N:     1.d0" DATA/FORCESOLUTION_sph
    sed -i "/^component dir vect source Z_UP: /c\component dir vect source Z_UP:  0.d0" DATA/FORCESOLUTION_sph 
  ;;
  E)
    sed -i "/^component dir vect source E: /c\component dir vect source E:     1.d0" DATA/FORCESOLUTION_sph
    sed -i "/^component dir vect source N: /c\component dir vect source N:     0.d0" DATA/FORCESOLUTION_sph
    sed -i "/^component dir vect source Z_UP: /c\component dir vect source Z_UP:  0.d0" DATA/FORCESOLUTION_sph
  ;;
esac
./bin/write_force_solution_file DATA/FORCESOLUTION_sph DATA/FORCESOLUTION .false. .false.
echo "  running solver for comp $comp on $NPROC processors..."
mpirun -np $NPROC ./bin/xspecfem3D
echo "  solver done for comp $comp"
mkdir -p ${seis_dir}/seis_sph_$comp
python3 rotate_seismogram.py --fn_matrix="DATA/rotate_nu" --rotate="XYZ->NEZ" --from_dir="OUTPUT_FILES" --to_dir="${seis_dir}/seis_sph_$comp" --from_template='${nt}.${sta}.BX${comp}.semd' --to_template='${nt}.${sta}.BX${comp}.sem.ascii'
cp DATA/FORCESOLUTION_sph DATA/STATIONS_sph DATA/FORCESOLUTION DATA/STATIONS ${seis_dir}/seis_sph_$comp
done

