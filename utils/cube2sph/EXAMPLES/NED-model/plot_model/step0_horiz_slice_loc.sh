#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:12:59
#SBATCH --job-name step0
#SBATCH --output=step0.txt 
#SBATCH --mem=0
#SBATCH --partition=debug

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load NiaEnv/2019b
specfem_dir="${HOME}/specfem3d-cube2sph/"
cub2sph_dir=$specfem_dir/utils/cube2sph/
module load intel openmpi hdf5 netcdf 

mkdir -p grdfolder pics  input

# parameters
NPROC=8

# horizontal planes 
depth=("dafaf" 7.5 22.5) #in km
# profile
lon0=-1.2 lon1=1.2
#lon0=-124.25; lon1=-122
lat0=-1.2; lat1=1.2
#lat0=45.5; lat1=45.5
for i in `seq 1 2`;
do
  name=$i
  :>input/profB$name.txt 
  python src/generate_plane.py $lon0 $lon1 $lat0 $lat1 128 profile.txt
  dep=${depth[$i]}
  awk  '{print "XZ ADF",$2,$1,0.,a*1000}' a=$dep profile.txt > temp.txt
  $cub2sph_dir/bin/write_stations_file temp.txt temp1.txt rotation_nu .false. .false.

  # merge data
  awk  '{print $4,$3,$6/1000}' temp.txt >  profB$name.txt.temp
  awk '{print $4,$3,$6}' temp1.txt >  temp.txt 
  paste temp.txt profB$name.txt.temp > input/profB$name.txt
  #exit 

  mpirun -np $NPROC  $specfem_dir/bin/xcreate_slice_loc input/profB$name.txt ../tomo.ckbd/DATABASES_MPI/ temp.txt.loc
  \rm temp1.txt rotation_nu
  mv temp.txt.loc input/profB$name.loc
  \rm -f temp* profB$name.txt.temp
done 
