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

# profile
lon0=-1.2 lon1=1.2
#lon0=-124.25; lon1=-122
lat0=0.25; lat1=0.25
#lat0=45.5; lat1=45.5
z0=-30
z1=0.
dist=5
nz=101

# gmt grdcut @earth_relief_03s -R-125/-119.5/44/45  -Gout.grd
# horizontal line
python src/generate_line.py $lon0 $lon1 $lat0 $lat1 300 profile.txt

for name in A; do 
  :>input/prof$name.txt
  for ((i=0;i<$nz;i++));
  do
    dep=`printf %g $(echo "scale=4; 1000*($z0 + ($z1 - $z0) / ($nz-1) * $i)"|bc)`
    awk  '{print "XZ ADF",$1,$2,a,$3}' a=$dep profile.txt >> prof$name.txt.temp
  done

  # get grd file and plot
  awk '{print $1,$2,$4,$3,0.0,-$5}' prof$name.txt.temp > temp.txt
  $cub2sph_dir/bin/write_stations_file temp.txt temp1.txt rotation_nu .false. .false.
  awk '{print $4,$3,$6}' temp1.txt > temp2.txt
  awk  '{print $5,$6}' prof$name.txt.temp > temp3.txt 
  paste temp2.txt temp3.txt > input/prof$name.txt

  mpirun -np $NPROC  $specfem_dir/bin/xcreate_slice_loc input/prof$name.txt ../tomo.ckbd/DATABASES_MPI/ temp.txt.loc

  \rm temp1.txt rotation_nu
  mv temp.txt.loc input/prof$name.loc
  \rm temp* prof$name.txt.temp

done 
