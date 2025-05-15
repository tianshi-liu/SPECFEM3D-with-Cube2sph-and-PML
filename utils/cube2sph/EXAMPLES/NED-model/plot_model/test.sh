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
mkdir -p grdfolder pics profiles

# run index
run_indx=`seq 29 29`
#param_set="dbulk dbeta drho"
#param_set="hess_kernel"
param_set="vp vs rho"

set -e 

# generate_grd
module load gcc gmt-6.0.0
for param in $param_set;do 
for name in B1 B2 B3 B4 B5;  do
  # get grd file and plot
  info=`gmt gmtinfo -C input/prof${name}.txt`
  xmin=`echo $info | awk '{print $7}'`
  xmax=`echo $info | awk '{print $8}'`
  ymin=`echo $info | awk '{print $9}'`
  ymax=`echo $info | awk '{print $10}'`

  echo $xmin $xmax $ymin $ymax 
  exit 

  dmax=`echo $info | awk '{print $10}'`
  z0=`echo $info | awk '{print -$8/1000}'`
  z1=`echo $info | awk '{print -$7/1000}'`
  bounds=-R0/$dmax/$z0/$z1 
  proj=-JX12c/6c
  dx=`echo "0. $dmax" | awk '{print (-$1 + $2) / 128.}'`
  dz=`echo "$z0 $z1" | awk '{print (-$1 + $2)/128.}'`

  for iter in $run_indx;
  do 
      jj=`printf %02d $iter`
      echo $jj 

      awk '{print $5,-$4/1000}' input/prof$name.txt > tmp.1 
      info=`gmt gmtinfo -C $param.$name.$jj.out`
      echo $info
      vmax=`echo $info |awk '{print $2}'`
      if [[ `echo $vmax 1 |awk '{print $1<=$2}'` == 1 ]]; then 
        awk -v a=$vmax '{print $1/a}' $param.$name.$jj.out > tmp.2 
      else
        awk '{print $1}' $param.$name.$jj.out > tmp.2 
      fi
      paste tmp.1 tmp.2 > tmp.3 

      # convert txt to grd
      gmt surface tmp.3 -Ggrdfolder/$param.iter$jj.prof$name.grd -I$dx/$dz $bounds
      mv tmp.3 profiles/$param.iter$jj.prof$name
      \rm tmp.1 tmp.2 $param.$name.$jj.out 
      gmt grdinfo -C grdfolder/$param.iter$jj.prof$name.grd
  done

done
done