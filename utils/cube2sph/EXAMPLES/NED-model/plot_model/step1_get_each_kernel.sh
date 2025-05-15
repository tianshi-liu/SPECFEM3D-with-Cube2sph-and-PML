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
specfem_dir="${HOME}/scratch/specfem3d-cube2sph/"
cub2sph_dir=$specfem_dir/utils/cube2sph/
mkdir -p grdfolder pics profiles

# run index
run_indx=`seq 0 0`
#param_set="dbulk dbeta drho"
param_set="beta_kernel"
#param_set="vp vs rho"
evt_indx=`seq 21 29`

set -e 
# interpolate
module load intel 
#for param in vp vs rho;do 
for param in $param_set ;do 
for name in A;  do
  for iter in $run_indx;
  do 
      ii=`printf %02d $iter`
      for jj in $evt_indx;
      do
        $specfem_dir/bin/xcreate_slice $param  ../../solver/M$ii/P$jj/GRADIENT/ ../../DATABASES_MPI/  \
            input/prof$name.txt input/prof$name.loc $param.$name.$ii.evt$jj.out .false.
      done
  done

done
done

# generate_grd
module load gcc gmt-6.0.0
for param in $param_set;do 
for name in A;  do
  # get grd file and plot
  info=`gmt gmtinfo -C input/prof${name}.txt`
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
    for ievt in $evt_indx; do 
      awk '{print $5,-$4/1000}' input/prof$name.txt > tmp.1 
      info=`gmt gmtinfo -C $param.$name.$jj.evt$ievt.out`
      echo $info
      v1=`echo $info |awk '{print sqrt($1^2)}'`
      v2=`echo $info |awk '{print sqrt($2^2)}'`
      vmax=$v1
      if [[ $v2 > $v1 ]]; then 
        vmax=$v2 
      fi
      awk -v a=$vmax '{print $1/a}' $param.$name.$jj.evt$ievt.out > tmp.2 
      paste tmp.1 tmp.2 > tmp.3 

      # convert txt to grd
      gmt surface tmp.3 -Ggrdfolder/$param.iter$jj.evt$ievt.prof$name.grd -I$dx/$dz $bounds
      mv tmp.3 profiles/$param.iter$jj.evt$ievt.prof$name
      \rm tmp.1 tmp.2 $param.$name.$jj.evt$ievt.out 
      gmt grdinfo -C grdfolder/$param.iter$jj.evt$ievt.prof$name.grd
    done
  done

done
done