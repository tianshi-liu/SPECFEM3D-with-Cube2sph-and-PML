#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=160
#SBATCH --time=00:25:59
#SBATCH --job-name=POST
#SBATCH --output=POST_%j.txt
#SBATCH --partition=compute 

set -e 
# include file
. parameters.sh
# load modules 
source module_env

# input vars
MODEL=M24
SOURCE_FILE=./src_rec/sources.dat.tele
NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`

echo "running fwat2 " > output_fwat2_log_${MODEL}.txt

# parfile changer script
change_par=$FWATLIB/change_par_file.sh

# get search direction
MODEL_START=`cat lbfgs.in |head -1`
iter=`echo "$MODEL" |awk -F'M' '{print $2}' | awk '{printf "%d",$1}'`
iter_start=`echo "$MODEL_START" |awk -F'M' '{print $2}' | awk '{printf "%d",$1}'` 

# copy model to optimize
LOCAL_PATH=./DATABASES_MPI
if [ $MODEL == "M00"  ];then
	rm -rf optimize/MODEL_M00; mkdir -p optimize/MODEL_M00
	\cp -r $LOCAL_PATH/* optimize/MODEL_M00
fi

# create log file
filename=output_fwat2_log_${MODEL}.txt
:> $filename

# sum kernels
PRECOND=`cat fwat_params/FWAT.PAR |grep '^PRECOND_TYPE' |awk -F: '{print $2}'`
mpirun -np $NPROC python $OPT_LIB/sum_kernel.py  $iter $SOURCE_FILE $PRECOND >> $filename

# set smoothing parameters
$change_par LOCAL_PATH $LOCAL_PATH ./DATA/Par_file
$change_par LOCAL_PATH $LOCAL_PATH ./DATA/meshfem3D_files/Mesh_Par_file
sigma_h=`grep '^SIGMA_H:' fwat_params/FWAT.PAR |awk -F: '{print $2}'`
sigma_v=`grep '^SIGMA_V:' fwat_params/FWAT.PAR |awk -F: '{print $2}'`

# smooth hess kernel if required
if [ $PRECOND == "default" ] && [ $MODEL == "M00"  ];then 
  param=hess_kernel
  mv optimize/SUM_KERNELS_$MODEL/*_$param.bin $LOCAL_PATH
  mpirun -np $NPROC $fksem/bin/xsmooth_sem_sph_pde 50000 25000 $param $LOCAL_PATH optimize/SUM_KERNELS_$MODEL/ .false.
  \rm $LOCAL_PATH/*_${param}.bin
  for i in `seq 1 $NPROC`;
  do
    ii=`echo $i |awk '{printf "%06d", $1-1}'`
    name=optimize/SUM_KERNELS_$MODEL/proc${ii}_$param
    mv ${name}_smooth.bin $name.bin 
  done 
fi 

# get search direction
mpirun -np $NPROC python $OPT_LIB/get_lbfgs_direc.py $iter $iter_start

# smooth search direction
kl_list=`cd $OPT_LIB; python -c "from tools import *; mod = get_direc_name_list(); print(' '.join(mod))"`
for param in $kl_list; 
do 
  mv optimize/SUM_KERNELS_$MODEL/*_$param.bin $LOCAL_PATH
  mpirun -np $NPROC $fksem/bin/xsmooth_sem_sph_pde $sigma_h $sigma_v $param $LOCAL_PATH optimize/SUM_KERNELS_$MODEL/ .false.
  \rm $LOCAL_PATH/*_$param.bin
  for i in `seq 1 $NPROC`;
  do
    ii=`echo $i |awk '{printf "%06d", $1-1}'`
    name=optimize/SUM_KERNELS_$MODEL/proc${ii}_$param
    mv ${name}_smooth.bin $name.bin 
  done 
done

# generate new model
LSDIR=./optimize/MODEL_${MODEL}_step01
mkdir -p $LSDIR
step_fac=`tail -1 lbfgs.in`
MIS_FILE=plots/$MODEL.mis
echo "python $OPT_LIB/get_lbfgs_step_fac.py $MODEL $LSDIR $step_fac $NPROC"
python $OPT_LIB/get_lbfgs_step_fac.py $MODEL $LSDIR $step_fac $NPROC  > $MIS_FILE
step_fac=`tail -1 $MIS_FILE |cut -d'=' -f2 |awk '{print $1}'`
dmax=`tail -1 $MIS_FILE |cut -d'=' -f2 |awk '{print $2}'`

# generate new model database
$change_par LOCAL_PATH $LSDIR DATA/Par_file
$change_par LOCAL_PATH $LSDIR  DATA/meshfem3D_files/Mesh_Par_file
$change_par SAVE_MESH_FILES .false. DATA/Par_file
# copy info to new 
\cp DATA/adepml_* .
\cp DATA/wavefield* .
\cp  $LOCAL_PATH/*Database $LSDIR/
\cp  $LOCAL_PATH/*adepml* $LSDIR/
\cp  $LOCAL_PATH/*undeformed_xyz.bin $LSDIR/
mpirun -np $NPROC $fksem/bin/xgenerate_databases 

# delete 
\rm adepml_* wavefield* 
echo " " >> $filename
echo "******************************************************" >> $filename
echo " Finished FWAT stage2 here!!!" >> $filename 
echo " " >> $filename
