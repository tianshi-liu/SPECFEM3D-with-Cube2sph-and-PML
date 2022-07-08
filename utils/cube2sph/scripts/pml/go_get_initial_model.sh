#!/usr/bin/bash

## job name and output file
#SBATCH --job-name combine_kernels
#SBATCH --output %j.o

###########################################################
# USER PARAMETERS

## 40 CPUs ( 10*4 ), walltime 5 hour
#SBATCH --nodes=10
#SBATCH --ntasks=400
#SBATCH --time=0:15:00

###########################################################
#model=M02
#model_prev=M01
#main_dir=$SCRATCH/ANAT_ani
old_model_dir=model_iso
new_model_dir=DATABASES_MPI
cd $SLURM_SUBMIT_DIR
module load intel openmpi
#mkdir -p ${main_dir}/mesh_$model/DATABASES_MPI
#ln -sf ~/specfem3d_subsample/bin/xupdate_model_by_perturbation bin/
#ln -sf ~/specfem3d_subsample/bin/xscale_rho_with_vs bin/
cp DATA/Par_file_sub_ani DATA/Par_file
ln -sf ~/specfem3d_subsample/bin/xget_isotropic_vbulk bin/
#ln -sf ~/specfem3d_subsample/bin/xwrite_zero_array bin/
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

#mpirun -np $NPROC ./bin/xscale_rho_with_vs ${old_model_dir} vsv vsh ckbd dvsv dvsh drho
#mpirun -np $NPROC ./bin/xupdate_model_by_perturbation dvsv,dvsh vsv,vsh TRUE,TRUE ${old_model_dir} ${new_model_dir} ckbd 0.03
mpirun -np $NPROC ./bin/xget_isotropic_vbulk ${old_model_dir} vp vs DATABASES_MPI vbulk
#mpirun -np $NPROC ./bin/xwrite_zero_array ${new_model_dir} Gc_nondim,Gs_nondim
#cd ${main_dir}/mesh_$model
#ln -sf ${main_dir}/specfem/DATA .
#mkdir bin OUTPUT_FILES
cp ${old_model_dir}/proc*rho.bin ${new_model_dir}
cp ${old_model_dir}/proc*vs.bin ${new_model_dir}
ls DATABASES_MPI/proc*vs.bin |
while read fn; do
  prefix=`echo $fn | cut -d'.' -f1`
  cp $fn ${prefix}v.bin
  cp $fn ${prefix}h.bin
done
ln -sf ~/specfem3d_subsample/bin/xgenerate_databases bin/
#cp ${main_dir}/specfem/OUTPUT_FILES/*.h OUTPUT_FILES
#ln -sf ${main_dir}/specfem/DATABASES_MPI/proc*_Database DATABASES_MPI/
#cp ${main_dir}/specfem/DATABASES_MPI/proc*_rho.bin DATABASES_MPI/
mpirun -np $NPROC ./bin/xgenerate_databases

