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

cd $SLURM_SUBMIT_DIR
module load intel openmpi
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
#mpirun -np $NPROC ./bin/xcombine_sem alpha_kernel,beta_kernel,rho_kernel,hess_kernel dir_list OUTPUT_SUM
#mpirun -np $NPROC ./bin/xcalc_kernel_abs_inv hess_kernel OUTPUT_SUM OUTPUT_SUM
#mpirun -np $NPROC ./bin/xcalc_kernel_abs_inv hess_kernel OUTPUT_SUM_set1 OUTPUT_SUM_set1
#for knm in alpha_kernel beta_kernel hess_kernel_abs; do
#for knm in alpha_kernel_precon beta_kernel_precon; do
#  mpirun -np $NPROC ./bin/xsmooth_sem_sph 50000 25000 $knm OUTPUT_SUM OUTPUT_SUM FALSE
#done
cp ~/specfem3d_subsample/bin/xsmooth_sem_sph_pde bin/
mpirun -np $NPROC ./bin/xsmooth_sem_sph_pde 30000 10000 betav_kernel ../event_kernel_Mani.set1/AK.ANM/OUTPUT_SUM ../event_kernel_Mani.set1/AK.ANM/OUTPUT_SUM FALSE
mpirun -np $NPROC ./bin/xsmooth_sem_sph_pde 30000 10000 betah_kernel ../event_kernel_Mani.set1/AK.ANM/OUTPUT_SUM ../event_kernel_Mani.set1/AK.ANM/OUTPUT_SUM FALSE
#mpirun -np $NPROC ./bin/xsmooth_sem_sph 30000 10000 dvsvs vs_ckbd_intg vs_ckbd_intg FALSE
#mpirun -np $NPROC ./bin/xsmooth_sem_sph 30000 10000 alpha_kernel,beta_kernel,rho_kernel,hess_kernel_abs_inv OUTPUT_SUM_set1 OUTPUT_SUM_set1 FALSE
#mpirun -np $NPROC ./bin/xcalc_kernel_div_hess alpha_kernel_smooth,beta_kernel_smooth,hess_kernel_abs_smooth OUTPUT_SUM OUTPUT_SUM
#mpirun -np $NPROC ./bin/xcalc_kernel_div_hess alpha_kernel,beta_kernel,hess_kernel OUTPUT_SUM OUTPUT_SUM
#for knm in alpha_kernel beta_kernel hess_kernel; do
#for knm in alpha_kernel_precon beta_kernel_precon; do
#  mpirun -np $NPROC ./bin/xsmooth_sem_sph 30000 15000 $knm OUTPUT_SUM OUTPUT_SUM FALSE
#done
#mpirun -np $NPROC ./bin/xcalc_kernel_div_hess alpha_kernel_smooth,beta_kernel_smooth,hess_kernel_smooth OUTPUT_SUM OUTPUT_SUM

