#!/usr/bin/bash
#./get_slice 400 vs ../optimize/SD_M01_slen0.14/DATABASES_MPI OUTPUT_FILES/slice_loc_45km.dat OUTPUT_FILES/slice_vs_M01_45km.val
#for id in 1 2 3 4
for id in 152
do
#./bin/xcreate_slice vs DATABASES_MPI DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_vs_init_vert${id}.val FALSE
#./bin/xcreate_slice beta_kernel_smooth_precon ../optimize/sum_kernels_M02/OUTPUT_SUM DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_kernel_M02_vert${id}.val
#./bin/xcreate_slice beta_kernel_smooth_precon ../optimize/sum_kernels_M08/OUTPUT_SUM DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_kernel_M08_vert${id}.val TRUE
./bin/xcreate_slice vs ../optimize/model_ckbd/DATABASES_MPI DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_vs_ckbd_vert${id}.val FALSE
./bin/xcreate_slice vs ../optimize/LBFGS_M24_slen0.01/DATABASES_MPI DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_vs_M24_vert${id}.val FALSE
#./bin/xcreate_slice dvsvs_smooth ../../subsample_test/specfem/vs_ckbd_intg DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_dvsvs_smooth_intg_ckbd_vert${id}.val TRUE
#./bin/xcreate_slice beta_kernel ../optimize/sum_kernels_Mckbd/OUTPUT_SUM_nia DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_beta_kernel_Mckbd_vert${id}.val FALSE
#./bin/xcreate_slice beta_kernel ../optimize/sum_kernels_Mfinal/OUTPUT_SUM DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_beta_kernel_final_vert${id}.val FALSE
#./bin/xcreate_slice dbeta ../optimize/sum_kernels_M23/INPUT_GRADIENT DATA/slice_vert${id}.xyz slices/slice_loc_vert${id}.dat slices/slice_dbeta_M23_vert${id}.val TRUE
done
