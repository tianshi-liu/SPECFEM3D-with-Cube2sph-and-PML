#!/usr/bin/bash
#./get_slice 400 vs ../optimize/SD_M01_slen0.14/DATABASES_MPI OUTPUT_FILES/slice_loc_45km.dat OUTPUT_FILES/slice_vs_M01_45km.val
#for depth in 10 15 20 25 30 35 40 45 50
#for depth in 60 70 80 90
#for depth in 10 25 45 60
#model=M19
#for depth in 10 25 45 60
for depth in 10 15 20 25 30 35 40 45 50 60 70 80 90
do
./bin/xcreate_slice vsv DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_vsv_init_${depth}km_0.02.val FALSE
./bin/xcreate_slice vs_unsmooth DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_vs_unsmooth_init_${depth}km_0.02.val FALSE
done
for lon in 144 152 160 1 2 3 4 5 E
do
./bin/xcreate_slice vsv DATABASES_MPI DATA/slice_vert${lon}.xyz slices/slice_loc_vert${lon}.dat slices/slice_vsv_init_${lon}.val FALSE
./bin/xcreate_slice vs_unsmooth DATABASES_MPI DATA/slice_vert${lon}.xyz slices/slice_loc_vert${lon}.dat slices/slice_vs_unsmooth_init_${lon}.val FALSE
done
#./bin/xcreate_slice vs DATABASES_MPI DATA/slice_${depth}km_0.02.xyz OUTPUT_FILES/slice_loc_${depth}km_0.02.dat OUTPUT_FILES/slice_vs_init_${depth}km_0.02.val
#./bin/xcreate_slice beta_kernel_smooth ../optimize/sum_kernels_M05/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_beta_kernel_smooth_M05_${depth}km_0.02.val TRUE
#./bin/xcreate_slice dbeta ../optimize/sum_kernels_Mckbd/INPUT_GRADIENT DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_dbeta_Mckbd_${depth}km_0.02.val TRUE
#./bin/xcreate_slice vs ../optimize/LBFGS_M24_slen0.01/DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_vs_M24_${depth}km_0.02.val FALSE
#./bin/xcreate_slice beta_kernel ../optimize/sum_kernels_Mckbd/OUTPUT_SUM_nia DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_beta_kernel_Mckbd_${depth}km_0.02.val FALSE


#./bin/xcreate_slice dbetav ../event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_dbetav_${model}_${depth}km_0.02.val TRUE
#./bin/xcreate_slice dbetah ../event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_dbetah_${model}_${depth}km_0.02.val TRUE
#./bin/xcreate_slice dGc_nondim ../event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_dGc_nondim_${model}_${depth}km_0.02.val TRUE
#./bin/xcreate_slice dGs_nondim ../event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_dGs_nondim_${model}_${depth}km_0.02.val TRUE
#./bin/xcreate_slice betav_kernel_smooth_norm ../event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_betav_kernel_${model}_${depth}km_0.02.val TRUE
#./bin/xcreate_slice betah_kernel_smooth_norm ../event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_betah_kernel_${model}_${depth}km_0.02.val TRUE

#################
#./bin/xcreate_slice vsv ../model_${model}/DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_vsv_${model}_${depth}km_0.02.val FALSE
#./bin/xcreate_slice vsh ../model_${model}/DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_vsh_${model}_${depth}km_0.02.val FALSE
#./bin/xcreate_slice Gc_nondim ../model_${model}/DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_Gc_nondim_${model}_${depth}km_0.02.val FALSE
#./bin/xcreate_slice Gs_nondim ../model_${model}/DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_Gs_nondim_${model}_${depth}km_0.02.val FALSE
#######################

#./bin/xcreate_slice betav_kernel_smooth $PROJECT/ANAT_ani/event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_betav_kernel_smooth_mini_${model}_${depth}km_0.02.val TRUE

#./bin/xcreate_slice betah_kernel_smooth $PROJECT/ANAT_ani/event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_betah_kernel_smooth_mini_${model}_${depth}km_0.02.val TRUE

#./bin/xcreate_slice Gc_kernel_smooth $PROJECT/ANAT_ani/event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_Gc_kernel_smooth_mini_${model}_${depth}km_0.02.val TRUE

#./bin/xcreate_slice Gs_kernel_smooth $PROJECT/ANAT_ani/event_kernel_${model}/OUTPUT_SUM DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ani/slice_Gs_kernel_smooth_mini_${model}_${depth}km_0.02.val TRUE

#./bin/xcreate_slice dvsvs ../../subsample_test/specfem/vsv_ckbd DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_dvsvvs_Mani.set1_${depth}km_0.02.val TRUE
#./bin/xcreate_slice dvsvs ../../subsample_test/specfem/vsh_ckbd DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_dvshvs_Mani.set1_${depth}km_0.02.val TRUE
#./bin/xcreate_slice dvsvs_smooth ../../subsample_test/specfem/vs_ckbd_intg DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_dvsvs_smooth_intg_ckbd_${depth}km_0.02.val TRUE
#./bin/xcreate_slice vs ../optimize/model_ckbd/DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_vs_ckbd_${depth}km_0.02.val FALSE
#./bin/xcreate_slice vs DATABASES_MPI DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices/slice_vs_init_${depth}km_0.02.val FALSE
#./bin/xcreate_slice betav_kernel_smooth ../../subsample_test/event_kernel_Mani.set1/OUTPUT_SUM_radial DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ckbd_radial/slice_betav_kernel_Mani.set1_${depth}km_0.02.val TRUE
#./bin/xcreate_slice betah_kernel_smooth ../../subsample_test/event_kernel_Mani.set1/OUTPUT_SUM_radial DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ckbd_radial/slice_betah_kernel_Mani.set1_${depth}km_0.02.val TRUE
#./bin/xcreate_slice Gc_kernel_smooth ../../subsample_test/event_kernel_Mani.set1/OUTPUT_SUM_radial DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ckbd_radial/slice_Gc_kernel_Mani.set1_${depth}km_0.02.val TRUE
#./bin/xcreate_slice Gs_kernel_smooth ../../subsample_test/event_kernel_Mani.set1/OUTPUT_SUM_radial DATA/slice_${depth}km_0.02.xyz slices/slice_loc_${depth}km_0.02.dat slices_ckbd_radial/slice_Gs_kernel_Mani.set1_${depth}km_0.02.val TRUE
#done
