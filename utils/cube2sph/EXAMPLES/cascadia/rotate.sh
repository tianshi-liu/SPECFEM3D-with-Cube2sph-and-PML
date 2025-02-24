#module load python/3.6.8
#eid='20181027_M5.1_d14.8'
#model_dir='M20_fwd'
BASE_DIR=.
mkdir -p ${BASE_DIR}/OUTPUT_FILES_sph
python rotate_seismogram.py --fn_matrix="rot_P27" --rotate="XYZ->NEZ" --from_dir="${BASE_DIR}/OUTPUT_FILES" --to_dir="${BASE_DIR}/OUTPUT_FILES_sph" --from_template='${nt}.${sta}.BX${comp}.semd' --to_template='${nt}.${sta}.BX${comp}.sem.ascii'
