from_dir='../../ANAT/optimize/LBFGS_M24_slen0.01/DATABASES_MPI'
to_dir='DATABASES_MPI_M24'
current_dir=`pwd`
cp ${from_dir}/*_vp.bin ${to_dir}
cp ${from_dir}/*_vs.bin ${to_dir}
cp ${from_dir}/*_rho.bin ${to_dir}
cd ${to_dir}
ls proc*vs.bin |
while read fn; do
  prefix=`echo $fn | cut -d'.' -f1`
  cp $fn ${prefix}v.bin
  cp $fn ${prefix}h.bin
done



