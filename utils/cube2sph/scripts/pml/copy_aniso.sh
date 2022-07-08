#cd vsv_ckbd
#ls proc*vs.bin | 
#while read fn; do
#  prefix=`echo $fn | cut -d'.' -f1`
#  cp $fn ../DATABASES_MPI/${prefix}v.bin
#done
#cd ..
#cd vsh_ckbd
#ls proc*vs.bin |              
#while read fn; do
#  prefix=`echo $fn | cut -d'.' -f1`
#  cp $fn ../DATABASES_MPI/${prefix}h.bin
#done
#cd ..

ls DATABASES_MPI/proc*vp.bin |
while read fn; do
  prefix=`echo $fn | cut -d'.' -f1`
  cp $fn ${prefix}v.bin
  cp $fn ${prefix}h.bin
done

ls DATABASES_MPI/proc*vs.bin |
while read fn; do
  prefix=`echo $fn | cut -d'.' -f1`
  cp $fn ${prefix}v.bin
  cp $fn ${prefix}h.bin
done
