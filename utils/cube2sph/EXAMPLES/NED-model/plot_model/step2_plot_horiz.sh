#!/bin/bash
set -e
module load gmt-6.0.0 

depth=("dafaf" 7.5 22.5) #in km

#for material in hess_kernel;
for material in vs vp rho; 
do
  for iter in `seq 1 1`; # model index
  do
    ii=`printf %02d $iter`
    echo $ii
    for i in `seq 1 2`; # profile index
    do
      # get model difference
      prof=B$i
      filename=grdfolder/$material.iter$ii.prof$prof.grd

      xmin=`gmt grdinfo $filename |grep x_min |awk '{print $3}'`
      xmax=`gmt grdinfo $filename |grep x_min |awk '{print $5}'`
      ymin=`gmt grdinfo $filename |grep y_min |awk '{print $3}'`
      ymax=`gmt grdinfo $filename |grep y_min |awk '{print $5}'`
      bounds=-R$xmin/$xmax/$ymin/$ymax
      echo $bounds

      # cut 
      gmt grdcut $filename -Gtmp.grd $bounds 
      filename=tmp.grd

      # makecpt
      info=`gmt grdinfo $filename -C`
      vmin=`echo $info| awk '{print $6}'`
      vmax=`echo $info| awk '{print $7}'`
      echo $vmin $vmax
      # vmin=-100
      # vmax=100
      gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cvik -I > out.cpt
      gmt grd2cpt $filename -Z -D -Cvik -I  > out.cpt

      # plot
      proj=-JM12c

      gmt begin pics/${material}_iter${ii}.prof$prof  jpg

        # gmt basemap $bounds -JX12c/1.5c -Bxaf+l"Longitude,Deg" -Byaf+l"Topo,m" \
        #     -BWNeb+t"${prof}-${prof}1" -Y12c

        gmt basemap $bounds $proj -Bxaf -Byaf -BWSen+t"Depth = ${depth[$i]}km"
        
        gmt grdimage $filename -Cout.cpt -E200
        gmt coast -A20 -W1p,black

        gmt plot station.lst $bounds $proj -St0.5c -Gred

        #awk '{print $1,$2/1000}' input/slab.bnd |gmt plot $bounds $proj -W1p,black
        #awk '{print $5,$3}' slab2.txt |gmt plot -W1.0p,black
        #gmt grdcontour $filename -C100
        gmt colorbar -G$vmin/$vmax -Cout.cpt -Bxaf+l"$material,km/s"
      gmt end

      \rm tmp.grd
    done
  done
done
