set -e
module load gmt-6.0.0 

# run index
run_indx=`seq 0 0`
#param_set="dbulk dbeta drho"
param_set="beta_kernel"
#param_set="vp vs rho"
evt_indx=`seq 21 29`

# plot 
#for material in beta_kernel;
for material in $param_set; 
do
  for iter in $run_indx; # model index
  do
    ii=`printf %02d $iter`
    for prof in A; # profile index
    do
      for jj in $evt_indx; do 
        echo "iter $ii evt $jj"
        # get model difference
        filename=grdfolder/$material.iter$ii.evt$jj.prof$prof.grd

        xmin=`gmt grdinfo $filename |grep x_min |awk '{print $3}'`
        xmax=`gmt grdinfo $filename |grep x_min |awk '{print $5}'`
        ymin=`gmt grdinfo $filename |grep y_min |awk '{print $3}'`
        ymax=`gmt grdinfo $filename |grep y_min |awk '{print $5}'`
        bounds=-R$xmin/$xmax/$ymin/$ymax
        #bounds=-R74.5/373.5/$ymin/$ymax

        # cut 
        gmt grdcut $filename -Gtmp.grd $bounds 
        filename=tmp.grd

        # makecpt
        info=`gmt grdinfo $filename -C`
        vmin=`echo $info| awk '{print $6}'`
        vmax=`echo $info| awk '{print $7}'`
        v0=$vmin 
        if [[ $vmin > $vmax ]]; then 
            v0=$vmax
        fi 
        # vmin=-$v0 
        # vmax=$v0
        echo $vmin $vmax $bounds
        # vmin=-100
        # vmax=100
        gmt makecpt -T$vmin/$vmax/50+n -Z -D -Cvik -I > out.cpt
        gmt grd2cpt $filename -Z -D -Cvik -I  > out.cpt

        # plot
        proj=-JX12c/-6c

        # plot topo
        info=`gmt gmtinfo topo.txt -C`
        xmin=`echo $info | awk '{print $1}'`
        xmax=`echo $info | awk '{print $2}'`
        hmin=`echo $info | awk '{print $3}'`
        hmax=`echo $info | awk '{print $4}'`
        hmax=`echo "$hmax*1.1" |bc`
        bounds1=-R$xmin/$xmax/$hmin/$hmax
        #bounds1=-R-124.1/-120.1/$hmin/$hmax	

        gmt begin pics/${material}_iter${ii}.evt$jj.prof$prof  jpg

          gmt basemap $bounds1 -JX12c/1.5c -Bxaf+l"Longitude,Deg" -Byaf+l"Topo,m" \
              -BWNeb+t"${prof}-${prof}1" -Y12c
          gmt plot topo.txt -W0.5p,black
          awk '{print $1,$2+100}' stations.txt |gmt plot -Si0.2c -Gred 

          gmt basemap $bounds $proj -Bxaf+l"Distance,km" -Byaf+l"Depth,km" -Y-6.2c -BWSet
          gmt grdimage $filename -Cout.cpt -E200

          #awk '{print $1,$2/1000}' input/slab.bnd |gmt plot $bounds $proj -W1p,black
          #awk '{print $5,$3}' slab2.txt |gmt plot -W1.0p,black
          #gmt grdcontour $filename -C100
          gmt colorbar -G$vmin/$vmax -Cout.cpt -Bxaf+l"$material,km/s"
        gmt end

        \rm tmp.grd
      done 
    done
  done
done

