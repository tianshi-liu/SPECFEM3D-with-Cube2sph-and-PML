#!/bin/bash
set -e 
module load gmt-6.0.0

# study region
lonmin=-126
lonmax=-118.5
latmin=41
latmax=47.75

# plot this model
bounds=-R/$lonmin/$lonmax/$latmin/$latmax
proj=-JM12c

# plot topography in study region
dx=`echo $lonmax $lonmin |awk '{print ($1-$2)/127.}'`
dy=`echo $latmax $latmin |awk '{print ($1-$2)/127.}'`
echo $dx $dy $bounds
gmt surface tapered_topo.txt -Gtopo_tapered.grd -I$dx/$dy  $bounds
gmt grdinfo -C topo_tapered.grd
#vmin=`gmt grdinfo -C topo_tapered.grd |awk '{print $6}'` 
vmax=`gmt grdinfo -C topo_tapered.grd |awk '{print $7}'`
gmt makecpt -Crelief -D -Z -G-2403/2403 > topo.cpt


# plot 
gmt begin topo_tapered jpg 

gmt basemap $bounds $proj -Bxaf -Byaf 
gmt grdimage topo_tapered.grd  -Ctopo.cpt $proj -E200
gmt colorbar -Ctopo.cpt -Bxaf+l"Topography,m" $bounds
gmt coast -A1000 -W0.5p $bounds

gmt end 

\rm topo_tapered.grd topo.cpt