#!/bin/bash
set -e 
module load gmt-6.0.0

# download domain
# set bounds bigger than study region
#  1.5 deg used here 
bounds_big=-R-126.5/-118/40.5/48.25

# study region
lonmin=-125
lonmax=-119.5
latmin=42
latmax=46.75

#### STOP HERE !!!!!! ### ##################

#SRTM30s, 1km spacing
# gmt grdcut $bounds -JM12c @earth_relief_30s -Gtopo.grd 
gmt grdcut $bounds_big  @earth_relief_01m -Gtopo.grd 
gmt grd2xyz topo.grd > external_topo.txt

# plot this model
bounds=-R/$lonmin/$lonmax/$latmin/$latmax
proj=-JM12c

# plot topography in study region
dx=`echo $lonmax $lonmin |awk '{print ($1-$2)/127.}'`
dy=`echo $latmax $latmin |awk '{print ($1-$2)/127.}'`
echo $dx $dy $bounds
gmt grdcut topo.grd $bounds -Gtopo_tapered.grd
gmt grd2cpt topo_tapered.grd -Crelief -D -Z > topo.cpt

# plot 
gmt begin topo jpg 

gmt basemap $bounds $proj -Bxaf -Byaf 
gmt grdimage topo_tapered.grd  -Ctopo.cpt $proj -E200
gmt colorbar -Ctopo.cpt -Bxaf+l"Topography,m" $bounds
gmt coast -A1000 -W0.5p $bounds

gmt end 

\rm topo.cpt topo_tapered.grd topo.grd