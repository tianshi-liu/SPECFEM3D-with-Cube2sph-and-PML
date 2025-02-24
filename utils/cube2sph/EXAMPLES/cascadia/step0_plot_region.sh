#!/bin/bash
set -e 
# set your study region here
lonmin=-125
lonmax=-119.5
latmin=42
latmax=46.75

# cube2sph region
xmin=-270000
xmax=270000
ymin=-300000
ymax=300000
rot=0.

# cube2sph region without pml
xmin1=-252000
xmax1=252000
ymin1=-281250
ymax1=281250
# xmin1=-249750.0
# xmax1=249750.0
# ymin1=-278906.25
# ymax1=278906.25

## do nothing .... #######################

# central 
clon=`echo $lonmin $lonmax | awk '{print 0.5*($1+$2)}'`
clat=`echo $latmin $latmax | awk '{print 0.5*($1+$2)}'`

## generate boundary gmt
module load intel hdf5 netcdf

info=(`echo "$xmin $xmax $ymin $ymax" |awk '{for (i=1;i<=4;i++) print $i/6371000. * 180. / 3.1415926}'`)
info1=(`echo "$xmin1 $xmax1 $ymin1 $ymax1" |awk '{for (i=1;i<=4;i++) print $i/6371000. * 180. / 3.1415926}'`)
./bin/cube2sph_boundary_gmt $clat $clon $rot ${info[*]} 0.
echo " "
mv boundary_2d.gmt boundary_2d.gmt.pml
./bin/cube2sph_boundary_gmt $clat $clon $rot ${info1[*]} 0.


# create study region file
: > studyregion.txt
echo $lonmin $latmin >> studyregion.txt
echo $lonmin $latmax >> studyregion.txt
echo $lonmax $latmax >> studyregion.txt
echo $lonmax $latmin >> studyregion.txt
echo $lonmin $latmin >> studyregion.txt

module load gmt-6.0.0

# larger region
x0=`echo "$lonmin-3"|bc -l`
x1=`echo "$lonmax+3"|bc -l`
y0=`echo "$latmin-3"|bc -l`
y1=`echo "$latmax+3"|bc -l`
bounds1=-R/$x0/$x1/$y0/$y1
proj=-JM12c

gmt begin region jpg 
gmt basemap $bounds1 $proj -Bxaf -Byaf 
gmt plot boundary_2d.gmt.pml -W0.5p,blue -l"Simu Region + PML" $bounds1 
gmt plot boundary_2d.gmt -W0.5p,black -l"Simu Region" $bounds1 
gmt plot studyregion.txt -W0.5p,red -l"Study Region" $bounds1 
gmt coast -A1000 -W0.5p $bounds1 
awk '{print $1,$2}' station.lst  |gmt plot -St0.4c -Gblack $bounds1
#awk '{print $4,$3}' ../cascadia-cube2sph/src_rec/STATIONS_P28_globe |gmt plot -St0.5c -Gblack
gmt end 

\rm studyregion.txt 
