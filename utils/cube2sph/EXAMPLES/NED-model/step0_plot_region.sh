#!/bin/bash
set -e 

cub2sph_dir=~/src//specfem3d-cube2sph/utils/cube2sph/
# set your study region here
lonmin=-1.5
lonmax=1.5
latmin=-1.5
latmax=1.5

# cube2sph region
xmin=-150000
xmax=150000
ymin=-150000
ymax=150000
rot=0.

# cube2sph region without pml
xmin1=-120000
xmax1=120000
ymin1=-120000
ymax1=120000
# xmin1=-249750.0
# xmax1=249750.0
# ymin1=-278906.25
# ymax1=278906.25

## do nothing .... #######################

# central 
clon=`echo $lonmin $lonmax | awk '{print 0.5*($1+$2)}'`
clat=`echo $latmin $latmax | awk '{print 0.5*($1+$2)}'`

## generate boundary gmt
module load mpi hdf5/serial netcdf/serial

info=(`echo "$xmin $xmax $ymin $ymax" |awk '{for (i=1;i<=4;i++) print $i/6371000. * 180. / 3.1415926}'`)
info1=(`echo "$xmin1 $xmax1 $ymin1 $ymax1" |awk '{for (i=1;i<=4;i++) print $i/6371000. * 180. / 3.1415926}'`)
$cub2sph_dir/bin/cube2sph_boundary_gmt $clat $clon $rot ${info[*]} 0.
echo " "
mv boundary_2d.gmt boundary_2d.gmt.pml
$cub2sph_dir/bin/cube2sph_boundary_gmt $clat $clon $rot ${info1[*]} 0.


# create study region file
: > studyregion.txt
echo $lonmin $latmin >> studyregion.txt
echo $lonmin $latmax >> studyregion.txt
echo $lonmax $latmax >> studyregion.txt
echo $lonmax $latmin >> studyregion.txt
echo $lonmin $latmin >> studyregion.txt

# module load gmt-6.0.0

# larger region
x0=`echo "$lonmin-1"|bc -l`
x1=`echo "$lonmax+1"|bc -l`
y0=`echo "$latmin-1"|bc -l`
y1=`echo "$latmax+1"|bc -l`
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
