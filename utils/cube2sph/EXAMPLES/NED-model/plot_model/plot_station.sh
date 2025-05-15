#!/bin/bash

set -e
module load gmt-6.0.0 

# 假定台站位于 105E，30N 处
stlo=-122.25
stla=44.42
# 使用 -JE 投影（四个参数: 中心经度/中心纬度/最大震中距/图片宽度）
J=E$stlo/$stla/110/10c

:> events.txt
for i in `seq 21 29`;
do
    saclst evlo evla f `ls ../../fwat_data/P$i/* |head -1` | awk '{print $2,$3}' >> events.txt
done 

gmt begin map jpg
gmt set FORMAT_GEO_MAP=+D
gmt coast -J$J -Rg -A10000 -Ggrey

# 绘制台站位置（三角形）
echo $stlo $stla | gmt plot -St0.4c -Gblack -Bya180
# 绘制地震位置（五角星）

cat events.txt | gmt plot -Sa0.25c -Gred

# 绘制等震中距线：30度（直径 60 度）、60 度（直径 120 度）
echo $stlo $stla 60d | gmt plot -SE- -W1p,red
echo $stlo $stla 182d | gmt plot -SE- -W1p,red
echo $stlo $stla 120d | gmt plot -SE- -W1p,red
# 添加文字：30、60、90 度处
gmt text -D0c/0.3c << EOF
$stlo 0 30\232
$stlo -30 60\232
$stlo -59.9 90\232
EOF

gmt end 