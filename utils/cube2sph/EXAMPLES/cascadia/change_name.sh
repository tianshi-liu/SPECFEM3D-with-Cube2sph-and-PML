#!/bin/bash
for fname in ./*_27
do
	main_name="${fname%_27}"
	cp ${main_name} ${main_name}_8
	echo "cp ${main_name} ${main_name}_8"
	cp ${main_name}_27 MESH-default/${main_name}
	echo "cp ${main_name}_27 MESH-default/${main_name}"
done
