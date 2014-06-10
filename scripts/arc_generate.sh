#!/bin/bash


dir="/net/data1/riklund/flux_output/nmax30states30contour0/"

if [ $# -ne 1 ]
then
	echo "Usage: $0 output_dir"
	exit 1
fi

rm -rf $1
mkdir -p $1


no=` ls $dir/gradients/ | sed -n -e 's@\([/a-zA-Z_]*\)\([-0-9\.e]\{2,20\}\)\(\.[a-z]*\)@\2@p' | sed -e 's@[-]\{0,1\}[0-9]\.[0-9]\{4,8\}e-[0-9]\{1,3\}@0.0@g' | sort -n`

for g in $no
do
	echo $g
	./Arc $dir"gradients/gradient_$g.dat" $dir"densities/density_$g.dat" -3 45 2000 $1/arc_$g.dat

done

