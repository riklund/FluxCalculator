#!/bin/bash

dir="/net/data1/riklund/flux_output/nmax30states30contour0/"


no=` ls $dir/gradients/ | sed -n -e 's@\([/a-zA-Z_]*\)\([-0-9\.e]\{2,20\}\)\(\.[a-z]*\)@\2@p' | sed -e 's@[-]\{0,1\}[0-9]\.[0-9]\{4,8\}e-[0-9]\{1,3\}@0.0@g' | sort -n`

count=0
rm -rf ARC_HIST
mkdir -p ARC_HIST

for g in $no
do
	echo $g
	./Arc $dir"gradients/gradient_$g.dat" $dir"densities/density_$g.dat" -3 45 2000 arc_tmp.dat

	nname=`printf %05d $count`
	count=$((count+1))
	gnuplot << EOF
set terminal pngcairo size 1600,1000 enhanced font 'Veranda,20'
set output 'ARC_HIST/$nname.png'
set title "Flux, g=$g"
set xlabel 'Angle /(rad)'
set ylabel 'Flux /(s^{-1})'
set format y '%g'
set yra [-5E-6:6E-5]
set nokey
plot 'arc_tmp.dat' u 1:2 w l
EOF

done

