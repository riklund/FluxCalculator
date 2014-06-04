#!/bin/bash

if [ $# -ne 1 ]
then
	echo "Usage: $0 densityfile"
	exit 1
	return 1
fi

no=`echo $1 | sed -n -e 's@\([/a-zA-Z_0-9]*\)\([-0-9\.e]\{2,20\}\)\(\.[a-z]*\)@\2@p' | sed -e 's@[-]\{0,1\}[0-9]\.[0-9]\{4,8\}e-[0-9]\{1,3\}@0.0@g'`

gnuplot -persistent << EOF
set terminal pngcairo size 1600,1000 enhanced font 'Veranda,20'
set size square
set title "Density, g=$no"
set xlabel 'x_1 /(cm)'
set ylabel 'x_2 /(cm)'
set xrange[1:15]
set yrange[1:15]
set cbrange[0:0.2]
set pm3d map
set palette defined (0 "white", 0.01 "blue", 0.02 "green", 0.04 "red")
splot "$1" notitle

EOF
