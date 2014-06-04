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
#set terminal wxt
set size square
set title "Density, g=$no"
set xlabel 'x_1 /(cm)'
set ylabel 'x_2 /(cm)'
set xrange[-2:25]
set yrange[-2:25]
set cbrange[-10:-3]
#set format cb '%te%T'
set pm3d map
set palette model RGB
set palette defined
#defined (-7 "white", -6 "blue", -5 "green", -4 "red", -3 "yellow")
splot "$1" using 1:2:(log10(\$3+\$4)) notitle 

EOF
