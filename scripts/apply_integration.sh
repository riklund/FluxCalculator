#!/bin/bash
if [ $# -ne 1 ]
then
	echo "Usage: $0 arc_dir"
	exit 1
fi

for i in `ls $1`
do
	g=`echo $i | sed -n -e 's:\(arc_\)\([-0-9\.]*\)\(\.dat\):\2:p'`
	val=`scripts/./integrate_flux.py $1/$i`
	if [ $? -ne 0 ]
	then
		echo "#FAIL at " $g
	else
		echo $g $val
	fi
done