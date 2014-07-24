#!/bin/bash
if [ $# -ne 1 ]
then
	echo "Usage: $0 input_directory"
	return 1
fi

A=`ls $1 | grep eigenTwo | sed -n -e 's@\(eigenTwo_\)\([-0-9\.]*\)\(\.dat\)@\2@p'`
for g in $A
do
	val=`cut -f2 -d' ' $1/eigenTwo_$g.dat`
	echo $g $val
done