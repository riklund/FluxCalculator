#!/bin/bash

echo "Don't use this script"
exit 1
return 1

if [ $# -ne 3 ]
then
	echo "Usage:" $0 "nmax states contourID"
	exit 1
fi

for i in $1 $2 $3
do
	if [[ ! $i =~ ^[0-9]+$ ]]
	then
		echo "Input must be integer."
		exit 2
	fi
done


num=`ls input/$1 | grep eigenTwo | sed -n -e 's@\(eigenTwo_\)\([-0-9e\.]*\)\(\.dat\)@\2@p'`

rm -rf output/$name
mkdir -p output/$name
mkdir -p output/$name/densities
mkdir -p output/$name/gradients

for i in $num
do
	name=output/$name/fluxConfig$i.conf
	cp fluxConfig.conf $name
	sed -i -e "s:@GV@:$i:g" $name
	sed -i -e "s:@INDIR:$name$g" $name
	sed -i -e "s:@OUTDIR:$name$g" $name
done
