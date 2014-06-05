#!/bin/bash

if [ $# -ne 3 ] && [ $# -ne 4 ]
then
	echo "Usage:" $0 "nmax states contourID [gval]"
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

dir="nmax"$1"states"$2"contour"$3

GV=`ls input/$dir | grep eigenTwo | sed -n -e 's@\(eigenTwo_\)\([-0-9e\.]*\)\(\.dat\)@\2@p'`

if ! [ $# -eq 4 ]
then
	rm -rf output/$dir
fi

mkdir -p output/$dir
mkdir -p output/$dir/densities
mkdir -p output/$dir/gradients

outputName="output/$dir/res.txt"

for g in $GV
do
	if [ "$g" != "$4" ]
	then
		continue
	fi

	echo "Running with g=$g"
	tmpFile=`mktemp`
	cp fluxConfig.conf $tmpFile
	sed -i -e "s:@GV@:$g:g" $tmpFile
	sed -i -e "s:@INDIR@:$dir:g" $tmpFile
	sed -i -e "s:@OUTDIR@:$dir:g" $tmpFile
	echo $g >> $outputName
	./Flux --configFile $tmpFile >> $outputName
	rm $tmpFile
done
