#!/bin/bash

if [ $# -ne 4 ]
then
	echo "Usage:" $0 "nmax states contourID dist"
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

preDir="/net/data1/riklund/flux_output/"

dir="nmax"$1"states"$2"contour"$3
outDir=$preDir/$dir"dist"$4

GV=`ls input/$dir | grep eigenTwo | sed -n -e 's@\(eigenTwo_\)\([-0-9e\.]*\)\(\.dat\)@\2@p'`

mkdir -p $outDir
mkdir -p $outDir/densities
mkdir -p $outDir/gradients

outputName="$outDir/res.txt"

for g in $GV
do
	echo "Running with g=$g"
	tmpFile=`mktemp`
	cp fluxConfig.conf $tmpFile
	sed -i -e "s:@GV@:$g:g" $tmpFile
	sed -i -e "s:@INDIR@:$dir:g" $tmpFile
	sed -i -e "s:@OUTDIR@:$outDir:g" $tmpFile
	sed -i -e "s:@STOP@:$4:g" $tmpFile
	echo $g >> $outputName
	./Flux --configFile $tmpFile >> $outputName
	rm $tmpFile
done
