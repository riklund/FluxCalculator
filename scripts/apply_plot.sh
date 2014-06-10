#!/bin/bash
if [ $# -lt 2 ] || [ $# -gt 3 ]
then
	echo "Usage: $0 density_dir output_dir do_flux"
	exit 1
	return 1
fi


rm -rf $2
mkdir $2

A=`ls $1 | grep '.dat' | grep -v ~`
for file in $A
do
	oname=`echo $file | sed -e 's@\([.]*\)\(\.dat\)@\1@g'`
	echo $oname
	if [ $# -eq 2 ]
	then
		scripts/./histogram.p $1/$file > $2/${oname}.png
	else
		scripts/./fluxHisto.p $1/$file > $2/${oname}.png
	fi
done

