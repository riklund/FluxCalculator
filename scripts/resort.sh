#!/bin/bash
if [ $# -ne 1 ]
then
	echo "Usage: $0 dir_name"
	exit 1
	return 1
fi


A=`mktemp`
ls $1 > $A
B=`mktemp`
cat $A | sed -e 's@\([a-zA-Z0-9_]*\)\([-0-9\.e]*\)\(\.[a-zA-Z0-9]*\)@\2@g' | sed -e 's@[-]\{0,1\}[0-9]\.[0-9]\{4,8\}e-[0-9]\{1,3\}@0.0@g'> $B
C=`mktemp`
paste $B $A > $C
rm $A
rm $B
D=`mktemp`
sort -n $C > $D
rm $C
A=`cut -f2 $D`
cat $D
echo $A
rm $D
count=0
for file in $A
do
	ext="${file##*.}"
	nname=`printf %05d $count`
	mv $1/$file $1/$nname.$ext
	count=$((count+1))
done
