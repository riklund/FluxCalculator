#!/bin/bash
num=$(awk 'BEGIN{for(i=-0.1;i<=0.1;i+=0.01)print i}')

for i in $num
do
	name=multiConfig/densityConfig$i.conf
	cp densityConfig.conf $name
	sed -i -e "s:@GV@:$i:g" $name
done
