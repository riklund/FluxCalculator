#!/bin/bash
cat $1 | \
grep '^[-0-9\.e]\|Total flux:' | \
grep -v 2VF | \
paste - - | \
tr '\t' ' ' | \
sed -n -e 's@\([-0-9\.e]*\)\( Total flux: p_1: (\)\([-0-9\.e]*\)\(,[-0-9e+\.]*)\)\( p_2: (\)\([-0-9\.e]*\)\(,[-0-9e+\.]*)\)@\1 \3@p' | \
sed -e 's@-1.04083e-17@0.0@g' \
| sort -n