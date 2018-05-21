#!/bin/bash
# usage: ./this 10-10, to build the heat-map around pixel (10, 10)

base="$1"

# 3 | 0
# -----
# 2 | 1


tac debug/$base-0 > heat-0
cp debug/$base-1 heat-1

cat heat-0 heat-1 > heat-right


tac debug/$base-3 > heat-3
cp debug/$base-2 heat-2

cat heat-3 heat-2 > heat-left-flipped

awk '{ for(i=NF; i>0; i--) printf("%f ", $i); printf("\n") }' heat-left-flipped > heat-left

paste -d" " heat-left heat-right > heat-all

../../utils/heat-map.sh heat-all > heat-$base.pgm

convert heat-$base.pgm -resize 500% heat-$base.png

rm heat-$base.pgm heat-{0,1,2,3,left,right,left-flipped}
