#!/bin/bash

heatmap=../../utils/heat-map.sh

mkdir view-debug

for i in $(ls debug)
do
    $heatmap debug/$i > view-debug/$i.pgm
    convert view-debug/$i.pgm -scale 500% view-debug/$i.png
    rm view-debug/$i.pgm
done

