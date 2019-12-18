#!/bin/bash

(cd ../.. ; mk)

source vars.sh

for period in $frequencies
do
    echo $period
    ../../sines -O $offset -p $period -s $s -w $w -h $h -o "sine-h-${period}-%03d.png"
    ../../sines -O $offset -p $period -s $s -w $w -h $h -v -o "sine-v-${period}-%03d.png"
done

# Mess it up with noise
for sin in sine-{v,h}-*.png
do
    convert "$sin" -evaluate Gaussian-noise $noise "$sin-out.pgm"
    convert "$sin-out.pgm" "$sin" # Be sure to keep png grayscale-mode
done
rm *pgm

