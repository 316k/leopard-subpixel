#!/bin/bash

(cd .. ; mk)

source vars.sh

for period in $frequencies
do
    echo $period
    ../../sines -p $period -s $s -w $w -h $h
    ../../sines -p $period -s $s -w $w -h $h -v
done

# Mess it up with noise
for sin in sine_*.png
do
    convert "$sin" -evaluate Gaussian-noise $noise "$sin-out.pgm"
    convert "$sin-out.pgm" "$sin" # Be sure to keep png grayscale-mode
done
rm *pgm

