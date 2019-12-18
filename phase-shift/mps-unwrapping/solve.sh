#!/bin/bash

(cd .. ; mk)

source vars.sh


for period in $frequencies
do
    echo solve $period
    ../../solve-phase sine-h-${period}-* phase-x-${period}.png
    ../../solve-phase sine-v-${period}-* phase-y-${period}.png
done


../../unwrap-phase-mps -O $offset \
    $(echo $frequencies | sed -r 's/([0-9.]+)/phase-x-\1.png/g') \
    $(echo $frequencies | sed -r 's/([0-9.]+)/phase-y-\1.png/g') \
    $frequencies \
    lut.png
