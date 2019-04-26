#!/bin/bash

(cd .. ; mk)

source vars.sh


for period_raw in $frequencies
do
    period=$(echo $period_raw | sed -r 's/\.0$//')
    echo solve $period $period_raw
    ../solve-phase sine_h_$w\_$h\_$period\_* phase-x-${period_raw}.png
    ../solve-phase sine_v_$w\_$h\_$period\_* phase-y-${period_raw}.png
done


../unwrap-phase-mps \
    $(echo $frequencies | sed -r 's/([0-9.]+)/phase-x-\1.png/g') \
    $(echo $frequencies | sed -r 's/([0-9.]+)/phase-y-\1.png/g') \
    $frequencies \
    lut.png
