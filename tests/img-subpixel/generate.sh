#!/bin/bash

for i in $(seq 1 25)
do

    for s in $(seq 0 6)
    do
        
        in=shifts/$s/$(printf "%06d" $i).pgm

        out=leo_256_256_$(printf "%03d" $[$i - 1])_0$s.pgm

        convert "$in" -resize 50% -depth 8 "$out"

    done
done

