#!/bin/bash
# Teste

(cd ../../; mk)

translation=../../translation

for x in $(9 seq -1 0.1 -0.1; echo 0; 9 seq 0.1 0.1 1)
do
    for y in $(9 seq -1 0.1 -0.1; echo 0; 9 seq 0.1 0.1 1)
    do
        ../../translation -x $x -y $y dirac.pgm > tmp.pgm
        convert tmp.pgm -scale 10000% -compress none $x%$y.png
    done
done

rm tmp.pgm
