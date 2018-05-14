#!/bin/bash

./clean.sh

(cd ../../ ; mk)

../../generate -w 40 -h 40 -n 30 -s 5

i=0
for img in leo*pgm
do
    ../../translation -x 0.35 -y 0 $img > $(printf "%03d" $i).pgm
    let i++
done

LSH_ITERATIONS=200

../../solve -i $LSH_ITERATIONS | tee solve.log

echo '======= Match error ======='
./match-error.sh

mkdir debug

../../subpixel -d -v matches-$[$LSH_ITERATIONS - 1].ppm | tee subpixel.log
../../utils/split-ppm -r x.pgm -g y.pgm -b /dev/null debug-subpixel.ppm

