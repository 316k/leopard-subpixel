#!/bin/bash

./clean.sh

(cd ../../ ; mk)

W=100
H=100

../../generate -w $W -h $H -n 100 -s 5
# Référence (sans sous-pixel)
../../subpixel-reference -w $W -h $H > ref.ppm

i=0
for img in leo*pgm
do
    ../../translation -x 0.1 -y 0 $img > $(printf "%03d" $i).pgm
    let i++
done

LSH_ITERATIONS=80

../../solve -i $LSH_ITERATIONS | tee solve.log

echo '======= Match error ======='
./match-error.sh matches-$[$LSH_ITERATIONS - 1].ppm

mkdir debug

../../subpixel -t 1 -d -v matches-$[$LSH_ITERATIONS - 1].ppm | tee subpixel.log
../../utils/split-ppm -r x.pgm -g y.pgm -b /dev/null debug-subpixel.ppm

echo '======= Match error (subpixel.ppm) ======='
./match-error.sh subpixel.ppm

