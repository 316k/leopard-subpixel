#!/bin/bash

./clean.sh

(cd ../../ ; mk)

W=300
H=300

../../generate -w $W -h $H -n 80 -s 5

../../subpixel-reference -w $W -h $H ref.png

i=0
for img in leo*
do
    ../../translation -x 0.15 -y -0.08 $img $(printf "%03d" $i).png
    let i++
done

LSH_ITERATIONS=201

../../solve -i $LSH_ITERATIONS | tee solve-cam.log
../../solve -p -i $LSH_ITERATIONS | tee solve-proj.log

echo '======= Match error ======='
echo Cam
./match-error.sh lutCam$[$LSH_ITERATIONS - 1].png
echo Proj
./match-error.sh lutProj$[$LSH_ITERATIONS - 1].png

# Subpixel LUTs
../../subpixel-grad -v -t 1 -d -v lutCam$[$LSH_ITERATIONS - 1].png | tee subpixel.log
../../utils/split-ppm -r x-cam.png -g y-cam.png -b /dev/null debug-subpixel.png

../../subpixel-grad -p -v -t 1 -d -v lutProj$[$LSH_ITERATIONS - 1].png | tee subpixel.log
../../utils/split-ppm -r x-proj.png -g y-proj.png -b /dev/null debug-subpixel.png

