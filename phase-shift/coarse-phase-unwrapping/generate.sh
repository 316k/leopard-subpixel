#!/bin/bash

(cd ../.. ; mk)

source vars.sh

# Lowest freq
../../sines -O $offset -p $lowperiod -s $s -w $w -h $h
../../sines -O $offset -p $lowperiod -s $s -w $w -h $h -v

# Highest freq
../../sines -O $offset -p $highperiod -s $s -w $w -h $h
../../sines -O $offset -p $highperiod -s $s -w $w -h $h -v

# Mess it up with noise
for sin in sine_*.png
do
    convert "$sin" -evaluate Gaussian-noise $noise "$sin-out.pgm"
    convert "$sin-out.pgm" "$sin" # Be sure to keep png grayscale-mode
done
rm *pgm


# # Fuckup -> make cubic
# for sin in sine_*.pgm
# do
#     ../cube-ref.sh "$sin" "$sin-out.pgm"
#     mv "$sin-out.pgm" "$sin"
# done
