#!/bin/bash

(cd ../.. ; mk)

source vars.sh

../../solve-phase sine_h_$w\_$h\_$lowperiod\_* low-x.png
../../solve-phase sine_v_$w\_$h\_$lowperiod\_* low-y.png

../../solve-phase sine_h_$w\_$h\_$highperiod\_* high-x.png
../../solve-phase sine_v_$w\_$h\_$highperiod\_* high-y.png


../../unwrap-phase-coarse -O $offset -h $highperiod low-x.png low-y.png high-x.png high-y.png lut.png
