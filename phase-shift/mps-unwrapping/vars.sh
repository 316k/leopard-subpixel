#!/bin/bash

# export w=1280 h=720 s=3 highperiod=30
export w=300 h=300 s=3 highperiod=30 noise=0.8 nb_frequencies=10

export maxperiod=$(echo $w $h | tr ' ' '\n' | sort -nr | sed 1q)
export step=$(echo "($w - $highperiod)/$nb_frequencies" | bc -l)
export frequencies=$(LC_ALL=C seq $highperiod $step $maxperiod | sed -r 's/00+$/0/')

