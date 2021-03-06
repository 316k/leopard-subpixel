#!/bin/bash

# export w=1280 h=720 s=3 highperiod=30
export w=300 h=300 s=3 highperiod=30 noise=0.8 exp_base=1.5

export maxperiod=$(echo $(echo $(echo $w $h | tr ' ' '\n' | sort -nr | sed 1q) "* 1.43" | bc)/1 | bc)
export frequencies=$(./steps.py $exp_base $highperiod $maxperiod | sed -r 's/00+$/0/')

export offset=$(echo $(echo "($maxperiod - $w) / 2" | bc)/1 | bc)
