#!/bin/bash

# export w=1280 h=720 s=3 highperiod=30
export w=300 h=300 s=3 noise=0.8
export highperiod=30 lowperiod=$(echo $(echo "1.43 * $w" | bc)/1 | bc)
export offset=$(echo $(echo "($lowperiod - $w) / 2" | bc)/1 | bc)

