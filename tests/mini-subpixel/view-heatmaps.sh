#!/bin/bash

heatmap=../../utils/heat-map.sh

rm -R view-heatmaps
mkdir view-heatmaps

ls debug | awk -F'-' '{ print $1 "-" $2}' | grep -Ev -- '-0$' | grep -Ev '^0-' | sort | uniq |
    while read base
    do
        ./build-region-heatmap.sh $base
        mv heat-$base.png view-heatmaps/heat-$base.png
    done
