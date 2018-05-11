#!/bin/bash

# Transforms a 2D dump of float values into a PGM heat map

w=$(awk '{ print NF }' $1 | sed 1q)
h=$(cat $1 | wc -l)

max=$(awk 'BEGIN { max=0 }
    {
        for(i=1; i<=NF; i++)
           if($i > max) max = $i
    } END { print max }' $1)

echo P2
echo $w $h
echo 255

awk '{ for(i=1; i<=NF; i++) printf("%d ", $i / '$max' * 255); printf("\n") }' $1
