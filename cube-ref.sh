#!/bin/bash
# Crée l'image de référence du "cube"

angle=-20
img=ref.ppm
sine=`convert xc: -format "%[fx:sin( $angle *pi/180)]" info:`
out=cube-ref.ppm

# tmp files
top=$(mktemp --suffix=.png)
left=$(mktemp --suffix=.png)
right=$(mktemp --suffix=.png)
top_shear=$(mktemp --suffix=.png)
left_shear=$(mktemp --suffix=.png)
right_shear=$(mktemp --suffix=.png)

if [ -n "$1" ]
then
    img="$1"
fi

if [ -n "$2" ]
then
    out="$2"
fi

# Create some square images for the cube
convert $img -resize 750x750^ -depth 16 -gravity center -extent 750x750 $top
convert $img -resize 750x750 -depth 16 $left
convert $img -resize 750x750 -depth 16 $right

# top image shear.
convert $top -resize  260x301! -alpha set -background none \
        -shear 0x30 -rotate -60 -gravity center -crop 1400x301+0+0 \
        -depth 16 $top_shear

# left image shear
convert $left  -resize  260x301! -alpha set -background none \
        -depth 16 -shear 0x30  $left_shear

# right image shear
convert $right  -resize  260x301! -alpha set -background none \
        -depth 16 -shear 0x-30  $right_shear

# combine them.
convert $left_shear $right_shear +append \
        \( $top_shear -repage +0-149 \) \
        -background none -layers merge +repage \
        -resize 120% \
        -matte -virtual-pixel Transparent \
        +distort AffineProjection "1,$sine,0,1,0,0" +repage \
        +distort barrel "0 0 -0.25" \
        -depth 16 "$out"

rm $top $left $right $top_shear $left_shear $right_shear

