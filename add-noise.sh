#!/bin/bash
angle=-20
sine=`convert xc: -format "%[fx:sin( $angle *pi/180)]" info:`
            #

            # -flip \

i=0
for img in leo*.pgm
do
    
    # Create some square images for the cube
    convert $img -resize 750x750^ -gravity center -extent 750x750 top.png
    convert $img -resize 750x750 left.png
    convert $img -resize 750x750 right.png

    # top image shear.
    convert top.png -resize  260x301! -alpha set -background none \
            -shear 0x30 -rotate -60 -gravity center -crop 1400x301+0+0 \
            top_shear.png

    # left image shear
    convert left.png  -resize  260x301! -alpha set -background none \
            -shear 0x30  left_shear.png

    # right image shear
    convert right.png  -resize  260x301! -alpha set -background none \
            -shear 0x-30  right_shear.png

    # combine them.
    convert left_shear.png right_shear.png +append \
            \( top_shear.png -repage +0-149 \) \
            -background none -layers merge +repage \
            -resize 120% \
            +noise poisson \
            -matte -virtual-pixel Transparent \
            +distort AffineProjection "1,$sine,0,1,0,0" +repage \
            +distort barrel "0 0 -0.25" \
            +noise Gaussian \
            $(printf '%03d' $i).pgm
    
    let i++
done
rm phase*.pgm
