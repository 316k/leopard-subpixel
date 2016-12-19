#!/bin/bash
angle=-20
img=ref.ppm
sine=`convert xc: -format "%[fx:sin( $angle *pi/180)]" info:`
            #

            # -flip \
                # Create some square images for the cube
convert $img -resize 750x750^ -depth 16 -gravity center -extent 750x750 top.png
convert $img -resize 750x750 -depth 16 left.png
convert $img -resize 750x750 -depth 16 right.png

# top image shear.
convert top.png -resize  260x301! -alpha set -background none \
        -shear 0x30 -rotate -60 -gravity center -crop 1400x301+0+0 \
        -depth 16 top_shear.png

# left image shear
convert left.png  -resize  260x301! -alpha set -background none \
        -depth 16 -shear 0x30  left_shear.png

# right image shear
convert right.png  -resize  260x301! -alpha set -background none \
        -depth 16 -shear 0x-30  right_shear.png

# combine them.
convert left_shear.png right_shear.png +append \
        \( top_shear.png -repage +0-149 \) \
        -background none -layers merge +repage \
        -resize 120% \
        -matte -virtual-pixel Transparent \
        +distort AffineProjection "1,$sine,0,1,0,0" +repage \
        +distort barrel "0 0 -0.25" \
        -depth 16 cube-ref.ppm
