# -*- shell-script -*-
NPROC=4
SUBPIXEL_X=0.1
SUBPIXEL_Y=0.1
W=1920
H=1080
ITERATIONS=66

.all:V: generate solve subpixel enhance translation subpixel-reference \
        error validate dump_pixels sines

reset:V: clean-cam
	rm -f {leo,phase_ref}*.pgm sines.txt tmp.log *.png sine_*.pgm

clean-cam:V:
    rm -f {phase_cam,matches,0,1,2,3,4,5,6,7,8,9,ref}*.{pgm,ppm} subpixel-{x,y}-*

clean-all:V: reset clean-cam
    rm -f generate solve subpixel enhance translation subpixel-reference \
       error validate dump_pixels sines
    rm -f cube-ref.ppm

sines.txt: generate
     ./generate -s 3 -n 20 -w $W -h $H

000.pgm: translation sines.txt
    # ./add-noise.sh
    count=0
    for i in leo*.pgm
    do
        file=$(printf "%03d.pgm\n" $count)
        echo $i to $file
        # convert $i -flip -flop $file
        # cp $i $file
        ./translation -x $SUBPIXEL_X -y $SUBPIXEL_Y $i > $file
        let "count = count + 1"
    done

matches-final.ppm: 000.pgm solve
    rm -f matches-*.pgm
    time ./solve -i $ITERATIONS
    cp matches-$[$ITERATIONS - 1].ppm matches-final.ppm

# Crée l'image de référence décalée en termes de sous-pixels
ref.ppm: subpixel-reference
    ./subpixel-reference -w $W -h $H -x $SUBPIXEL_X -y $SUBPIXEL_Y > ref.ppm

cube-ref.ppm: ref.ppm
    ./cube-ref.sh

presentation.pdf: presentation.md
	pandoc -t beamer -V theme:Berkeley -V colortheme:dolphin presentation.md -o presentation.pdf

&: &.c helpers.c
	gcc -Wall -g $stem.c -fopenmp -lm -o $stem
