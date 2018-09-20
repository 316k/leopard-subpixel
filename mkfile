# -*- shell-script -*-
NPROC=4

.all:V: generate solve subpixel translation subpixel-reference \
        error validate dump_gray_pixels sines dump_subpixel subpixel-grad

clean:V:
    rm -f generate solve subpixel translation subpixel-reference \
       error validate dump_gray_pixels sines dump_subpixel subpixel-grad

presentation.pdf: presentation.md
	pandoc -t beamer -V theme:Berkeley -V colortheme:dolphin presentation.md -o presentation.pdf

&: &.c helpers.c
	gcc -std=c99 -Wall -g $stem.c -fopenmp -lm -o $stem
