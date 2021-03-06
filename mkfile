# -*- shell-script -*-
NPROC=4

programs=generate solve subpixel translation subpixel-reference \
     error validate dump_gray_pixels sines dump_subpixel subpixel-grad \
     solve-phase unwrap-phase-coarse unwrap-phase-mps solve-bin

.all:V: $programs

clean:V:
    rm -f $programs

presentation.pdf: presentation.md
    pandoc -t beamer -V theme:Berkeley -V colortheme:dolphin presentation.md -o presentation.pdf

solve-bin: solve-bin.c helpers.c
    gcc -std=c99 -Wall -g solve-bin.c -fopenmp -lgmp -lm -o solve-bin

&: &.c helpers.c
    gcc -std=c99 -Wall -g $stem.c -fopenmp -lm -o $stem
