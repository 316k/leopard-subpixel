/*
  Dump pixel values from a PGM image in the format :

     x    y    pix_value_a [pix_value_b ...]
     ...
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

float squared(float val) {
    return val * val;
}

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, w, h, phase, color;

    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;
    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] a.pgm [b.pgm ...]\n", argv0);
        exit(1);

    ARGEND

    if(argc < 1)
        goto usage;

    float*** imgs = malloc(sizeof(float**) * argc);

    for(i=0; i<argc; i++)
        imgs[i] = load_pgm(argv[i], &w, &h);

    for(i=0; i < h; i++)
        for(j=0; j<w; j++) {
            printf("% 3d % 3d ", j, i);

            for(k=0; k<argc; k++)
                printf("%03d ", (int) imgs[k][i][j]);

            putchar(8); // backspace last char
            putchar('\n');
        }

    return 0;
}
