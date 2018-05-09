/*
  Join three PGM images in a single PPM image
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../args.h"
#include "../helpers.c"

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, foo, shift, w, h;

    float x = 0.0;
    float y = 0.0;

    int depth = 16;

    // Args parsing
    ARGBEGIN

    ARG_CASE('d')
        depth = ARGI;

    WRONG_ARG
        usage:
        printf("usage: %s [-d depth] red.pgm green.pgm blue.pgm out.ppm\n", argv0);
        exit(1);

    ARGEND

    if(argc < 4) goto usage;

    float** r = load_pgm(argv[0], &w, &h);
    float** g = load_pgm(argv[1], &w, &h);
    float** b = load_pgm(argv[2], &w, &h);

    float** out[3];

    out[0] = r;
    out[1] = g;
    out[2] = b;

    save_ppm(argv[3], out, w, h, depth);

    return 0;
}
