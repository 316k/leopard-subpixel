/*
  Join three PGM images in a single PPM image
  TODO: use PNG instead
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

    int depth = 8;

    // Args parsing
    ARGBEGIN

    ARG_CASE('d')
        depth = ARGI;

        if(depth != 8 && depth != 16)
            goto usage;

    WRONG_ARG
        usage:
        printf("usage: %s [-d depth=%d] r.png g.png b.png out.png\n", argv0, depth);
        exit(1);

    ARGEND

    if(argc < 4) goto usage;

    // TODO : auto-choose depth from input images

    float** r = load_gray(argv[0], &w, &h);
    float** g = load_gray(argv[1], &w, &h);
    float** b = load_gray(argv[2], &w, &h);

    float** out[3];

    out[0] = r;
    out[1] = g;
    out[2] = b;

    save_color_png(argv[3], out, w, h, depth);

    return 0;
}
