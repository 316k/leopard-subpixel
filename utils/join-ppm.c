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

    int depth = -1, force_depth = 0;

    // Args parsing
    ARGBEGIN

    ARG_CASE('d')
        depth = ARGI;
        force_depth = 1;

    WRONG_ARG
        usage:
        printf("usage: %s [-d depth] red.pgm green.pgm blue.pgm out.ppm\n", argv0);
        exit(1);

    ARGEND

    if(argc < 4) goto usage;

    if(!force_depth) {

        int size_r, size_g, size_b;

        FILE *file = fopen(argv[0], "r");
        read_image_header(file, &w, &h, &size_r);
        fclose(file);

        depth = size_r == 65535 ? 16 : 8;

        // R/G/B channels should have the same depth
        fopen(argv[1], "r");
        read_image_header(file, &w, &h, &size_g);
        fclose(file);

        if(size_g != size_r || size_b != size_r) {
            printf("error: all three pgm images should have the same depth\n"
                   "\t(use -d ... to force a depth)\n");
            exit(1);
        }
    }

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
