/*
  Split a PPM image in three PGM images
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, foo, shift, w, h;

    float x = 0.0;
    float y = 0.0;

    char* rname = "split-r.pgm";
    char* gname = "split-g.pgm";
    char* bname = "split-b.pgm";

    int depth = 16;

    // Args parsing
    ARGBEGIN
    ARG_CASE('r')
        rname = ARGS;

    ARG_CASE('g')
        gname = ARGS;

    ARG_CASE('b')
        bname = ARGS;

    ARG_CASE('d')
        depth = ARGI;

    WRONG_ARG
        usage:
        printf("usage: %s [-r split-r.pgm] [-g split-g.pgm] [-b split-b.pgm] [-d depth] filename\n", argv0);
        exit(1);

    ARGEND

    if(argc < 1) goto usage;

    float*** in = load_ppm(argv[0], &w, &h);

    save_pgm(rname, in[0], w, h, depth);
    save_pgm(gname, in[1], w, h, depth);
    save_pgm(bname, in[2], w, h, depth);

    return 0;
}
