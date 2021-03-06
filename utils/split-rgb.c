/*
  Split an RGB image in three grayscale images
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

    char* rname = "split-r.png";
    char* gname = "split-g.png";
    char* bname = "split-b.png";

    // Args parsing
    ARGBEGIN
    ARG_CASE('r')
        rname = ARGS;

    ARG_CASE('g')
        gname = ARGS;

    ARG_CASE('b')
        bname = ARGS;

    WRONG_ARG
        usage:
        printf("usage: %s [-r %s] [-g %s] [-b %s] rgb.png\n", argv0, rname, gname, bname);
        exit(1);

    ARGEND

    if(argc != 1) goto usage;

    int size;

    FILE *f = fopen(argv[0], "r");

    require_file(f, argv[0]);

    read_image_header(f, &w, &h, &size);
    fclose(f);

    int depth = size == 65535 ? 16 : 8;

    float*** in = load_color(argv[0], &w, &h);

    save_gray_png(rname, in[0], w, h, depth);
    save_gray_png(gname, in[1], w, h, depth);
    save_gray_png(bname, in[2], w, h, depth);

    return 0;
}
