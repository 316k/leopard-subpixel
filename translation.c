/*
  Apply sub-pixel translation to a grayscale image
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

float get(float** matrix, int i, int j, int w, int h) {

    if(i < 0)
        i = 0;
    if(i > h - 1)
        i = h - 1;
    if(j < 0)
        j = 0;
    if(j > w - 1)
        j = w - 1;

    return matrix[i][j];
}

int sign(float x) {
    return x >= 0 ? 1 : -1;
}

int main(int argc, char** argv) {

    int nthreads = 4, i, j, w, h;

    float x = 0.0;
    float y = 0.0;


    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('x')
        x = ARGF;

        if(!(x <= 1 && x >= -1))
            goto usage;

    ARG_CASE('y')
        y = ARGF;

        if(!(y <= 1 && y >= -1))
            goto usage;

    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] [-x delta_x=%f]\n"
               "\t[-y delta_y=%f] gray-in.png gray-out.png\n",
               argv0, nthreads, x, y);
        exit(1);

    ARGEND

    if(argc != 2) goto usage;

    float** in = load_gray(argv[0], &w, &h);
    char* out_fname = argv[1];

    float** out = malloc_f32matrix(w, h);

    int decalage_x = -sign(x), decalage_y = -sign(y);

    for(i=0; i < h; i++)
        for(j=0; j < w; j++) {
            out[i][j] = get(in, i, j, w, h) * (1 - fabs(x))*(1 - fabs(y))
                + get(in, i, j + decalage_x, w, h) * (1 - fabs(y))*fabs(x)
                + get(in, i + decalage_y, j, w, h) * (1 - fabs(x))*fabs(y)
                + get(in, i + decalage_y, j + decalage_x, w, h) * fabs(x*y);
        }

    save_gray_png(out_fname, out, w, h, 8);

    return 0;
}
