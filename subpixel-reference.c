#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

float get(float** matrix, int i, int j, int w, int h) {

    if(i < 0 || i >= h || j < 0 || j >= w) {
        return 0;
    }

    return matrix[i][j];

}

int main(int argc, char** argv) {

    int nthreads = 4, i, j, w = 1920, h = 1080;

    float x = 0.0, y = 0.0;

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

    ARG_CASE('w')
        w = ARGI;

    ARG_CASE('h')
        h = ARGI;

    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] [-x del-x=%f] [-y del-y=%f]\n"
               "\t[-w %d] [-h %d]\n",
               argv0, nthreads, x, y, w, h);
    exit(1);

    ARGEND

    float*** out = malloc_f32cube(3, w, h);

    for(i=0; i < h; i++)
        for(j=0; j<w; j++) {
            out[X][i][j] = j + x;
            out[Y][i][j] = i + y;
            out[DIST][i][j] = 0;
        }

    save_color_map(NULL, out, w, h, w, h, 1);

    return 0;
}
