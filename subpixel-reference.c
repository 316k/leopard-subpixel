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

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, foo, shift, w = 1920, h = 1080;

    float decalage = 0.0;

    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('d')
        decalage = ARGF;
        if(!(decalage <= 1 && decalage >= 0))
            goto usage;

    ARG_CASE('w')
        w = ARGI;
    
    ARG_CASE('h')
        h = ARGI;
    
    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] [-d decalage=%f]\n"
               "\t[-w %d] [-h %d]\n",
               argv0, nthreads, decalage, w, h);
    exit(1);
    
    ARGEND;
    
    float*** out = malloc_f32cube(3, w, h);
    
    for(i=0; i < h; i++)
        for(j=0; j<w; j++) {
            out[X][i][j] = j;
            out[Y][i][j] = i + decalage;
            out[DIST][i][j] = 0;
        }
    
    save_color_map(NULL, out, w, h, w, h, 1);
    
    return 0;
}
