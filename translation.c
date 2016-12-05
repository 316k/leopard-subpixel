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

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, foo, shift, w, h;

    float x = 0.0;
    float y = 0.0;


    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('x')
        x = ARGF;
    
        if(!(x <= 1 && x >= 0))
            goto usage;

    ARG_CASE('y')
        y = ARGF;
    
        if(!(y <= 1 && y >= 0))
            goto usage;

    ARG_CASE('w')
        w = ARGI;
    
    ARG_CASE('h')
        h = ARGI;
    
    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] [-x del-x=%f]\n"
               "\t[-y del-y=%f] filename\n"
               "\twith x,y element of -0.5..0.5\n",
               argv0, nthreads, x, y);
        exit(1);
    
    ARGEND
    
    if(argc < 1) goto usage;
    
    float** in = load_pgm(argv[argc - 1], &w, &h);

    float** out = malloc_f32matrix(w, h);

    for(i=0; i < h; i++)
        for(j=0; j < w; j++) {
            out[i][j] = get(in, i, j, w, h) * (1 - x)*(1 - y)
                + get(in, i, j + 1, w, h) * (1 - y)*x
                + get(in, i + 1, j, w, h) * (1 - x)*y
                + get(in, i + 1, j + 1, w, h) * x*y;
        }
            
    save_pgm(NULL, out, w, h, 8);
    
    return 0;
}
