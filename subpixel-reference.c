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

    int nthreads = 4, i, j, k, shift, w = 1920, h = 1080;

    float x = 0.0, y = 0.0;

    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('x')
        x = 0.5 + ARGF / 2;
    
        if(!(x <= 1 && x >= 0))
            goto usage;

    ARG_CASE('y')
        y = 0.5 + ARGF / 2;
    
        if(!(y <= 1 && y >= 0))
            goto usage;

    ARG_CASE('w')
        w = ARGI;
    
    ARG_CASE('h')
        h = ARGI;
    
    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] [-x del-x=%f] [-y del-y=%f]\n",
               "\t[-w %d] [-h %d]\n"
               "\twith x,y element of -0.5..0.5\n",
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
    
    FILE* vals;
    
    for(k=0; k<2; k++) {
        float min = 100000;
        float dec = 0;
        float max = -100000;

        if(k == 0)
            vals = fopen("subpixel-x-ref", "w");
        else
            vals = fopen("subpixel-y-ref", "w");
        
        for(i=1; i<h-1; i++)
            for(j=1; j<w-1; j++) {
                dec += out[k][i][j] - floor(out[k][i][j]);
                min = fmin(min, out[k][i][j] - floor(out[k][i][j]));
                max = fmax(max, out[k][i][j] - floor(out[k][i][j]));
                fprintf(vals, "%f\n", out[k][i][j] - floor(out[k][i][j]));
            }
        
        dec /= (h - 2) * (w - 2);
        fprintf(stderr, "avg=%f min=%f max=%f\n", dec, min, max);
        fclose(vals);
    }
    
    save_color_map(NULL, out, w, h, w, h, 1);
    
    return 0;
}
