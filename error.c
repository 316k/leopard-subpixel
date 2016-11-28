#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

float squared(float val) {
    return val * val;
}

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, w, h;

    float decalage = 0.0;

    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;
    
    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] ref.ppm test.ppm\n",
               argv0, nthreads);
    exit(1);
    
    ARGEND

    if(argc < 2)
        goto usage;
    
    float*** ref = load_ppm(argv[0], &w, &h);
    float*** test = load_ppm(argv[1], &w, &h);

    unsigned long bad = 0;
    double val, error = 0, total = 0;
    
    for(k=0; k < 2; k++) {
        error = 0;
        bad = 0;
        for(i=0; i < h; i++)
            for(j=0; j<w; j++) {
                
                // Quadratic error
                val = squared(ref[k][i][j] - test[k][i][j]);
                // printf("%f %f %f\n", test[k][i][j], ref[k][i][j], val);
                
                if(val < 10000) {
                    error += val;
                } else {
                    bad++;
                }
            }

        error /= (float)(w*h - bad);
        
        total += error;
        
        printf("Error %d : %f\n    Bad matches : %d\n\n", k, error, bad);
    }
    
    printf("Total error : %f\n", total);
    
    return 0;
}
