/*
  Compare two LUT that should encode the same thing and show the error
  between them
  (TODO : remove threshold and bad/ok stuff, use awk instead)
 */
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

    int nthreads = 4, i, j, k, w, h, center = 0;

    float decalage = 0.0, threshold = 1.3;

    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('c')
        center = 1;

    ARG_CASE('l')
        threshold = ARGF;

    WRONG_ARG
        usage:
        printf("usage: %s [-c] [-t nb_threads=%d] [-l threshold=%f] ref.ppm test.ppm\n",
               argv0, nthreads, threshold);
        exit(1);

    ARGEND

    if(argc < 2)
        goto usage;

    float*** ref = load_ppm(argv[0], &w, &h);
    float*** test = load_ppm(argv[1], &w, &h);

    double val;

    for(k=0; k < 2; k++) {
        for(i=0; i < h; i++)
            for(j=0; j<w; j++) {

                ref[k][i][j] = ref[k][i][j] / 65535.0 * (k == 0 ? w : h) + (center ? 0.5 : 0);
                test[k][i][j] = test[k][i][j] / 65535.0 * (k == 0 ? w : h);

                val = fabsl(ref[k][i][j] - test[k][i][j]);

                if(val > threshold)
                    printf("bad %c %f %f %f\n", k == 0 ? 'X' : 'Y', ref[k][i][j], test[k][i][j], val);
                else
                    printf("ok %c %f %f %f\n", k == 0 ? 'X' : 'Y', ref[k][i][j], test[k][i][j], val);
            }
    }

    return 0;
}
