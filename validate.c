/*
  Indicates if pixels in a LUT encode themselves
  (ex.: lut[0,0] => colors: x=0 y=0)
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

    int nthreads = 4, i, j, k, w, h;

    float x = 0.0, y = 0.0;

    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;
    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] lut.ppm\n",
               argv0, nthreads);
    exit(1);

    ARGEND

    if(argc < 1) {
        goto usage;
    }

    float*** matches = load_ppm(argv[argc - 1], &w, &h);

    #pragma omp parallel for private(i, j)
    for(i=0; i<h; i++)
        for(j=0; j<w; j++) {

            if(matches[X][i][j] == 65535.0) {
                matches[X][i][j] = matches[Y][i][j] = matches[DIST][i][j] = -1.0;
            } else {
                matches[X][i][j] = round(matches[X][i][j] / 65535.0 * w);
                matches[Y][i][j] = round(matches[Y][i][j] / 65535.0 * h);
            }
        }

    for(k=0; k<2; k++) {

        for(i=0; i<h; i++)
            for(j=0; j<w; j++) {
                printf("%c %d %d %s\n", k == 0 ? 'X' : 'Y', k == 0 ? j : i,
                       (int) matches[k][i][j], (k == 0 ? j : i) == (int) matches[k][i][j] ? "yes" : "no");
            }
    }

    return 0;
}
