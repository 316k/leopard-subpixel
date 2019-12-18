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

int main(int argc, char** argv) {

    int nthreads = 4, i, j, k, w, h;
    int xy_at_once = 0;

    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;

    LARG_CASE("xy")
        xy_at_once = 1;

    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] [--xy xy at once] lut.png\n",
               argv0, nthreads);
    exit(1);

    ARGEND

    if(argc != 1)
        goto usage;

    float*** matches = load_color(argv[0], &w, &h);

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

    if(xy_at_once) {

        for(i=0; i<h; i++)
            for(j=0; j<w; j++) {
                printf("%d %d %s\n", j, i, fabs(j - matches[0][i][j]) <= 1 && fabs(i - matches[1][i][j]) <= 1 ? "yes" : "no");
            }

    } else {
        for(k=0; k<2; k++) {

            for(i=0; i<h; i++)
                for(j=0; j<w; j++) {
                    printf("%c %d %d %s\n", k == 0 ? 'X' : 'Y', k == 0 ? j : i,
                           (int) matches[k][i][j], (k == 0 ? j : i) == (int) matches[k][i][j] ? "yes" : "no");
                }

        }

    }

    return 0;
}
