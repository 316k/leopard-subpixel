/*
  This program is used with traditional phase-shifting

  Unwraps high-frequencies phases using an unambiguous (but with poor
  precision) low-frequency reference

  Input 16bits high/low freq phase-maps, outputs the LUT
*/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

int correct_with_neighbourhood = 0;

float chop(float x, float low, float high) {
    return fmin(fmax(x, low), high);
}

float avg_around(float** phasemap, int x, int y, int w, int h) {
    int nb = 0;
    float total = 0;

    for(int k=-correct_with_neighbourhood; k<correct_with_neighbourhood; k++) {
        for(int l=-correct_with_neighbourhood; l<correct_with_neighbourhood; l++) {
            if(y + k < 0 || y + k >= h || x + l < 0 || x + l >= w)
                continue;

            total += phasemap[y + k][x + l];
            nb++;
        }
    }

    return total / nb;
}

int main(int argc, char** argv) {

    int nthreads = 4, high_period = 10, offset = 0;

    // Args parsing
    ARGBEGIN

    ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('h')
        high_period = ARGI;

    ARG_CASE('O')
        offset = ARGI;

    ARG_CASE('n')
        fprintf(stderr, "***correct-with-neighbourhood not implemented yet\n");
        correct_with_neighbourhood = ARGI;

    WRONG_ARG
        usage:
        printf("usage: %s [-t nthreads=%d] [-h high_period (pixels)=%d]\n"
               "\t[-n correct-with-neighbourhood=%d] [-O offset=%d]\n"
               "\tlow-x.png low-y.png high-x.png high-y.png lut.png\n",
               argv0, nthreads, high_period, correct_with_neighbourhood, offset);
        exit(1);

    ARGEND

    omp_set_num_threads(nthreads);

    srand(time(NULL));

    if(argc != 5) {
        fprintf(stderr, "Wrong number of args\n");
        goto usage;
    }

    int w, h;

    float** lowx = load_gray(argv[0], &w, &h);
    float** lowy = load_gray(argv[1], &w, &h);
    float** highx = load_gray(argv[2], &w, &h);
    float** highy = load_gray(argv[3], &w, &h);

    float*** lut = malloc_f32cube(3, w, h);

    for(int i=0; i<h; i++) {
        for(int j=0; j<w; j++) {

            float coarse_x = lowx[i][j] / 65535.0 * w + 0.5;
            float coarse_y = lowy[i][j] / 65535.0 * h + 0.5;

            float fine_x = highx[i][j] / 65535.0;
            float fine_y = highy[i][j] / 65535.0;

            int fringenum_x = (int) (coarse_x / high_period);
            int fringenum_y = (int) (coarse_y / high_period);

            // TODO : use correct_with_neighbourhood to fix the fine
            // fringe according to what makes the most sense in the
            // neigbourhood

            lut[X][i][j] = chop((fringenum_x + fine_x) * high_period - offset, 0, w - 1);
            lut[Y][i][j] = chop((fringenum_y + fine_y) * high_period - offset, 0, h - 1);
            lut[DIST][i][j] = 0;
        }
    }

    // TODO : Allow different sizes of cam/proj pairs
    save_color_map(argv[4], lut, w, h, w, h, 1);

    return EXIT_SUCCESS;
}
