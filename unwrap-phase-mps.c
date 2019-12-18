/*
  This program is used with traditional phase-shifting

  Unwraps high-frequencies phases using gradually higher frequencies

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

float chop(float x, float low, float high) {
    return fmin(fmax(x, low), high);
}

int main(int argc, char** argv) {

    int nthreads = 4, offset = 0;

    // Args parsing
    ARGBEGIN

    ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('O')
        offset = ARGI;

    WRONG_ARG
        usage:
        printf("usage: %s [-t nthreads=%d] [-O offset=%d]\n"
               "\thighest-x.png ... lowest-x.png highest-y.png ... lowest-y.png lowest-period ... highest-period lut.png\n",
               argv0, nthreads, offset);
        exit(1);

    ARGEND

    omp_set_num_threads(nthreads);

    srand(time(NULL));

    if(argc < 7 || (argc - 1) % 3 != 0) {
        fprintf(stderr, "Wrong number of args\n");
        goto usage;
    }

    int w, h;

    int nb_freqs = (argc - 1) / 3;

    float*** phases_x = malloc(sizeof(float**) * nb_freqs);
    float*** phases_y = malloc(sizeof(float**) * nb_freqs);
    float* frequencies = malloc(sizeof(float) * nb_freqs);

    for(int i=0; i<nb_freqs; i++) {
        phases_x[i] = load_gray(argv[i], &w, &h);
        phases_y[i] = load_gray(argv[nb_freqs + i], &w, &h);
        frequencies[i] = 1.0/atof(argv[nb_freqs * 2 + i]);
    }

    float*** lut = malloc_f32cube(3, w, h);

    // Start with the lowest frequency & iteratively update next ones

    // Normalize lowest frequency
    for(int i=0; i<h; i++) {
        for(int j=0; j<w; j++) {
            phases_x[nb_freqs - 1][i][j] = phases_x[nb_freqs - 1][i][j] / 65535.0;
            phases_y[nb_freqs - 1][i][j] = phases_y[nb_freqs - 1][i][j] / 65535.0;
        }
    }

    for(int k=nb_freqs - 2; k >= 0; k--) {
        for(int i=0; i<h; i++) {
            for(int j=0; j<w; j++) {

                float prev_phi_x = phases_x[k + 1][i][j];
                float prev_phi_y = phases_y[k + 1][i][j];

                float phi_x = phases_x[k][i][j] / 65535.0;
                float phi_y = phases_y[k][i][j] / 65535.0;

                float prev_freq = frequencies[k + 1];
                float freq = frequencies[k];

                phases_x[k][i][j] = phi_x - round(phi_x - freq/prev_freq * prev_phi_x);
                phases_y[k][i][j] = phi_y - round(phi_y - freq/prev_freq * prev_phi_y);

                lut[X][i][j] = chop(phases_x[k][i][j] / frequencies[k], 0, w);
                lut[Y][i][j] = chop(phases_y[k][i][j] / frequencies[k], 0, h);
                lut[DIST][i][j] = 0;
            }
        }
        /* char filename[FNAME_MAX_LEN]; */
        /* sprintf(filename, "phases-%d.png", k); */
        /* save_color_map(filename, lut, w, h, w, h, 1); */
    }

    for(int i=0; i<h; i++) {
        for(int j=0; j<w; j++) {

            float unwrapped_x = phases_x[0][i][j];
            float unwrapped_y = phases_y[0][i][j];

            lut[X][i][j] = chop(unwrapped_x / frequencies[0] - offset, 0, w - 1);
            lut[Y][i][j] = chop(unwrapped_y / frequencies[0] - offset, 0, h - 1);
            lut[DIST][i][j] = 0;
        }
    }

    // TODO : Allow different sizes of cam/proj pairs
    save_color_map(argv[argc - 1], lut, w, h, w, h, 1);

    return EXIT_SUCCESS;
}
