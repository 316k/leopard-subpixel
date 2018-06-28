#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>

#include "args.h"
#include "helpers.c"

int main(int argc, char** argv) {
    int i, j;
    int nthreads = 4, w = 1920, h = 1080, nb_shifts = 1;
    float freq = 0.001, vertical = 0;

    ARGBEGIN
        ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('w')
        w = ARGI;

    ARG_CASE('h')
        h = ARGI;

    ARG_CASE('f')
        freq = ARGF;

    ARG_CASE('s')
        nb_shifts = ARGI;

    ARG_CASE('v')
        vertical = 1;

    WRONG_ARG
            printf("usage: %s [-t nb_threads=%d] [-w width=%d] [-h height=%d]\n"
                   "\t\t[-f frequency=%f] [-s nb_shifts=%d] [-v vertical sine]\n",
                   argv0, nthreads, w, h, freq, nb_shifts);
            exit(0);

    ARGEND
    omp_set_num_threads(nthreads);

    // Init image matrix
    float** image = malloc_f32matrix(w, h);

    for(int shift=0; shift < nb_shifts; shift++) {
        memset(image[0], 0, sizeof(float) * h * w);

        // Cos addition
        #pragma omp parallel for private(i, j)
        for(i=0; i < h; i++) {
            for(j=0; j < w; j++) {

                float fx = 0, fy = 0;
                float phase = shift * 2.0 * PI / (float) nb_shifts;

                if(vertical)
                    fy = 2 * PI * freq;
                else
                    fx = 2 * PI * freq;

                image[i][j] = (int) (127 + 128 * sinf(fy * i + fx * j + phase) + 0.5);
            }
        }

        char filename[20];
        sprintf(filename, "sine_%c_%d_%d_%d_%02d.pgm", (vertical ? 'v' : 'h'), w, h, (int) (freq * 1000), shift);
        save_pgm(filename, image, w, h, 8);
    }

    return EXIT_SUCCESS;
}
