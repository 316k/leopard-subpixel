#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

int main(int argc, char** argv) {
    int i, j;
    int nthreads = 4, w = 1920, h = 1080,
        nb_waves = 32, nb_shifts = 3,
        nb_patterns = 40,
        init_seed = 1337;

    float base_freq = 0.02;

    ARGBEGIN
        ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('w')
        w = ARGI;

    ARG_CASE('h')
        h = ARGI;

    ARG_CASE('q')
        nb_waves = ARGI;

    ARG_CASE('s')
        nb_shifts = ARGI;

    ARG_CASE('n')
        nb_patterns = ARGI;

    ARG_CASE('f')
        base_freq = ARGF;

    LARG_CASE("seed")
        init_seed = ARGI;

    WRONG_ARG
            printf("usage: %s [-t nb_threads=%d] [-w width=%d] [-h height=%d]\n"
                   "\t\t[-q nb_waves=%d] [-s nb_shifts=%d] [-n nb_patterns=%d]\n"
                   "\t\t[-f base_freq=%f] [--seed init_seed=%d]\n", argv0,
                   nthreads, w, h, nb_waves, nb_shifts, nb_patterns, base_freq, init_seed);
            exit(0);

    ARGEND
    omp_set_num_threads(nthreads);

    // Init image matrix
    float** image = malloc_f32matrix(w, h);

    float* phases = malloc(sizeof(float) * nb_waves);
    float* angles = malloc(sizeof(float) * nb_waves);
    float* freqs  = malloc(sizeof(float) * nb_waves);

    FILE* info = fopen("sines.txt", "w+");
    fprintf(info, "%d %d %d %d %d\n", w, h, nb_waves, nb_patterns, nb_shifts);

    for(int n=0; n < nb_patterns; n++) {
        printf("Rendering pattern %03d\n", n);

        for(i=0; i<nb_waves; i++) {
            srand(init_seed + n * 900 + i);
            phases[i] = rand()/(float)RAND_MAX * 2 * PI;
            angles[i] = rand()/(float)RAND_MAX * PI;
            freqs[i] = base_freq * pow(2, 2 * rand()/(float)RAND_MAX - 1);
            fprintf(info, "%f %f %f ", phases[i], angles[i], freqs[i]);
        }

        fprintf(info, "\n");

        float*** intensities = malloc_f32cube(nb_shifts, w, h);

        for(int shift=0; shift < nb_shifts; shift++) {
            memset(image[0], 0, sizeof(float) * h * w);

            // Cos addition
            for(int wave = 0; wave < nb_waves; wave++) {
                #pragma omp parallel for private(i, j)
                for(i=0; i < h; i++) {
                    for(j=0; j < w; j++) {

                        float fx, fy, phase;

                        fx = 2 * PI * freqs[wave] * sin(angles[wave]);
                        fy = 2 * PI * freqs[wave] * cos(angles[wave]);

                        phase = phases[wave] + shift * 2.0 * PI / (float) nb_shifts;

                        image[i][j] += cosf(fy * i + fx * j + phase);
                    }
                }
            }

            #pragma omp parallel for private(i, j)
            for(i=0; i < h; i++) {
                for(j=0; j < w; j++) {
                    intensities[shift][i][j] = image[i][j] = grey_scale_erfc(image[i][j] / sqrt(nb_waves));
                }
            }

            char filename[50];
            sprintf(filename, "leo_%d_%d_%03d_%02d.pgm", w, h, n, shift);
            save_pgm(filename, image, w, h, 8);
        }

        char filename[50];
        sprintf(filename, "phase_ref_%d_%d_%03d.pgm", w, h, n);
        save_phase(intensities, filename, nb_shifts, w, h);

        free_f32cube(intensities, nb_shifts);
    }

    fclose(info);

    return EXIT_SUCCESS;
}
