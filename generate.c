/*
  Generates leopard patterns to be projected during the scan
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
    int i, j;
    int nthreads = 4, w = 1920, h = 1080,
        nb_waves = 32, nb_shifts = 3,
        nb_patterns = 40,
        init_seed = 1337, identify_filenames = 0;

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

    LSARG_CASE('i', "identify-filenames")
        identify_filenames = 1;

    LARG_CASE("seed")
        init_seed = ARGI;

    WRONG_ARG
            printf("usage: %s [-t nb_threads=%d] [-w width=%d] [-h height=%d]\n"
                   "\t\t[-q nb_waves=%d] [-s nb_shifts=%d] [-n nb_patterns=%d]\n"
                   "\t\t[-f base_freq=%f] [--seed init_seed=%d]\n"
                   "\t\t[-i|--identify-filenames]\n", argv0,
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
        }

        float** arrays[3] = {&phases, &angles, &freqs};
        for(int arr=0; arr<3; arr++) {
            for(i=0; i<nb_waves; i++) {
                fprintf(info, "%f", (*arrays[arr])[i]);

                if(i < nb_waves - 1)
                    fputc(' ', info);
            }
            fputc('\n', info);
        }

        float*** intensities = malloc_f32cube(nb_shifts, w, h);

        for(int shift=0; shift < nb_shifts; shift++) {
            memset(image[0], 0, sizeof(float) * h * w);

            // Cos addition
            for(int wave = 0; wave < nb_waves; wave++) {

                float fx, fy, phase;

                fx = 2 * PI * freqs[wave] * sinf(angles[wave]);
                fy = 2 * PI * freqs[wave] * cosf(angles[wave]);

                phase = phases[wave] + shift * 2.0 * PI / (float) nb_shifts;

                #pragma omp parallel for private(i, j)
                for(i=0; i < h; i++) {
                    for(j=0; j < w; j++) {
                        image[i][j] += cosf(fy * i + fx * j + phase);
                    }
                }
            }

            #pragma omp parallel for private(i, j)
            for(i=0; i < h; i++) {
                for(j=0; j < w; j++) {
                    intensities[shift][i][j] = image[i][j] = gray_scale_erfc(image[i][j] / sqrt(nb_waves));
                }
            }

            char filename[50];
            // sprintf(filename, "leo_%d_%d_%03d_%02d.pgm", w, h, n, shift);

            int idx = nb_shifts * n + shift;

            if(identify_filenames) {
                // leo_width_height_nb-shifts_idx
                sprintf(filename, "leo_%d_%d_%d_%03d.png", w, h, nb_shifts, idx);
            } else {
                sprintf(filename, "leo_%03d.png", idx);
            }

            save_gray_png(filename, image, w, h, 8);
        }

        char filename[50];
        sprintf(filename, "phase_ref_%d_%d_%03d.png", w, h, n);
        save_phase(intensities, filename, nb_shifts, w, h);

        free_f32cube(intensities, nb_shifts);
    }

    fclose(info);

    return EXIT_SUCCESS;
}
