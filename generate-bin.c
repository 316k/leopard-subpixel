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
        nb_waves = 32, nb_patterns = 20,
        init_seed = 1337;

    float base_freq = 0.05;

    int user_defined_band = 0;
    float freq_min, freq_max;

    ARGBEGIN
        ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('w')
        w = ARGI;

    ARG_CASE('h')
        h = ARGI;

    ARG_CASE('q')
        nb_waves = ARGI;

    ARG_CASE('n')
        nb_patterns = ARGI;

    ARG_CASE('f')
        base_freq = ARGF;

    LARG_CASE("seed")
        init_seed = ARGI;

    LARG_CASE("custom-band")
        user_defined_band = 1;
    freq_min = ARGF;
    freq_max = ARGF;

    WRONG_ARG
        printf("usage: %s [-t nb_threads=%d] [-w width=%d] [-h height=%d]\n"
               "\t\t[-q nb_waves=%d] [-n nb_patterns=%d]\n"
               "\t\t[-f base_freq=%f] [--seed init_seed=%d]\n"
               "\t\t[--custom-band min max]\n",
               argv0, nthreads, w, h, nb_waves, nb_patterns,
               base_freq, init_seed);
    exit(0);

    ARGEND
        omp_set_num_threads(nthreads);

    // Init image matrix
    float** image = malloc_f32matrix(w, h);

    float* phases = malloc(sizeof(float) * nb_waves);
    float* angles = malloc(sizeof(float) * nb_waves);
    float* freqs  = malloc(sizeof(float) * nb_waves);

    // XXX : delete this
    FILE* info = fopen("sines.txt", "w+");
    fprintf(info, "%d %d %d %d %d\n", w, h, nb_waves, nb_patterns, 0);

    for(int n=0; n < nb_patterns; n++) {
        printf("Rendering pattern %03d\n", n);

        for(i=0; i<nb_waves; i++) {
            srand(init_seed + n * 900 + i);
            phases[i] = rand()/(float)RAND_MAX * 2 * PI;
            angles[i] = rand()/(float)RAND_MAX * PI;

            if(user_defined_band) {
                freqs[i] = rand()/(float)RAND_MAX * (freq_max - freq_min) + freq_min;
            } else {
                freqs[i] = base_freq * pow(2, 2 * rand()/(float)RAND_MAX - 1);
            }
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

        float** intensities = malloc_f32matrix(w, h);

        memset(image[0], 0, sizeof(float) * h * w);

        // Cos addition
        for(int wave = 0; wave < nb_waves; wave++) {

            float fx, fy, phase;

            fx = 2 * PI * freqs[wave] * sinf(angles[wave]);
            fy = 2 * PI * freqs[wave] * cosf(angles[wave]);
            phase = phases[wave];

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
                intensities[i][j] = image[i][j] = gray_scale_erfc(image[i][j] / sqrt(nb_waves));
            }
        }

        char filename[FNAME_MAX_LEN];

        sprintf(filename, "leo_%03d.png", n);

        save_gray_png(filename, image, w, h, 8);

        free_f32matrix(intensities);
    }

    fclose(info);

    return EXIT_SUCCESS;
}
