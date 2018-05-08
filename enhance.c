/*
  Utils to post-process the generated LUT
  (TODO : remove median filter, use imagemagick instead)
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "helpers.c"

#define ALL 0
#define REMOVE_NOISE 1
#define MEDIAN 2

// http://stackoverflow.com/questions/1787996/c-library-function-to-do-sort
int floatcmp(const void *a,const void *b) {
    float *x = (float *) a;
    float *y = (float *) b;

    if (*x < *y) return -1;
    else if (*x > *y) return 1;

    return 0;
}


float median(float** matrix, int x, int y, int w, int h, int size) {

    int i_start = fmax(y - size, 0);
    int i_end = fmin(y + size, h - 1);

    int j_start = fmax(x - size, 0);
    int j_end = fmin(x + size, w - 1);

    int nb_values = (i_end - i_start + 1) * (j_end - j_start + 1);

    float* values = malloc(nb_values * sizeof(float));

    int i, j, k = 0;

    for(i=i_start; i<=i_end; i++) {
        for(j=j_start; j<=j_end; j++) {
            values[k] = matrix[i][j];
            k++;
        }
    }

    qsort(values, nb_values, sizeof(float), floatcmp);

    float median = values[nb_values/2];

    free(values);

    return median;
}

float*** median_filter(float** image[3], int w, int h, int size) {
    float*** medians = malloc_f32cube(3, w, h);

    #pragma omp parallel for
    for(int i=0; i<h; i++)
        for(int j=0; j<w; j++)
            for(int k=0; k<3; k++) {
                medians[k][i][j] = median(image[k], j, i, w, h, size);
            }

    return medians;
}

int main(char argc, char** argv) {

    int i, j, k, foo, shift;

    int nthreads = 4, mode = 0, w, h;

    char* filename;
    char* patched_file;

    // Args parsing
    for(i=1; i < argc - 1; i++) {
        if(strcmp(argv[i], "-t") == 0) {
            nthreads = atoi(argv[i + 1]); i++;
        } else if(strcmp(argv[i], "-m") == 0) {
            mode = atoi(argv[i + 1]); i++;
        } else {
        usage:
            printf("usage: %s [-t nb_threads=%d] [-m mode=%d] filename.ppm\n"
                   "\n"
                   "Enhanced file is printed on stdout\n"
                   "Available modes\n"
                   "\t1: remove noise (histogram check)\n"
                   "\t2: median filter\n"
                   "\tdefault: all filters\n",
                   argv[0], mode, nthreads);
            exit(0);
        }
    }
    if(i != argc - 1) goto usage;

    filename = argv[argc - 1];

    omp_set_num_threads(nthreads);

    float*** image = load_ppm(filename, &w, &h);
    float*** tmp;

    char loop = 0;

    // TODO : noise threshold for REMOVE_NOISE & MEDIAN
    do {
        switch(mode) {
        case REMOVE_NOISE:
            fprintf(stderr, "Remove noise\n");
            break;
        case MEDIAN:
            tmp = median_filter(image, w, h, 3);
            free_f32cube(image, 3);
            image = tmp;
            break;

        default: // loop on each mode
            if(loop == 0) {
                mode = 1;
                loop = 1;
            } else {
                // Stop looping on each mode
                loop = 0;
            }
        }

        mode++;
    } while(loop);

    save_ppm(NULL, image, w, h, 16);

    return EXIT_SUCCESS;
}
