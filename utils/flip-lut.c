/*
  "Flips" an image point-of-view according to a given LUT
*/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../args.h"
#include "../helpers.c"

int main(char argc, char** argv) {

    int i, j, k;
    int lut_w, lut_h, img_w, img_h;

    int nthreads = 4;
    float threshold = 0.5;

    // Args parsing
    ARGBEGIN
    ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('p')
        threshold = ARGF;

    WRONG_ARG
        usage:
        fprintf(stderr, "usage: %s [-t nthreads=%d] [-p err_threshold=%0.3f] cam-img.png proj-lut.png out.png\n",
                argv0, nthreads, threshold);
        exit(1);
    ARGEND

    omp_set_num_threads(nthreads);
    srand(time(NULL));

    if(argc != 3)
        goto usage;

    char* cam_img_name = argv[0];
    char* proj_lut_name = argv[1];
    char* out_name = argv[2];

    // Lecture de l'image pour trouver le from_w, from_h
    float** img = load_gray(cam_img_name, &img_w, &img_h);
    float*** lut = load_color(proj_lut_name, &lut_w, &lut_h);
    float** out = malloc_f32matrix(lut_w, lut_h);

    // Pour chaque pixel de la LUT, trouver l'équivalent dans l'image
    #pragma omp parallel for private(i, j)
    for(i=0; i<lut_h; i++)
        for(j=0; j<lut_w; j++) {
            float x = lut[0][i][j] / 65535.0;
            float y = lut[1][i][j] / 65535.0;

            int img_y = (int) (y * img_h + 0.5);
            int img_x = (int) (x * img_w + 0.5);

            if(lut[2][i][j] / 65535.0 < threshold) {
                out[i][j] = img[img_y][img_x];
            }
        }

    save_gray_png(out_name, out, lut_w, lut_h, 8);

    return EXIT_SUCCESS;
}
