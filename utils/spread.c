/*
  Randomly swap pixels around a given area
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../args.h"
#include "../helpers.c"

int randint(int low, int high) {
    return (int) rand()/(float)RAND_MAX * (high - low) + low;
}

int main(char argc, char** argv) {

    int i, j, k;
    int lut_w, lut_h;

    int area = 5;
    float probability = 1;

    // Args parsing
    ARGBEGIN
    ARG_CASE('a')
        area = ARGI;

    ARG_CASE('p')
        probability = ARGF;

    WRONG_ARG
        usage:
        printf("usage: %s [-a area=%d] [-p probability=%f] in.png spread.png\n",
               argv0, area, probability);
        exit(1);
    ARGEND

    srand(time(NULL));

    if(argc != 2)
        goto usage;

    char* proj_lut_name = argv[0];

    // Lecture de l'image pour trouver le from_w, from_h
    float*** lut = load_color(proj_lut_name, &lut_w, &lut_h);
    float*** out = malloc_f32cube(3, lut_w, lut_h);
    float** swapped = malloc_f32matrix(lut_w, lut_h);

    for(i=0; i<lut_h; i++)
        for(j=0; j<lut_w; j++) {
            out[0][i][j] = lut[0][i][j];
            out[1][i][j] = lut[1][i][j];
            out[2][i][j] = lut[2][i][j];
        }

    for(i=0; i<lut_h; i++)
        for(j=0; j<lut_w; j++) {

            if(rand()/(float)RAND_MAX > probability)
                continue;

            // Prevents double-swapping
            if(swapped[i][j])
                continue;

            int new_i = randint(fmax(i - area, 0), fmin(i + area, lut_h - 1));
            int new_j = randint(fmax(j - area, 0), fmin(j + area, lut_w - 1));

            // Try another swap if the chosen pixel has already moved
            for(int try=0; swapped[new_i][new_j] && try<5; try++) {
                new_i = randint(fmax(i - area, 0), fmin(i + area, lut_h - 1));
                new_j = randint(fmax(j - area, 0), fmin(j + area, lut_w - 1));
            }

            // No unswapped pixel randomly chosen
            if(swapped[new_i][new_j])
                continue;

            swapped[i][j] = 1;
            swapped[new_i][new_j] = 1;

            float tmp = out[0][new_i][new_j];
            out[0][new_i][new_j] = out[0][i][j];
            out[0][i][j] = tmp;

            tmp = out[1][new_i][new_j];
            out[1][new_i][new_j] = out[1][i][j];
            out[1][i][j] = tmp;

            tmp = out[2][new_i][new_j];
            out[2][new_i][new_j] = out[2][i][j];
            out[2][i][j] = tmp;
        }

    save_color_png(argv[1], out, lut_w, lut_h, 16);

    return EXIT_SUCCESS;
}
