/*
  Dumps subpixel matching-distance for each pixel
*/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helpers.c"
#include "args.h"

/**
 * 0,0       1,0
 *   |---|---|
 *   | a | b | |
 *   |---|---| | v
 *   | c | d | v
 *   |---|---|
 * 0,1  ---> 1,1
 *       u
 */
static inline float subpixel_value(float u, float v, float a, float b, float c, float d) {
    return (1 - v) * ((1 - u) * a + u * b) + v * ((1 - u) * c + u * d);
}

static inline float fsquare(float val) {
    return val * val;
}

/**
 * Add a multiple of 2pi to val to make it closer to ref.
 * Note : both ref and val should be in [-pi, pi]
 * Ex.:
 *   ref = 3.121938, val = -3.112159
 *   => fixed_val = 3.171026
 */
float unwrap(float ref, float val) {

    for(int i = -1; i <= 1; i++) {
        float fixed_val = val + i * 2 * PI;

        float diff = fabs(ref - fixed_val);

        if(diff < PI)
            return fixed_val;
    }

    // Shouldn't happen...
    exit(-1);
}

void dump_cost(
    float* match, float*** ref_codes, int nb_patterns, int x, int y,
    int dx, int dy, int q, int precision) {

    float** costs = malloc_f32matrix(precision + 1, precision + 1);

    for(int u=0; u<=precision; u++)
        for(int v=0; v<=precision; v++)
            costs[v][u] = 0.0;

    for(int k=0; k<nb_patterns; k++) {

        float m = match[k],
            a = unwrap(m, ref_codes[k][y][x]),
            b = unwrap(m, ref_codes[k][y][x + dx]),
            c = unwrap(m, ref_codes[k][y + dy][x]),
            d = unwrap(m, ref_codes[k][y + dy][x + dx]);

        for(int u=0; u<=precision; u++) {
            for(int v=0; v<=precision; v++) {

                /* Interpolation linéaire à l'intérieur du pixel de
                   référence et calcul de la distance au pixel
                   de caméra */
                float vall = subpixel_value(u/(float)precision, v/(float)precision, a, b, c, d);

                costs[v][u] += fsquare(m - vall);
            }
        }
    }

    char debug_filename[FNAME_MAX_LEN];
    sprintf(debug_filename, "debug/%d-%d-%d", x, y, q);
    FILE *debug = fopen(debug_filename, "w");

    for(int u=0; u<precision + 1; u++) {
        for(int v=0; v<precision + 1; v++) {
            fprintf(debug, "%f ", costs[u][v]);
        }
        fputc('\n', debug);
    }
    fclose(debug);
}

int main(int argc, char** argv) {

    int nthreads = 4, i, j, k, from_w, from_h, to_w, to_h, foo,
        nb_shifts, nb_patterns;

    char* ref_format = "leo_%d_%d_%03d_%02d.pgm";
    char* cam_format = "%03d.pgm";

    // Args parsing
    ARGBEGIN

    ARG_CASE('t')
        nthreads = ARGI;

    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] filename x y\n",
               argv0, nthreads);
    exit(1);

    ARGEND

    if(argc < 3) goto usage;

    omp_set_num_threads(nthreads);

    srand(time(NULL));

    FILE* info = fopen("sines.txt", "r");

    // Check file size to avoid problems if sines.txt is empty
    fseek(info, 0, SEEK_END);
    if(!ftell(info)) {
        printf("error: empty sines.txt\n");
        exit(-1);
    }
    fseek(info, 0, SEEK_SET);

    fscanf(info, "%d %d %d %d %d", &to_w, &to_h, &foo, &nb_patterns, &nb_shifts);
    fclose(info);

    char* ref_phase_format = "phase_ref_%d_%d_%03d.pgm";
    char* cam_phase_format = "phase_cam_%d_%d_%03d.pgm";

    float*** matches = load_ppm(argv[0], &from_w, &from_h);

    int ref_x = atoi(argv[1]), ref_y = atoi(argv[2]);

    #pragma omp parallel for private(i, j)
    for(i=0; i<from_h; i++)
        for(j=0; j<from_w; j++) {

            if(matches[X][i][j] == 65535.0) {
                matches[X][i][j] = matches[Y][i][j] = matches[DIST][i][j] = -1.0;
            } else {
                matches[X][i][j] = round(matches[X][i][j] / 65535.0 * to_w);
                matches[Y][i][j] = round(matches[Y][i][j] / 65535.0 * to_h);
                matches[DIST][i][j] = matches[DIST][i][j] / 65535.0 * (nb_patterns * PI / 2.0);
            }
        }

    float*** cam_codes = load_codes(cam_phase_format, cam_format, 1, nb_patterns, nb_shifts, from_w, from_h);
    float*** ref_codes = load_codes(ref_phase_format, ref_format, 0, nb_patterns, nb_shifts, to_w, to_h);

    // Up-right, Down-right, Down-left, Up-left
    int pos_x[] = {+1, +1, -1, -1};
    int pos_y[] = {-1, +1, +1, -1};

    i=ref_y;
    j=ref_x;

    int x = matches[X][i][j];
    int y = matches[Y][i][j];

    // Undefined matches stay undefined
    if(x < 0) {
        fprintf(stderr, "Undefined match\n");
        return EXIT_FAILURE;
    }

    float* match = malloc(sizeof(float) * nb_patterns);

    for(k=0; k<nb_patterns; k++) {
        match[k] = cam_codes[k][i][j];
    }

    // Find best subpixel value

    int q; // current quadrant

    for(q=0; q<4; q++) {

        // Assure que les bornes ne sont pas dépassées
        if((pos_x[q] == -1 && x == 0) ||
           (pos_y[q] == -1 && y == 0) ||
           (pos_x[q] == +1 && x == to_w - 1) ||
           (pos_y[q] == +1 && y == to_h - 1)) {
            fprintf(stderr, "Undefined match\n");
            return EXIT_FAILURE;
        }

        dump_cost(match, ref_codes, nb_patterns, x, y, pos_x[q], pos_y[q], q, 50);
    }
    free(match);

    return 0;
}
