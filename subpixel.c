/*
  Compute a subpixel-accurate LUT from a pixel-accurate LUT and phases
  references
*/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helpers.c"
#include "args.h"

// Divide each pixel in PRECISION possible sub-pixels (with both X and Y)
#define PRECISION 10

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

float unwrap(float ref, float val) {
    float diff = INFINITY;
    float new_diff;
    // FIXME ? : modulo 2PI seulement
    for(int i = -2; i <= 2; i++) {
        new_diff = fabs(ref - (val + i * PI));

        if(new_diff > diff) {
            return val + (i - 1) * PI;
        }

        diff = new_diff;
    }

    return val + 2 * PI;
}

float f32matrix_min(float** costs, float *u, float *v, int w, int h) {
    float min = INFINITY;

    for(int i=0; i<h; i++)
        for(int j=0; j<w; j++) {
            if(costs[i][j] < min) {
                *u = j;
                *v = i;
                min = costs[i][j];
            }
        }

    return min;
}

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, from_w, from_h, to_w, to_h, foo,
        nb_shifts, nb_patterns, debug_surface = 0, verbose = 0;

    char* ref_format = "leo_%d_%d_%03d_%02d.pgm";
    char* cam_format = "%03d.pgm";

    // Args parsing
    ARGBEGIN

    LSARG_CASE('d', "debug-surface")
        debug_surface = 1;

    LSARG_CASE('v', "verbose")
        verbose = 1;

    ARG_CASE('t')
        nthreads = ARGI;

    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] [-d --debug-surface] [-v --verbose] filename\n",
               argv0, nthreads);
    exit(1);

    ARGEND

    if(argc < 1) goto usage;

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

    float*** matches = load_ppm(argv[argc - 1], &from_w, &from_h);

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

    float*** subpixel = malloc_f32cube(3, from_w, from_h);

    float*** cam_codes = load_codes(cam_phase_format, cam_format, 1, nb_patterns, nb_shifts, from_w, from_h);
    float*** ref_codes = load_codes(ref_phase_format, ref_format, 0, nb_patterns, nb_shifts, to_w, to_h);

    float** colormap = malloc_f32matrix(from_w, from_h);

    // Up-right, Down-right, Down-left, Up-left
    int pos_x[] = {+1, +1, -1, -1};
    int pos_y[] = {-1, +1, +1, -1};

    int progress_bar_increment = from_h / 50;

    if(verbose && progress_bar_increment) {
        // Progress-bar
        for(i=0; i<from_h; i += progress_bar_increment) {
            fprintf(stderr, ".", i);
        }
        fprintf(stderr, "\n");
    }

    #pragma omp parallel for private(i, j, k)
    for(i=0; i<from_h; i++) {

        if(verbose && progress_bar_increment && i % progress_bar_increment == 0)
            fprintf(stderr, ".");

        for(j=0; j<from_w; j++) {

            int x = matches[X][i][j];
            int y = matches[Y][i][j];

            // Undefined matches stay undefined
            if(x < 0) {
                colormap[i][j] = -1;
                subpixel[X][i][j] = subpixel[Y][i][j] = subpixel[DIST][i][j] = -1;
                continue;
            }

            float* match = malloc(sizeof(float) * nb_patterns);

            for(k=0; k<nb_patterns; k++) {
                match[k] = cam_codes[k][i][j];
            }

            float** costs = malloc_f32matrix(PRECISION, PRECISION);

            float quadrant_best[4];
            float decalage_x[4], decalage_y[4];

            for(k=0; k<4; k++)
                quadrant_best[k] = INFINITY;

            // Minimize subpixel value

            int q; // current quadrant

            for(q=0; q<4; q++) {

                // Assure que les bornes ne sont pas dépassées
                if((pos_x[q] == -1 && x == 0) ||
                   (pos_y[q] == -1 && y == 0) ||
                   (pos_x[q] == +1 && x == to_w - 1) ||
                   (pos_y[q] == +1 && y == to_h - 1))
                    continue;

                for(int u=0; u<PRECISION; u++)
                    for(int v=0; v<PRECISION; v++)
                        costs[u][v] = 0.0;

                for(k=0; k<nb_patterns; k++) {

                    // TODO : Check this
                    float m = match[k],
                        a = unwrap(m, ref_codes[k][y][x]),
                        b = unwrap(m, ref_codes[k][y][x + pos_x[q]]),
                        c = unwrap(m, ref_codes[k][y + pos_y[q]][x]),
                        d = unwrap(m, ref_codes[k][y + pos_y[q]][x + pos_x[q]]);

                    for(int u=0; u<PRECISION; u++) {
                        for(int v=0; v<PRECISION; v++) {

                            /* Interpolation linéaire à l'intérieur du pixel de
                               référence et calcul de la distance au pixel
                               de caméra */
                            float vall = subpixel_value(u/(float)(PRECISION - 1), v/(float)(PRECISION - 1), a, b, c, d);

                            // TODO : fabs (L1) semble donner des plus
                            // beaux résultats que fsquare (L2), à vérifier
                            costs[v][u] += fabs(m - vall);
                        }
                    }
                }

                quadrant_best[q] = f32matrix_min(costs, &decalage_x[q], &decalage_y[q], PRECISION, PRECISION);

                /* decalage_x[q] = pos_x[q] * decalage_x[q] / (2.0 * PRECISION) + 0.5; */
                /* decalage_y[q] = pos_y[q] * decalage_y[q] / (2.0 * PRECISION) + 0.5; */

                decalage_x[q] = pos_x[q] * decalage_x[q]/((float)PRECISION - 1) + 0.5;
                decalage_y[q] = pos_y[q] * decalage_y[q]/((float)PRECISION - 1) + 0.5;

                if(debug_surface && rand()/(float)RAND_MAX < 0.0001) {
                    char debug_filename[50];
                    sprintf(debug_filename, "debug/%d-%d-%d", x, y, q);
                    FILE *debug = fopen(debug_filename, "w");

                    for(int u=0; u<PRECISION; u++) {
                        for(int v=0; v<PRECISION; v++) {
                            fprintf(debug, "%f ", costs[u][v]);
                        }
                        fputc('\n', debug);
                    }
                    fclose(debug);
                }
            }

            free(match);
            free_f32matrix(costs);

            // Trouve le meilleur match de sous-pixel
            float min = quadrant_best[0];
            int index = 0, equals = 0;
            for(k=1; k<4; k++) {
                if(quadrant_best[k] < min) {
                    min = quadrant_best[k];
                    index = k;
                    equals = 0;
                } else if(quadrant_best[k] == 0) { // TODO : fabs(quadrant_best[k] - min) < 1e-6) { // good epsilon ?
                    equals++;
                }

                /* if(quadrant_best[k] != 0) { */
                /*     printf("%f\n", quadrant_best[k]); */
                /* } */
            }

            if(index != -1) {
                colormap[i][j] = index;

                if(equals) {
                    colormap[i][j] = 4;
                }

                subpixel[X][i][j] = decalage_x[index];
                subpixel[Y][i][j] = decalage_y[index];
            }

            // Keep distance information
            subpixel[DIST][i][j] = matches[DIST][i][j];

            for(k=0; k<2; k++)
                subpixel[k][i][j] += matches[k][i][j];
        }
    }

    if(verbose)
        fprintf(stderr, "\n");

    FILE* vals;

    // Statistics
    for(k=0; k<2; k++) {
        float min = INFINITY; // TODO : ?
        float avg = 0, abs_avg = 0;
        float max = -INFINITY;

        if(k == 0)
            vals = fopen("subpixel-x-vals", "w");
        else
            vals = fopen("subpixel-y-vals", "w");

        for(i=1; i<from_h-1; i++)
            for(j=1; j<from_w-1; j++) {
                avg += subpixel[k][i][j] - matches[k][i][j];
                abs_avg += fabs(subpixel[k][i][j] - matches[k][i][j]);
                min = fmin(min, subpixel[k][i][j] - matches[k][i][j]);
                max = fmax(max, subpixel[k][i][j] - matches[k][i][j]);
                fprintf(vals, "%d %f\n", (int)colormap[i][j], subpixel[k][i][j] - matches[k][i][j]);
            }

        avg /= (from_h - 2) * (from_w - 2);
        abs_avg /= (from_h - 2) * (from_w - 2);
        printf("avg=%f abs_avg=%f min=%f max=%f\n", avg, abs_avg, min, max);
        fclose(vals);
    }

    if(debug_surface) {
        // TODO : Base la colormap sur la position du sous-pixel dans le
        // pixel plutôt que sur le diff avec l'image originale : la
        // méthode peut possiblement corriger un mauvais match, donc
        // certains pixels se retrouvent avec un déplacement > 1 pixel, ce
        // qui rend la colormap incompréhensible
        float*** out_colormap = malloc_f32cube(3, from_w, from_h);
        for(i=0; i<from_h; i++) {
            for(j=0; j<from_w; j++) {

                /* Code :
                   -1 = No match = White
                   0 = Up-right = R
                   1 = Down-right = G
                   2 = Down-left = B
                   3 = Up-left = Cyan
                   4 = Multiple choices (center) = Black
                */
                switch((int) colormap[i][j]) {
                case -1:
                    out_colormap[0][i][j] = out_colormap[1][i][j] = out_colormap[2][i][j] = 255;
                    break;

                case 0:
                    out_colormap[1][i][j] = out_colormap[2][i][j] = 0;
                    out_colormap[0][i][j] = 255;
                    break;

                case 1:
                    out_colormap[0][i][j] = out_colormap[2][i][j] = 0;
                    out_colormap[1][i][j] = 255;
                    break;

                case 2:
                    out_colormap[0][i][j] = out_colormap[1][i][j] = 0;
                    out_colormap[2][i][j] = 255;
                    break;

                case 3:
                    out_colormap[1][i][j] = out_colormap[2][i][j] = 255;
                    out_colormap[0][i][j] = 0;
                    break;

                case 4:
                    out_colormap[0][i][j] = out_colormap[1][i][j] = out_colormap[2][i][j] = 0;
                    break;
                }
            }
        }

        save_ppm("debug-subpixel.ppm", out_colormap, from_w, from_h, 8);
    }

    save_color_map("subpixel.ppm", subpixel, from_w, from_h, to_w, to_h, nb_patterns * PI/2.0);

    return 0;
}
