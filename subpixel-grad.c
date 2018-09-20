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

/**
 * Bilinear interpolation
 *
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

float cost(float* match, float*** ref_codes, float x, float y, int nb_patterns) {
    int base_x = (int) x;
    int base_y = (int) y;

    float offset_x = x - base_x, offset_y = y - base_y;

    float total_cost = 0;

    for(int k=0; k<nb_patterns; k++) {

        float m = match[k],
            a = ref_codes[k][base_y][base_x],
            b = ref_codes[k][base_y][base_x + 1],
            c = ref_codes[k][base_y + 1][base_x],
            d = ref_codes[k][base_y + 1][base_x + 1];

        float match_x = cosf(m), match_y = sinf(m);

        float a_x = cosf(a), a_y = sinf(a);
        float b_x = cosf(b), b_y = sinf(b);
        float c_x = cosf(c), c_y = sinf(c);
        float d_x = cosf(d), d_y = sinf(d);

        // Interpolate vectors (instead of directly interpolating
        // angles) to avoid discontinuities in cost function
        float val_x = subpixel_value(offset_x, offset_y, a_x, b_x, c_x, d_x);
        float val_y = subpixel_value(offset_x, offset_y, a_y, b_y, c_y, d_y);

        total_cost += sqrt(fsquare(val_x - match_x) + fsquare(val_y - match_y));
    }

    return total_cost / nb_patterns;
}

int main(int argc, char** argv) {

    int nthreads = 4, i, j, k, from_w, from_h, to_w, to_h, foo,
        nb_shifts, nb_patterns, debug_surface = 0, verbose = 0, proj_lut = 0;

    char* ref_format = "leo_%d_%d_%03d_%02d.pgm";
    char* cam_format = "%03d.pgm";

    float step_decrease_rate = 0.8;
    int max_iters = 100;

    // Gradient descent constants
    const float epsilon = 0.001;
    float precision = 1e-3;

    // Args parsing
    ARGBEGIN

    ARG_CASE('p')
        proj_lut = 1;

    LSARG_CASE('d', "debug-surface")
        debug_surface = 1;

    LSARG_CASE('v', "verbose")
        verbose = 1;

    LSARG_CASE('s', "slow-gradient")
        // Switch to slower gradient descent (more precision)
        max_iters = 1000;
        step_decrease_rate = 0.99;

    ARG_CASE('t')
        nthreads = ARGI;

    WRONG_ARG
        usage:
        printf("usage: %s [-t nb_threads=%d] [-p proj_lut] [-s|--slow-gradient]\n"
               "\t[-d --debug-surface] [-v --verbose] filename\n\n"
               "Enabling --slow-gradient will take more time, but might fix\n"
               "some slightly wrong matches\n",
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

    // Generating projector LUT => read captured images as projected
    // images
    if(proj_lut) {
        ref_format = "%03d.pgm";
        cam_format = "leo_%d_%d_%03d_%02d.pgm";
        ref_phase_format = "phase_cam_%d_%d_%03d.pgm";
        cam_phase_format = "phase_ref_%d_%d_%03d.pgm";

        int tmp;
        tmp = from_w;
        from_w = to_w;
        to_w = tmp;

        tmp = from_h;
        from_h = to_h;
        to_h = tmp;
    }

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

    float*** cam_codes = load_codes(cam_phase_format, cam_format, !proj_lut, nb_patterns, nb_shifts, from_w, from_h);
    float*** ref_codes = load_codes(ref_phase_format, ref_format, proj_lut, nb_patterns, nb_shifts, to_w, to_h);

    int progress_bar_increment = from_h / 50;

    if(verbose && progress_bar_increment) {
        // Progress-bar
        for(i=0; i<from_h; i += progress_bar_increment) {
            fprintf(stderr, ".");
        }
        fprintf(stderr, "\n");
    }

    int itermax_reached = 0, max_iters_reached = 0;

    #pragma omp parallel for private(i, j, k)
    for(i=0; i<from_h; i++) {

        if(verbose && progress_bar_increment && i % progress_bar_increment == 0)
            fprintf(stderr, ".");

        for(j=0; j<from_w; j++) {

            float x = matches[X][i][j];
            float y = matches[Y][i][j];

            // Undefined matches stay undefined
            if(x < 0) {
                subpixel[X][i][j] = subpixel[Y][i][j] = subpixel[DIST][i][j] = -1;
                continue;
            }

            if(i == 0 || j == 0) {
                subpixel[X][i][j] = x;
                subpixel[Y][i][j] = y;
                subpixel[DIST][i][j] = matches[DIST][i][j];
                continue;
            }

            float* match = malloc(sizeof(float) * nb_patterns);

            for(k=0; k<nb_patterns; k++) {
                match[k] = cam_codes[k][i][j];
            }

            float previous_step_size = INFINITY;

            int iter = 0;

            float step = 1;

            for(iter=0; previous_step_size > precision && iter < max_iters; iter++) {

                if(x < 1 || y < 1 || x > to_w - 2 || y > to_h - 2)
                    break;

                float current_cost = cost(match, ref_codes, x, y, nb_patterns);
                float grad_x = (cost(match, ref_codes, x + epsilon, y, nb_patterns) - current_cost) / epsilon;
                float grad_y = (cost(match, ref_codes, x, y + epsilon, nb_patterns) - current_cost) / epsilon;

                float prev_x = x, prev_y = y;

                // Decrease speed over time to prevent loops or going
                // too far away from starting point
                step *= step_decrease_rate;

                x -= fmax(fmin(step * grad_x, 0.5), -0.5);
                y -= fmax(fmin(step * grad_y, 0.5), -0.5);

                previous_step_size = fmax(fabs(prev_x - x), fabs(fsquare(prev_y - y)));
            }

            max_iters_reached = fmax(max_iters_reached, iter);
            if(iter == max_iters) {
                itermax_reached++;
            }

            free(match);

            subpixel[X][i][j] = x + 0.5;
            subpixel[Y][i][j] = y + 0.5;

            // Keep distance information
            subpixel[DIST][i][j] = matches[DIST][i][j];
        }
    }

    if(verbose)
        fprintf(stderr, "\n");

    if(debug_surface) {
        float*** out_colormap = malloc_f32cube(3, from_w, from_h);

        for(i=0; i<from_h; i++) {
            for(j=0; j<from_w; j++) {
                int pixel_x = (int) subpixel[X][i][j];
                int pixel_y = (int) subpixel[Y][i][j];

                if(pixel_x != -1) {
                    float dx = subpixel[X][i][j] - pixel_x;
                    float dy = subpixel[Y][i][j] - pixel_y;

                    out_colormap[0][i][j] = 255 * dx;
                    out_colormap[1][i][j] = 255 * dy;

                } else {

                    out_colormap[0][i][j] = 0;
                    out_colormap[1][i][j] = 0;
                    out_colormap[2][i][j] = 255;

                }
            }
        }

        save_ppm("debug-subpixel.ppm", out_colormap, from_w, from_h, 8);
    }

    save_color_map(
        proj_lut ? "lutSubProj.ppm" : "lutSubCam.ppm",
        subpixel, from_w, from_h, to_w, to_h, nb_patterns * PI/2.0);

    if(verbose) {
        fprintf(stderr, "Itermax reached %d\n", itermax_reached);
        fprintf(stderr, "maximum iterations %d\n", max_iters_reached);
    }

    return 0;
}
