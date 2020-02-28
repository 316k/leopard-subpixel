/*
  Compute a subpixel-accurate LUT from a pixel-accurate LUT and
  preprocessed dimensions
*/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

#define DEBUG 0


int nb_patterns, from_w, from_h, to_w, to_h, nb_bits;

uint16_t n_best_dimensions = 0, total_best_dimensions;

static inline float fmax4(float a, float b, float c, float d) {
    float ab = fmax(a, b);
    float cd = fmax(c, d);
    return fmax(ab, cd);
}

static inline float fmin4(float a, float b, float c, float d) {
    float ab = fmin(a, b);
    float cd = fmin(c, d);
    return fmin(ab, cd);
}

static inline int fargmin4(float a, float b, float c, float d) {
    if(fmin(a, b) < fmin(c, d)) {
        return a < b ? 0 : 1;
    } else {
        return c < d ? 3 : 4;
    }
}

/**
 * Bilinear interpolation of four values
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
float billy(float u, float v, float a, float b, float c, float d) {
    return (1 - v) * ((1 - u) * a + u * b) + v * ((1 - u) * c + u * d);
}

/**
 * Compute the subpixel cost between a camera code and four reference
 * pixels interpolated with (0, 0) <= (dx, dy) <= (1, 1)
 *
 * See Eq (9) of Subpixel Unsynchronized Unstructured Light, Chaima El
 * Asmi, Sébastien Roy
 */
float subpixel_cost(float dx, float dy,
                    float* cam_code,
                    float* a, float* b, float* c, float* d) {

    float dot_product = 0;
    float norm_cam = 0;
    float norm_ref = 0;

    for(int k=0; k<nb_bits; k++) {
        float interpolated_code = billy(dx, dy, a[k], b[k], c[k], d[k]);
        //total_cost += SQUARE(interpolated_code - cam_code[k]);

        dot_product += interpolated_code * cam_code[k];
        norm_cam += SQUARE(cam_code[k]);
        norm_ref += SQUARE(interpolated_code);
    }

    norm_cam = sqrt(norm_cam);
    norm_ref = sqrt(norm_ref);

    // return 1 - dot_product / (norm_ref * norm_cam);
    return acosf(dot_product / (norm_ref * norm_cam));
}

/* Start of a chunk of size (total_best_dimensions * 2) containing
 *  [dim #0 a] [dim #0 b] [dim #1 a] [dim #1 b] ... [dim #9 a] [dim #9 b]
 */
static inline size_t best_dims_start(int x, int y) {
    return y * (to_w - 1) * total_best_dimensions * 2 + x * total_best_dimensions * 2;
}

struct two_solutions {
    float dx1;
    float dy1;
    float dx2;
    float dy2;
};
typedef struct two_solutions two_solutions;

int find_intersections(
    float bi, float ci, float di, float si,
    float bj, float cj, float dj, float sj, two_solutions* sols) {

    float A1 = bj * (ci - di) + bi * (-cj + dj);
    float A2 = cj * (-bi + di) + ci * (bj - dj);

    if(A1 == 0 || A2 == 0)
        return 0;

    float B1 = -bj * ci + bi * cj;
    float B2 = si * (bj + cj - dj) + sj * (-bi - ci + di);
    float K = -cj * si + ci * sj;

    float delta = SQUARE(B1 + B2) - 4 * A1 * K;

    if(delta < 0)
        return 0;

    sols->dx1 = (-(B1 + B2) - sqrt(delta)) / (2 * A1);
    sols->dy1 = (-(B1 - B2) - sqrt(delta)) / (2 * A2);

    sols->dx2 = (-(B1 + B2) + sqrt(delta)) / (2 * A1);
    sols->dy2 = (-(B1 - B2) + sqrt(delta)) / (2 * A2);

    return 1;
}

// Gradient descent parameters
float step_decrease_rate = 0.8;
int max_iters = 100;

// Gradient descent constants
const float epsilon = 0.001;
float precision = 1e-4;

int biggest_iter_reached = 0, max_iters_reached = 0, nb_gradient_descents = 0;

/**
 * Find subpixel through a gradient descent
 */
int gradient_descent_solution(float* m, float* a, float* b, float* c, float* d,
                              float* dx, float* dy) {
    nb_gradient_descents++;

    float previous_step_size = INFINITY;

    int iter = 0;
    float step = 1;

    float x = 0.5, y = 0.5;

    for(iter=0; previous_step_size > precision && iter < max_iters; iter++) {

        // For the pixel to stay in the borders
        x = fmax(fmin(x, 1), 0);
        y = fmax(fmin(y, 1), 0);

        /* if(x < 0 || y < 0 || x > 1 || y > 1) { */
        /*     x = rand()/(float)RAND_MAX; */
        /*     y = rand()/(float)RAND_MAX; */
        /*     break; */
        /* } */

        float current_cost = subpixel_cost(x, y, m, a, b, c, d);
        float grad_x = (subpixel_cost(x + epsilon, y, m, a, b, c, d) - current_cost) / epsilon;
        float grad_y = (subpixel_cost(x, y + epsilon, m, a, b, c, d) - current_cost) / epsilon;

        float prev_x = x, prev_y = y;

        // Decrease speed over time to prevent loops or going
        // too far away from starting point
        step *= step_decrease_rate;

        x -= fmax(fmin(step * grad_x, 0.25), -0.25);
        y -= fmax(fmin(step * grad_y, 0.25), -0.25);

        previous_step_size = fmax(fabs(prev_x - x), fabs(SQUARE(prev_y - y)));
    }

    biggest_iter_reached = fmax(biggest_iter_reached, iter);

    if(iter == max_iters) {
        max_iters_reached++;
    }

    *dx = x;
    *dy = y;

    return 1;
}

int main(int argc, char** argv) {

    int nthreads = 4, quadratic = 0;
    int disable_gradient_descent = 0;
    int verbose = 1;

    char* ref_format = "leo_%03d.png";
    char* cam_format = "%03d.png";

    char* output_fname = "lutCamSubpixel.png";
    char* preproc_fname = "preproc";
    char* lut_fname = "lutCamPixel.png";

    char* dump_prefix = "dump-subpix-";

    char filename[FNAME_MAX_LEN]; // Generic filename buffer

    // Args parsing
    ARGBEGIN
        ARG_CASE('t')
        nthreads = ARGI;

    LSARG_CASE('p', "preprocessed-data")
        preproc_fname = ARGS;

    ARG_CASE('f')
        lut_fname = ARGS;

    ARG_CASE('R')
        ref_format = ARGS;

    ARG_CASE('C')
        cam_format = ARGS;

    ARG_CASE('O')
        output_fname = ARGS;

    ARG_CASE('D')
        dump_prefix = ARGS;

    ARG_CASE('L')
        lut_fname = ARGS;

    ARG_CASE('q')
        quadratic = 1;

    ARG_CASE('n')
        n_best_dimensions = ARGI;

    LARG_CASE("skip-discontinuities")

        disable_gradient_descent = 1;

    WRONG_ARG
        printf("usage: %s [-t nb_threads=%d]\n"
               "\t[-O output=%s] [-L lut=%s]\n"
               "\t[-R ref_format=\"%s\"] [-C cam_format=\"%s\"]\n"
               "\t[-D dump-prefix=\"%s\"]\n"
               "\t[-p|--preprocessed-data=\"%s\"] [-q quadratic]\n"
               "\t[-n number of dimensions to use=%d]\n"
               "\t[--skip-discontinuities]\n",
               argv0, nthreads, output_fname, lut_fname,
               ref_format, cam_format, dump_prefix, preproc_fname, (int) n_best_dimensions);
        exit(1);

    ARGEND

    omp_set_num_threads(nthreads);

    srand(time(NULL));

    // Check file size to avoid problems if sines.txt is empty
    FILE* info = fopen("sines.txt", "r");

    if(info != NULL)
        fseek(info, 0, SEEK_END);

    if(info == NULL || !ftell(info)) {
        printf("error: empty sines.txt\n");
        exit(-1);
    }
    fseek(info, 0, SEEK_SET);

    int nb_waves, nb_shifts;

    fscanf(info, "%d %d %d %d %d",
           &to_w, &to_h, &nb_waves, &nb_patterns, &nb_shifts);

    fclose(info);

    // Lecture d'une image pour trouver le from_w, from_h
    sprintf(filename, cam_format, 0);
    free_f32matrix(load_gray(filename, &from_w, &from_h));

    nb_bits = nb_patterns - 1;

    float*** cam_codes;
    float*** ref_codes;

    if(quadratic) {
        nb_bits = nb_patterns * (nb_patterns - 1) / 2;

        printf("Loading cam codes\n");
        cam_codes = load_continuous_codesq(cam_format, from_w, from_h, nb_patterns);

        printf("Loading ref codes\n");
        ref_codes = load_continuous_codesq(ref_format, to_w, to_h, nb_patterns);
    } else {
        cam_codes = malloc_f32cube(from_h, nb_bits, from_w);
        ref_codes = malloc_f32cube(to_h, nb_bits, to_w);

        // Camera *float* codes
        printf("Loading cam codes\n");
        load_continuous_codes(cam_codes, cam_format, from_w, from_h, nb_patterns, quadratic);

        // Ref *float* codes
        printf("Loading ref codes\n");
        load_continuous_codes(ref_codes, ref_format, to_w, to_h, nb_patterns, quadratic);
    }

    // Load pixel-precise LUT
    float*** matches = load_color(lut_fname, &to_w, &to_h);

    #pragma omp parallel for
    for(int i=0; i<from_h; i++)
        for(int j=0; j<from_w; j++) {

            if(matches[X][i][j] == 65535.0) {
                matches[X][i][j] = matches[Y][i][j] = matches[DIST][i][j] = -1.0;
            } else {
                matches[X][i][j] = round(matches[X][i][j] / 65535.0 * (to_w - 1));
                matches[Y][i][j] = round(matches[Y][i][j] / 65535.0 * (to_h - 1));
                matches[DIST][i][j] = PI; // matches[DIST][i][j] / 65535.0 * (nb_bits);
            }
        }

    float*** subpix_matches = malloc_f32cube(3, from_w, from_h);
    float*** discontinuity_map = malloc_f32cube(3, from_w, from_h);

    // Fill with pixel matches
    #pragma omp parallel for
    for(int i=0; i<from_h; i++) {
        for(int j=0; j<from_w; j++) {

            subpix_matches[X][i][j] = matches[X][i][j];
            subpix_matches[Y][i][j] = matches[Y][i][j];
            subpix_matches[DIST][i][j] = matches[DIST][i][j];
        }
    }

    // Load preproc file
    FILE *fp = fopen(preproc_fname, "rb");

    if(fp == NULL) {
        fprintf(stderr, "%s: can't open %s\n", argv0, preproc_fname);
        exit(-1);
    }

    // Read number of best dimensions
    fread(&total_best_dimensions, sizeof(uint16_t), 1, fp);

    if(n_best_dimensions == 0)
        n_best_dimensions = total_best_dimensions;

    /* printf("%d -- %d\n", n_best_dimensions, total_best_dimensions); */

    size_t preproc_size = (to_h - 1) * (to_w - 1) * total_best_dimensions * 2;

    // Read best_dimensions
    uint16_t* best_dimensions = (uint16_t*) malloc(sizeof(uint16_t) * preproc_size);
    fread(best_dimensions, sizeof(uint16_t), preproc_size, fp);

    assert(fgetc(fp) == -1);
    fclose(fp);

    printf("Computing subpixel LUT...\n");

    int progress_bar_increment = from_h / 50;

    if(verbose && progress_bar_increment) {
        // Progress-bar
        for(int i=0; i<from_h; i += progress_bar_increment) {
            fprintf(stderr, ".");
        }
        fprintf(stderr, "\n");
    }

    long last = time(NULL);

    #pragma omp parallel for
    for(int i=0; i<from_h; i++) {

        if(verbose && progress_bar_increment && i % progress_bar_increment == 0)
            fprintf(stderr, ".");

        for(int j=0; j<from_w; j++) {

            int match_x = (int) matches[X][i][j];
            int match_y = (int) matches[Y][i][j];
            float cost = matches[DIST][i][j];

            // Pour la bordure, on skip
            if(// Actual
               cost == -1 ||
               i == 0 || i == from_h - 1 ||
               j == 0 || j == from_w - 1 ||
               // Match sur la bordure
               match_x == 0 || match_x >= to_w - 1 ||
               match_y == 0 || match_y >= to_h - 1) {
                /* subpix_matches[X][i][j] = match_x; */
                /* subpix_matches[Y][i][j] = match_y; */
                /* subpix_matches[DIST][i][j] = cost / nb_bits; */
                continue;
            }

            /* printf("%d %d => %d %d (%f)\n", j, i, match_x, match_y, cost); */

            /* float central_cost = pixel_cost(cam_codes[i][j], ref_codes[match_y][match_x], nb_bits); */

            /* int min_neix = -1, min_neiy = -1; */
            /* float min_nei_cost = 10000; */

            /* for(int neiy=-1; neiy<=1; neiy++) */
            /*     for(int neix=-1; neix<=1; neix++) { */
            /*         if(neiy == 0 && neix == 0) */
            /*             continue; */

            /*         float c = pixel_cost(cam_codes[i][j], ref_codes[match_y + neiy][match_x + neix], nb_bits); */
            /*         printf("  %d %d -> %f\n", neix, neiy, c); */

            /*         if(c < min_nei_cost) { */
            /*             min_neix = neix; */
            /*             min_neiy = neiy; */
            /*             min_nei_cost = c; */
            /*         } */

            /*         if(c < central_cost) { */
            /*             if(1) { */
            /*                 printf("(%d %d) => (%d %d)\n", j, i, match_x, match_y); */
            /*                 printf("  %d %d -> %f\n", neix, neiy, c); */
            /*                 /\* assert(c >= central_cost); *\/ */
            /*             } */
            /*         } */
            /*     } */

            /* printf("  Min: %d %d -> %f\n", min_neix, min_neiy, min_nei_cost); */

            /* // Top left (-1, -1) */
            /* float topleft = pixel_cost(cam_codes[i][j], ref_codes[match_y + -1][match_x + -1], nb_bits) + */
            /*     pixel_cost(cam_codes[i][j], ref_codes[match_y + 0][match_x + -1], nb_bits) + */
            /*     pixel_cost(cam_codes[i][j], ref_codes[match_y + -1][match_x + 0], nb_bits), */

            /*     // Top right (+1, -1) */
            /*     topright = pixel_cost(cam_codes[i][j], ref_codes[match_y + -1][match_x + +1], nb_bits) + */
            /*     pixel_cost(cam_codes[i][j], ref_codes[match_y + 0][match_x + +1], nb_bits) + */
            /*     pixel_cost(cam_codes[i][j], ref_codes[match_y + -1][match_x + 0], nb_bits), */

            /*     // Bottom left (-1, +1) */
            /*     botleft = pixel_cost(cam_codes[i][j], ref_codes[match_y + +1][match_x + -1], nb_bits) + */
            /*     pixel_cost(cam_codes[i][j], ref_codes[match_y + 0][match_x + -1], nb_bits) + */
            /*     pixel_cost(cam_codes[i][j], ref_codes[match_y + +1][match_x + 0], nb_bits), */

            /*     // Bottom right (+1, +1) */
            /*     botright = pixel_cost(cam_codes[i][j], ref_codes[match_y + +1][match_x + +1], nb_bits) + */
            /*     pixel_cost(cam_codes[i][j], ref_codes[match_y + 0][match_x + +1], nb_bits) + */
            /*     pixel_cost(cam_codes[i][j], ref_codes[match_y + +1][match_x + 0], nb_bits); */

            /* float ca = pixel_cost(cam_codes[i][j], ref_codes[match_y + -1][match_x + -1], nb_bits), */
            /*     cb = pixel_cost(cam_codes[i][j], ref_codes[match_y + -1][match_x + +1], nb_bits), */
            /*     cc = pixel_cost(cam_codes[i][j], ref_codes[match_y + +1][match_x + -1], nb_bits), */
            /*     cd = pixel_cost(cam_codes[i][j], ref_codes[match_y + +1][match_x + +1], nb_bits); */

            /* printf("corners costs: %f %f %f %f\n", */
            /*        ca, cb, cc, cd); */

            /* float minc = fmin4(ca, cb, cc, cd); */

            /* printf("quadrants costs: %f %f %f %f\n", */
            /*        topleft, topright, botleft, botright */
            /*     ); */

            /* float minval = fmin4(topleft, topright, botleft, botright); */

            /* printf("lowest c: %d %d %d %d\n", */
            /*        ca == minc, cb == minc, cc == minc, cd == minc */
            /*     ); */

            /* printf("lowest q: %d %d %d %d\n", */
            /*        topleft == minval, topright == minval, botleft == minval, botright == minval */
            /*     ); */

            /* printf("dist: %f %f %f %f\n", */
            /*        topleft - minval, topright - minval, botleft - minval, botright - minval */
            /*     ); */

            // ---------------------------------------------------------------------------------------------------------------------------

            // code, mais avec .5 dans le sous-pixel du quadrant

            /* // Top left (-1, -1) */
            /* float middle = 0.25; */
            /* float stopleft = subpixel_cost(middle, middle, cam_codes[i][j], */
            /*                               ref_codes[match_y][match_x], */
            /*                               ref_codes[match_y][match_x - 1], */
            /*                               ref_codes[match_y - 1][match_x], */
            /*                               ref_codes[match_y - 1][match_x - 1]); */
            /* // (+1, -1) */
            /* float stopright = subpixel_cost(middle, middle, cam_codes[i][j], */
            /*                                ref_codes[match_y][match_x], */
            /*                                ref_codes[match_y][match_x + 1], */
            /*                                ref_codes[match_y - 1][match_x], */
            /*                                ref_codes[match_y - 1][match_x + 1]); */
            /* // (-1, +1) */
            /* float sbotleft = subpixel_cost(middle, middle, cam_codes[i][j], */
            /*                               ref_codes[match_y][match_x], */
            /*                               ref_codes[match_y][match_x - 1], */
            /*                               ref_codes[match_y + 1][match_x], */
            /*                               ref_codes[match_y + 1][match_x - 1]); */
            /* // (+1, +1) */
            /* float sbotright = subpixel_cost(middle, middle, cam_codes[i][j], */
            /*                                ref_codes[match_y][match_x], */
            /*                                ref_codes[match_y][match_x + 1], */
            /*                                ref_codes[match_y + 1][match_x], */
            /*                                ref_codes[match_y + 1][match_x + 1]); */

            /* float minsub = fmin4(stopleft, stopright, sbotleft, sbotright); */
            /* printf("lowest s: %d %d %d %d\n", */
            /*        stopleft == minsub, stopright == minsub, sbotleft == minsub, sbotright == minsub */
            /*     ); */


            /* if(i >= 10) */
            /*     exit(-1); */

            float best_subpix_x = match_x;
            float best_subpix_y = match_y;
            int is_best_discontinuity = 0;
            int best_quadrant = 0;
            float newcost = INFINITY;

            // Match complet (pixel et sous-pixel)
            float* solutions_x = (float*) malloc(sizeof(float) * n_best_dimensions * 2 * 4);
            float* solutions_y = (float*) malloc(sizeof(float) * n_best_dimensions * 2 * 4);
            // int is_discontinuity[] = {0, 0, 0, 0};
            size_t nb_solutions = 0;

            // quadrant = 0, 1, 2, 3
            //=> (x,y) = (-1,-1 ; +1,-1 ; -1,+1 ; +1,+1)
            for(int quadrant=0; quadrant<4; quadrant++) {

                int current_x = j, current_y = i;
                int current_match_x = match_x, current_match_y = match_y;

                int delta_x = quadrant % 2 == 0 ? -1 : +1;
                int delta_y = quadrant < 2 ? -1 : +1;

                // Is this quadrant representing a discontinuity?
                int match_ax = matches[X][i][j],
                    match_bx = matches[X][i][j + delta_x],
                    match_cx = matches[X][i + delta_y][j],
                    match_dx = matches[X][i + delta_y][j + delta_x],
                    match_ay = matches[Y][i][j],
                    match_by = matches[Y][i][j + delta_x],
                    match_cy = matches[Y][i + delta_y][j],
                    match_dy = matches[Y][i + delta_y][j + delta_x];

                int max_distance_x = fmax4(match_ax, match_bx, match_cx, match_dx)
                    - fmin4(match_ax, match_bx, match_cx, match_dx);
                int max_distance_y = fmax4(match_ay, match_by, match_cy, match_dy)
                    - fmin4(match_ay, match_by, match_cy, match_dy);

                if(max_distance_x > 1 || max_distance_y > 1) {
                    if(!disable_gradient_descent) {
                        // Yup, discontinuity
                        // is_discontinuity[quadrant] = 1;

                        // Solution par descente de gradient (ordonner a,b,c,d correctement selon delta_{x,y})
                        float dx, dy;

                        float* m = cam_codes[i][j];
                        float* a = ref_codes[match_ay][match_ax];
                        float* b = ref_codes[match_by][match_bx];
                        float* c = ref_codes[match_cy][match_cx];
                        float* d = ref_codes[match_dy][match_dx];

                        gradient_descent_solution(m, a, b, c, d, &dx, &dy);

                        nb_solutions++;

                        best_subpix_x = current_match_x + dx * delta_x;
                        best_subpix_y = current_match_y + dy * delta_y;
                        is_best_discontinuity = 1;
                        best_quadrant = quadrant;
                        newcost = subpixel_cost(dx, dy, m, a, b, c, d);
                    }

                    continue;
                }


                // Else: it's smooth here

                if(quadrant < 2)
                    current_match_y--;

                if(quadrant % 2 == 0)
                    current_match_x--;

                uint16_t* dims = best_dimensions + best_dims_start(current_x, current_y);

                for(int n=0; n<n_best_dimensions; n++) {
                    int dim_a = dims[n * 2 + 0], dim_b = dims[n * 2 + 1];

                    assert(dim_a < nb_bits);
                    assert(dim_b < nb_bits);

                    assert(current_match_x >= 0);
                    assert(current_match_y >= 0);
                    assert(current_match_x < to_w - 1);
                    assert(current_match_y < to_h - 1);

                    // Extraire les valeurs ai, bi, ci, di, zi pour chaque
                    // dimension
                    float zi = cam_codes[current_y][current_x][dim_a];
                    float zj = cam_codes[current_y][current_x][dim_b];

                    float ai = ref_codes[current_match_y][current_match_x][dim_a];
                    float bi = ref_codes[current_match_y][current_match_x + 1][dim_a];
                    float ci = ref_codes[current_match_y + 1][current_match_x][dim_a];
                    float di = ref_codes[current_match_y + 1][current_match_x + 1][dim_a];

                    if(zi > fmax4(ai, bi, ci, di) || zi < fmin4(ai, bi, ci, di)) {
                        // printf("No intersection (i)\n");
                        continue;
                    }

                    float aj = ref_codes[current_match_y][current_match_x][dim_b];
                    float bj = ref_codes[current_match_y][current_match_x + 1][dim_b];
                    float cj = ref_codes[current_match_y + 1][current_match_x][dim_b];
                    float dj = ref_codes[current_match_y + 1][current_match_x + 1][dim_b];

                    if(zj > fmax4(aj, bj, cj, dj) || zj < fmin4(aj, bj, cj, dj)) {
                        // printf("No intersection (j)\n");
                        continue;
                    }

                    // Avoir un petit tableau de 4 éléments (ou une struct?) pour chaque thread
                    // Solver pour (bi - ai, ci - ai, di - ai, zi - ai), (*j) pour chaque dimension
                    two_solutions sols;

                    if(find_intersections(bi - ai, ci - ai, di - ai, zi - ai,
                                          bj - aj, cj - aj, dj - aj, zj - aj, &sols)) {

                        // printf("%d %d => %f %f %f %f\n", dim_a, dim_b, sols.dx1, sols.dy1, sols.dx2, sols.dy2);

                        /* printf("%f = %f\n", billy(sols.dx1, sols.dy1, ai, bi, ci, di), zi); */
                        /* printf("%f = %f\n-\n", billy(sols.dx1, sols.dy1, aj, bj, cj, dj), zj); */
                        /* printf("%f = %f\n", billy(sols.dx2, sols.dy2, ai, bi, ci, di), zi); */
                        /* printf("%f = %f\n-\n", billy(sols.dx2, sols.dy2, aj, bj, cj, dj), zj); */

                        if(fmin(sols.dx1, sols.dy1) >= 0 && fmax(sols.dx1, sols.dy1) <= 1) {

                            /* assert(0 <= sols.dx1 && sols.dx1 <= 1); */
                            /* assert(0 <= sols.dy1 && sols.dy1 <= 1); */

                            solutions_x[nb_solutions * 2 + 0] = current_match_x + sols.dx1;
                            solutions_y[nb_solutions * 2 + 1] = current_match_y + sols.dy1;
                            nb_solutions++;

                            float cost = subpixel_cost(sols.dx1, sols.dy1,
                                                       cam_codes[current_y][current_x],

                                                       ref_codes[current_match_y][current_match_x],
                                                       ref_codes[current_match_y][current_match_x + 1],
                                                       ref_codes[current_match_y + 1][current_match_x],
                                                       ref_codes[current_match_y + 1][current_match_x + 1]);

                            if(DEBUG)
                                printf("----> %.6f %.6f <=> %.6f %c\n", sols.dx2, sols.dy2, cost, cost < newcost ? '*' : ' ');

                            if(cost < newcost) {
                                best_subpix_x = current_match_x + sols.dx1;
                                best_subpix_y = current_match_y + sols.dy1;
                                is_best_discontinuity = 0;
                                best_quadrant = quadrant;
                                newcost = cost;
                            }
                        }

                        if(fmin(sols.dx2, sols.dy2) >= 0 && fmax(sols.dx2, sols.dy2) <= 1) {

                            /* assert(0 <= sols.dx2 && sols.dx2 <= 1); */
                            /* assert(0 <= sols.dy2 && sols.dy2 <= 1); */

                            solutions_x[nb_solutions * 2 + 0] = current_match_x + sols.dx2;
                            solutions_y[nb_solutions * 2 + 1] = current_match_y + sols.dy2;
                            nb_solutions++;

                            float cost = subpixel_cost(sols.dx2, sols.dy2,
                                                       cam_codes[current_y][current_x],

                                                       ref_codes[current_match_y][current_match_x],
                                                       ref_codes[current_match_y][current_match_x + 1],
                                                       ref_codes[current_match_y + 1][current_match_x],
                                                       ref_codes[current_match_y + 1][current_match_x + 1]);

                            if(DEBUG)
                                printf("----> %.6f %.6f <=> %.6f %c\n", sols.dx2, sols.dy2, cost, cost < newcost ? '*' : ' ');

                            if(cost < newcost) {
                                best_subpix_x = current_match_x + sols.dx2;
                                best_subpix_y = current_match_y + sols.dy2;
                                is_best_discontinuity = 0;
                                best_quadrant = quadrant;
                                newcost = cost;
                            }
                        }
                    } else {
                        if(DEBUG)
                            printf("no sol\n");
                    }

                    /* #define DUMP_SUB 1 */
                    /* // S'il y a une solution, ajouter la solution à la liste de solutions potentielles */
                    /* if(DUMP_SUB && i == match_y && j == match_x) { */
                    /*     for(int iii=0; iii<nb_solutions; iii++) { */
                    /*         // printf("subpix_x %f\n", solutions_x[iii * 2] - current_match_x); */
                    /*         printf("subpix_both %f %f\n", solutions_x[iii * 2] - match_x, solutions_y[iii * 2+1] - match_y); */
                    /*     } */
                    /* } */
                }

                // TODO: Choisir la meilleur solution
                // best_subpix_x=... best_subpix_y=... selon (best) ou (mean) ou (median)
                // parmi les solutions possibles calculées

            }

            if(DEBUG)
                printf("GRU (%f %f) <=> (%f)\n", best_subpix_x, best_subpix_y, newcost);

            // assert(newcost <= cost);

            free(solutions_x);
            free(solutions_y);

            if(newcost <= 1) {
                /* printf("%f\n", newcost); */
                // Écrit la valeur de la carte sous-pixel
                subpix_matches[X][i][j] = best_subpix_x;
                subpix_matches[Y][i][j] = best_subpix_y;
                subpix_matches[DIST][i][j] = newcost;

                discontinuity_map[X][i][j] = is_best_discontinuity * 255;
                discontinuity_map[Y][i][j] = is_best_discontinuity * (best_quadrant / 4.0 * 255);
                discontinuity_map[DIST][i][j] = 0;
            }
        }
    }

    if(verbose)
        putchar('\n');

    printf("time: %ld\n", time(NULL) - last);

    sprintf(filename, "%s-x.dat", dump_prefix);
    FILE *fdebug = fopen(filename, "w");

    for(int i=0; i<from_h; i++) {
        for(int j=0; j<from_w; j++) {
            fprintf(fdebug, "%f ", subpix_matches[X][i][j]);
        }

        putc('\n', fdebug);
    }

    fclose(fdebug);

    sprintf(filename, "%s-y.dat", dump_prefix);
    fdebug = fopen(filename, "w");

    for(int i=0; i<from_h; i++) {
        for(int j=0; j<from_w; j++) {
            fprintf(fdebug, "%f ", subpix_matches[Y][i][j]);
        }

        putc('\n', fdebug);
    }

    fclose(fdebug);

    save_color_png("discontinuities.png", discontinuity_map, from_w, from_h, 8);

    if(verbose) {
        printf("Gradient descent:\n"
               "\tbiggest iteration reached: %d\n"
               "\tmax iteration reached %d times\n"
               "\t(%d gradient descents total)\n",
               biggest_iter_reached, max_iters_reached, nb_gradient_descents);
    }

    // Write subpixel-precise LUT
    save_color_map(
        output_fname,
        subpix_matches, from_w, from_h, to_w, to_h, PI);

    return EXIT_SUCCESS;
}
