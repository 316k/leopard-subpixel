#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

int nb_patterns, img_w, img_h, nb_bits;


float billy(float u, float v, float a, float b, float c, float d) {
    /*
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
    return (1 - v) * ((1 - u) * a + u * b) + v * ((1 - u) * c + u * d);
}

/*
 * Fit a plane on the 4 points :
 * (0, 0, a),
 * (1, 0, b),
 * (0, 1, c),
 * (1, 1, d)
 *
 * Answer is given as a vector [_, _, _] . [1, x, y] => implicit equation
 */
static inline float fit_plane_1(float a, float b, float c, float d) {
    return 0.75*a + 0.25*c + 0.25*b - 0.25*d;
}

static inline float fit_plane_x(float a, float b, float c, float d) {
    return -0.5*a - 0.5*c + 0.5*b + 0.5*d;
}

static inline float fit_plane_y(float a, float b, float c, float d) {
    return -0.5*a + 0.5*c - 0.5*b + 0.5*d;
}

/*
planishTestABCD[a_, b_, c_, d_] := Module[{},
  Min[Abs[1 - {
      (a - c)/(b - d),
      (a - b)/(c - d)
      }]]
  ]
*/
static inline float fit_error(float a, float b, float c, float d) {
    float ratio_x1 = (a - b) / (c - d);
    float ratio_x2 = (c - d) / (a - b);

    float ratio_y1 = (a - c) / (b - d);
    float ratio_y2 = (b - d) / (a - c);

    return fmin(
        fmin(fabs(ratio_x1 - 1), fabs(ratio_x2 - 1)),
        fmin(fabs(ratio_y1 - 1), fabs(ratio_y2 - 1)));
}

int compare_idx_with_ref_asc(const void *ptr_a, const void *ptr_b, void *arg_ref) {
    float* ref = (float*) arg_ref;
    uint16_t* p_a = (uint16_t*) ptr_a;
    uint16_t* p_b = (uint16_t*) ptr_b;
    int idx_a = *p_a;
    int idx_b = *p_b;

    float a = ref[idx_a];
    float b = ref[idx_b];

    return a > b ? 1 : -1;
}

int compare_idx_with_ref_desc(const void *ptr_a, const void *ptr_b, void *arg_ref) {
    float* ref = (float*) arg_ref;
    uint16_t* p_a = (uint16_t*) ptr_a;
    uint16_t* p_b = (uint16_t*) ptr_b;
    int idx_a = *p_a;
    int idx_b = *p_b;

    float a = ref[idx_a];
    float b = ref[idx_b];

    return a < b ? 1 : -1;
}

int main(int argc, char** argv) {

    int verbose = 1;

    int nthreads = 4, quadratic = 0;
    uint16_t n_best_dimensions = 10;

    float error_threshold = 0.33;

    char* ref_format = "leo_%03d.png";
    char* fname = "preproc-sort-v3";

    // Args parsing
    ARGBEGIN
        ARG_CASE('t')
        nthreads = ARGI;

    LSARG_CASE('o', "out")
        fname = ARGS;

    ARG_CASE('R')
        ref_format = ARGS;

    ARG_CASE('q')
        quadratic = 1;

    ARG_CASE('n')
        n_best_dimensions = ARGI;

    ARG_CASE('v')
        verbose = 1;

    /* ARG_CASE('t') */
    /*     error_threshold = ARGF; */

    WRONG_ARG
        printf("usage: %s [-t nb_threads=%d] [-R ref_format=\"%s\"]\n"
               "\t\t[-n n_best_dimensions=%d]\n"
               "\t\t[-o|--out filename=%s] [-q quadratic] [-v verbose]\n",
               argv0, nthreads, ref_format, n_best_dimensions, fname);
        exit(1);

    ARGEND

    omp_set_num_threads(nthreads);

    srand(time(NULL));

    // Check file size to avoid problems if sines.txt is empty
    FILE* info = fopen("sines.txt", "r");
    nthreads = 1;
    if(info != NULL)
        fseek(info, 0, SEEK_END);

    if(info == NULL || !ftell(info)) {
        printf("error: empty sines.txt\n");
        exit(-1);
    }
    fseek(info, 0, SEEK_SET);

    int foo;
    fscanf(info, "%d %d %d %d %d",
           &img_w, &img_h, &foo, &nb_patterns, &foo);

    fclose(info);

    nb_bits = nb_patterns - 1;

    float*** codes;

    printf("Loading codes\n");
    if(quadratic) {
        nb_bits = nb_patterns * (nb_patterns - 1) / 2;
        codes = load_continuous_codesq(ref_format, img_w, img_h, nb_patterns);
    } else {
        // TODO : Réordonner en w, h, nb_bits
        codes = malloc_f32cube(img_h, nb_bits, img_w);

        // --- Load codes ---
        load_continuous_codes(codes, ref_format, img_w, img_h, nb_patterns, 0);
    }

    int nb_combinations = nb_bits * (nb_bits - 1) / 2;

    int max_nb_pairs = 4 * nb_bits;

    assert(n_best_dimensions <= max_nb_pairs);

    printf("%d %d\n", n_best_dimensions, nb_combinations);

    // (i j nb_best {0,1})
    size_t output_size = (img_h - 1) * (img_w - 1) * n_best_dimensions * 2;
    uint16_t* output = (uint16_t*) malloc(sizeof(uint16_t) * output_size);

    printf("Computing best dimensions...\n");

    int progress_bar_increment = img_h / 50;

    if(verbose && progress_bar_increment) {
        // Progress-bar
        for(int i=0; i<img_h; i += progress_bar_increment) {
            fprintf(stderr, ".");
        }
        fprintf(stderr, "\n");
    }

    #pragma omp parallel for
    for(int i=0; i<img_h - 1; i++) {

        if(verbose && progress_bar_increment && i % progress_bar_increment == 0)
            fprintf(stderr, ".");

        // Allocate buffers for this thread
        float* dim_error = malloc(sizeof(float) * nb_bits);
        uint16_t* dim_error_idx = malloc(sizeof(uint16_t) * nb_bits);

        float* dim_slope = malloc(sizeof(float) * nb_bits);
        uint16_t* dim_slope_idx = malloc(sizeof(uint16_t) * nb_bits);

        // Trigonometric 2D plane normal angle (between 0 and 2pi)
        float* dim_angles = malloc(sizeof(float) * nb_bits);
        uint16_t* dim_angles_idx = malloc(sizeof(uint16_t) * nb_bits);

        uint16_t* considered_dimensions = calloc(nb_bits, sizeof(uint16_t));

        // Sorted pairs
        uint16_t* all_pairs = malloc(sizeof(uint16_t) * 2 * max_nb_pairs);

        float* paires_scores = malloc(sizeof(float) * max_nb_pairs); // [low_idx_pos] = score;
        uint16_t* paires_scores_idx = malloc(sizeof(uint16_t) * max_nb_pairs);

        for(int j=0; j<img_w - 1; j++) {

            // POUR CHAQUE PIXEL

            // Pour chaque dimension du code du pixel
            // Filter out les dimensions qui donnent des plans trop mauvais
            //
            // -> Peu inclinés
            // -> Pas vraiment des plans
            for(uint16_t dim = 0; dim < nb_bits; dim++) {

                // Init index arrays
                dim_error_idx[dim] = dim;
                dim_slope_idx[dim] = dim;
                dim_angles_idx[dim] = dim;

                float a = codes[i][j][dim];
                float b = codes[i][j + 1][dim];
                float c = codes[i + 1][j][dim];
                float d = codes[i + 1][j + 1][dim];

                // Fit error (how close is the fitted plane to the
                // actual bilinear surface)
                dim_error[dim] = fit_error(a, b, c, d);

                float fitted_plane_x = fit_plane_x(a, b, c, d);
                float fitted_plane_y = fit_plane_y(a, b, c, d);

                // Plane slope
                float normal_x = -fitted_plane_x;
                float normal_y = -fitted_plane_y;
                const float normal_z = 1;

                float norm = sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z);

                normal_x /= norm;
                normal_y /= norm;
                // normal_z /= norm;

                float squared_slope = SQUARE(normal_x) + SQUARE(normal_y);

                dim_slope[dim] = squared_slope;

                float norm2D = sqrt(squared_slope);
                float min_angle = acosf(normal_x/norm2D) / (2 * PI) * 360;
                // float min_angle = arccos(normal_x/norm2D) / (2 * PI) * 360;

                float angle = min_angle;

                if(-normal_y < 0) {
                    angle = 360 - min_angle;
                }

                dim_angles[dim] = angle;
            }

            qsort_r(dim_angles_idx, nb_bits, sizeof(uint16_t), compare_idx_with_ref_asc, dim_angles);
            qsort_r(dim_error_idx, nb_bits, sizeof(uint16_t), compare_idx_with_ref_asc, dim_error);
            qsort_r(dim_slope_idx, nb_bits, sizeof(uint16_t), compare_idx_with_ref_desc, dim_slope);

            int nb_considered_dimensions = 0;

            // Keep only the dimensions below error_threshold
            for(int i=0; i<nb_bits; i++) {
                uint16_t idx = dim_error_idx[i];

                if(dim_error[idx] > error_threshold)
                    break;

                considered_dimensions[idx] = 1;
                nb_considered_dimensions++;
            }

            /* // Remove the {1 - slope_threshold}% planes with the worst slopes */
            /* for(int i=0; i<=(1 - slope_threshold) * nb_bits; i++) { */
            /*     uint16_t idx = dim_error_idx[nb_bits - i - 1]; */
            /*     if(considered_dimensions[idx]) { */
            /*         nb_considered_dimensions--; */
            /*         considered_dimensions[idx] = 0; */
            /*     } */
            /* } */


            // Number of pairs considered for the final selection of
            // the n best
            int nb_pairs = 0;

            // low index position, high index position
            int low_idx_pos = 0, high_idx_pos = 0;

            // Ensure low/high_idx_pos are in considered_dimensions[]
            while(!considered_dimensions[dim_angles_idx[low_idx_pos]])
                low_idx_pos = (low_idx_pos + 1) % nb_bits;

            while(!considered_dimensions[dim_angles_idx[high_idx_pos]])
                high_idx_pos = (high_idx_pos + 1) % nb_bits;

            // Increment high_idx_pos until [high - low] ~= 90deg

            float diff = fabs(dim_angles[dim_angles_idx[high_idx_pos + 1]] - dim_angles[dim_angles_idx[low_idx_pos]]);

            while(diff < 90 && high_idx_pos < nb_bits - 1) {

                uint16_t low_idx = dim_angles_idx[low_idx_pos];
                uint16_t high_idx = dim_angles_idx[high_idx_pos + 1];

                diff = fabs(dim_angles[low_idx] - dim_angles[high_idx]);
                high_idx_pos++;
            }

            while(low_idx_pos < nb_bits) {

                uint16_t low_idx = dim_angles_idx[low_idx_pos];
                uint16_t high_idx = dim_angles_idx[high_idx_pos];

                // Si low
                float low_angle = dim_angles[low_idx];
                float high_angle = dim_angles[high_idx];

                // TODO : le score devrait être
                float score = fmin(dim_slope[low_idx], dim_slope[high_idx]);


                // Add (low_idx, high_idx) to the list of pairs
                // Format: [low_0, high_0, low_1, high_1, ... ]
                /* printf("guiglu #%d [%d %d] -> (%f %f) => %f\n", low_idx_pos, low_idx, high_idx, low_angle, high_angle, score); */
                /* printf("nb_pairs=%d\n", nb_pairs); */
                /* assert(nb_pairs * 2 < 2 * 2 * nb_bits); */
                /* assert(nb_pairs * 2 + 1 < 2 * 2 * nb_bits); */

                if(low_idx != high_idx) {
                    all_pairs[nb_pairs * 2 + 0] = low_idx;
                    all_pairs[nb_pairs * 2 + 1] = high_idx;

                    // Indexes are sorted wrt. their scores afterwards
                    paires_scores_idx[nb_pairs] = nb_pairs;
                    paires_scores[nb_pairs] = score;

                    nb_pairs++;
                }

                // Incrémenter soit low, soit high, selon ce qui donne
                // un écart le plus proche de 90deg possible (en mod
                // 360)
                /* printf("guiglu #%d [%d %d] -> (%f %f) => %f\n", low_idx_pos, low_idx, high_idx, low_angle, high_angle, score); */

                // Find the next considered_dimension
                uint16_t low_next = low_idx_pos + 1;
                uint16_t high_next = (high_idx_pos + 1)  % nb_bits;

                while(low_next < nb_bits && !considered_dimensions[dim_angles_idx[low_next]])
                    low_next++;

                // No more bits available
                if(low_next == nb_bits)
                    break;

                while(!considered_dimensions[dim_angles_idx[high_next]])
                    high_next = (high_next + 1) % nb_bits;

                // ...
                // Utiliser low_next ou high_next pour le prochain?
                // dim_angles_idx[dim_angles_idx[high_next]];
                float low_next_angle = dim_angles[dim_angles_idx[low_next]];
                float high_next_angle = dim_angles[dim_angles_idx[high_next]];

                // Distance from 90deg of increasing the low index
                float cost_low = fabs(high_angle - low_next_angle - 90);
                // Distance from 90deg of increasing the high index
                float cost_high = fabs(high_next_angle - low_angle - 90);

                // Angle differences are in mod 360
                cost_low = fmod(cost_low, 360);
                cost_high = fmod(cost_high, 360);

                cost_low = fmin(cost_low, 360 - cost_low);
                cost_high = fmin(cost_high, 360 - cost_high);

                if(cost_low < cost_high) {
                    /* printf("%f < %f\n", cost_low, cost_high); */
                    low_idx_pos = low_next;
                } else {
                    /* printf("%f >= %f\n", cost_low, cost_high); */
                    high_idx_pos = high_next;
                }
            }

            assert(nb_pairs >= n_best_dimensions);

            assert(nb_pairs <= max_nb_pairs);

            qsort_r(paires_scores_idx, nb_pairs, sizeof(uint16_t), compare_idx_with_ref_desc, paires_scores);

            /* for(int i=0; i<2 * nb_bits; i++) { */
            /*     printf("%d, %d\n", i, paires_scores_idx[i]); */
            /* } */
            /* int best_idx = paires_scores_idx[0]; */
            /* for(int i=0; i<nb_pairs - 1; i++) { */
            /*     printf("#% 3d %f >= %f\n", i, paires_scores[paires_scores_idx[i]], paires_scores[paires_scores_idx[i + 1]]); */
            /*     // assert(paires_scores[paires_scores_idx[i]] >= paires_scores[paires_scores_idx[i + 1]]); */
            /* } */
            /* printf("Best pair : [%d %d] %f\n", */
            /*        all_pairs[best_idx * 2], all_pairs[best_idx * 2 + 1], */
            /*        paires_scores[best_idx] */
            /*     ); */

            // TODO : Forcer de la diversité dans les paires, aka,
            // s'assurer qu'on a des dimensions indépendantes le plus
            // possible

            // FIXME : la même paire de dimensions peut revenir deux fois...
            // (107 108) ( 68  69) ( 68  69) (185 186) (187 188) (186 187) (184 185) (117 118) (112 113) ( 84  85)

            // Write to output[...]
            for(int pair_num=0; pair_num<n_best_dimensions; pair_num++) {
                int idx = paires_scores_idx[pair_num];

                int dim0 = all_pairs[idx * 2];
                int dim1 = all_pairs[idx * 2 + 1];

                size_t output_idx =
                    i * (img_w - 1) * n_best_dimensions * 2 +
                    j * n_best_dimensions * 2 +
                    pair_num * 2;

                /* printf("%d %d*%d %d)\n", i, j, dim0, dim1); */

                output[output_idx] = dim0;
                output[output_idx + 1] = dim1;
            }
        }

        free(dim_error);
        free(dim_error_idx);

        free(dim_slope);
        free(dim_slope_idx);

        free(dim_angles);
        free(dim_angles_idx);

        free(considered_dimensions);

        free(all_pairs);

        free(paires_scores);
        free(paires_scores_idx);
    }

    printf("Writing to file\n");

    // Write to file once everything has been computed in parallel
    FILE *fp = fopen(fname, "wb");

    // Store number of best dimensions
    fwrite(&n_best_dimensions, sizeof(uint16_t), 1, fp);

    // Write everything
    fwrite(output, sizeof(uint16_t), output_size, fp);

    fclose(fp);

    return EXIT_SUCCESS;
}
