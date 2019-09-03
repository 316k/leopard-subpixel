/*
  Compute camera-projector LUTs from matching patterns/pictures
*/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#include "args.h"
#include "helpers.c"

mpz_t** malloc_f32matrixul(int w, int h, int nb_bits) {
    mpz_t** matrix = malloc(sizeof(mpz_t*) * h);
    matrix[0] = (mpz_t*) calloc(h * w, sizeof(mpz_t));

    for(int i=0; i < h; i++) {
        matrix[i] = (*matrix + w * i);

        for(int j=0; j<w; j++) {
            mpz_init2(matrix[i][j], nb_bits);
        }
    }

    return matrix;
}


void free_f32matrixul(mpz_t** matrix) {
    free(matrix[0]);
    free(matrix);
}


// --- Global variables ---
// Input-related parameters
int nb_patterns, from_w, from_h, to_w, to_h, skip = 1, nb_bits;

// Algorithm-related parameters
int nb_values = -1;

// Should be about 3*sigmas of the normal distribution error
float* offsets;
/* float noise_error = 0.1; */
// float quantizing_range = 0.4;
// int boxes_mid, boxes_range;

size_t hash_table_size;
// float hash_divider;

// Threshold to consider a pixel as part of the scanned zone
// -1 == disable mask, hash all pixels
float mask_threshold = 20;
float** mask;

long time_since_start;

#define HASH_TO_CODES   hash_to_codes(to_codes, hash_table, bits_used, i, j, &nb_collisions)

#define HASH_FROM_CODES hash_from_codes(matches, from_codes, to_codes, \
                                          hash_table, bits_used, i, j, \
                                          &nb_new_matches, &nb_better_matches)

int bitCount(unsigned long n) {
    n = ((0xaaaaaaaaaaaaaaaaL & n) >>  1) + (0x5555555555555555L & n);
    n = ((0xccccccccccccccccL & n) >>  2) + (0x3333333333333333L & n);
    n = ((0xf0f0f0f0f0f0f0f0L & n) >>  4) + (0x0f0f0f0f0f0f0f0fL & n);
    n = ((0xff00ff00ff00ff00L & n) >>  8) + (0x00ff00ff00ff00ffL & n);
    n = ((0xffff0000ffff0000L & n) >> 16) + (0x0000ffff0000ffffL & n);
    n = ((0xffffffff00000000L & n) >> 32) + (0x00000000ffffffffL & n);
    return (int) n;
}

void hash_to_codes(mpz_t** to_codes, int* hash_table[2], int* bits_used,
                   int i, int j, int* nb_collisions) {
    size_t hash = 0;

    for(int k=0; k < nb_values; k++) {
        // int bit = 1 << bits_used[k];
        // hash |= (!!(to_codes[i][j] & bit)) << k;

        hash |= mpz_tstbit(to_codes[i][j], bits_used[k]) << k;
    }

    // Premier arrivé prend la place
    if(hash_table[X][hash] == -1) {
        hash_table[X][hash] = j;
        hash_table[Y][hash] = i;
    } else {
        (*nb_collisions)++;
    }
}

/**
 * Attempts to improve an existing match by comparing the neighborhood of
 * the matched *to* pixel with the matched *from* pixel
 */
void forward_matching(float*** matches, float*** from_codes, float*** to_codes,
                      float* from_code, int from_x, int from_y, int to_x, int to_y) {

    float* to_code = malloc(sizeof(float) * nb_bits);

    for(int i=fmax(to_y - 1, 0); i < fmin(to_y + 1, to_h - 1); i++) // <= ?
        for(int j=fmax(to_x - 1, 0); j < fmin(to_x + 1, to_w - 1); j++) {

            for(int k=0; k < nb_bits; k++) {
                to_code[k] = to_codes[k][i][j];
            }

            float distance = distance_modulo_pi(from_code, to_code, nb_bits);

            // When new distance is lower, update the match
            if(distance < matches[DIST][from_y][from_x]) {
                matches[X][from_y][from_x]    = j;
                matches[Y][from_y][from_x]    = i;
                matches[DIST][from_y][from_x] = distance;

                assert(matches[DIST][from_y][from_x] > 0);
            }
        }

    free(to_code);
}

/**
 * Attempts to create new matches by looking in the neighborhood of
 * the matched *from* pixel
 */
void backward_matching(float*** matches, float*** from_codes, float*** to_codes,
                       float* to_code, int from_x, int from_y,
                       int to_x, int to_y) {

    float* from_code = malloc(sizeof(float) * nb_bits);

    for(int i=fmax(from_y - 1, 0); i < fmin(from_y + 1, from_h - 1); i++) // <= ?
        for(int j=fmax(from_x - 1, 0); j < fmin(from_x + 1, from_w - 1); j++) {

            // Skip {current pixel -> current match}
            if(i == from_y && j == from_x)
                continue;

            for(int k=0; k < nb_bits; k++) {
                from_code[k] = from_codes[k][i][j];
            }

            float distance = distance_modulo_pi(from_code, to_code, nb_bits);

            // When it's a new match or new distance is lower, update the match
            if(matches[DIST][i][j] == -1.0 || distance < matches[DIST][i][j]) {
                matches[X][i][j]    = to_x;
                matches[Y][i][j]    = to_y;
                matches[DIST][i][i] = distance;
                assert(matches[DIST][from_y][from_x] > 0);
            }
        }

    free(from_code);
}

void hash_from_codes(float*** matches, mpz_t** from_codes, mpz_t** to_codes,
                     int* hash_table[2], int* bits_used, int i, int j,
                     int* nb_new_matches, int* nb_better_matches) {

    size_t hash = 0;

    for(int k=0; k < nb_values; k++) {

        hash |= mpz_tstbit(from_codes[i][j], bits_used[k]) << k;
    }

    // Collision = match
    if(hash_table[X][hash] != -1) {
        int x = hash_table[X][hash];
        int y = hash_table[Y][hash];

        int distance = mpz_hamdist(from_codes[i][j], to_codes[y][x]); // bitCount(from_code ^ to_code);

        // Si la nouvelle distance est plus petite, on update le match
        if(matches[DIST][i][j] == -1.0 || distance < matches[DIST][i][j]) {
            if(matches[DIST][i][j] == -1.0)
                (*nb_new_matches)++;
            else
                (*nb_better_matches)++;

            matches[X][i][j] = x;
            matches[Y][i][j] = y;
            matches[DIST][i][j] = distance;
        }
    }
}

/**
 * Peforms one iteration of the LSH algorithm.
 *
 * This function is NOT thread-safe and should NOT be called multiple
 * times in parallel. All the parallel work is done from inside the
 * function.
 */
void lsh(float*** matches, mpz_t** from_codes, mpz_t** to_codes) {

    int i, j;
    int nb_collisions = 0, nb_new_matches = 0, nb_better_matches = 0;

    hash_table_size = (size_t) 1 << nb_values;

    printf("nb_values=%d\n", nb_values);

    printf("hash_table_size=%d (%f times bigger than %d)\n",
           (int) hash_table_size, hash_table_size /(float)(to_w * to_h),
           to_w * to_h);

    int* hash_table[2];

    for(i=0; i<2; i++) {
        hash_table[i] = malloc(hash_table_size * sizeof(int));

        for(j=0; j<hash_table_size; j++)
            hash_table[i][j] = -1;
    }

    int* bits_used = random_phases(nb_values, nb_bits);

    int iteration = rand() % 4;

    // Hash des to_codes
    switch(iteration) {
    case 0:
        #pragma omp parallel for private(i, j)
        for(i=0; i < to_h; i++)
            for(j=0; j < to_w; j++)
                HASH_TO_CODES;
        break;

    case 1:
        #pragma omp parallel for private(i, j)
        for(i=0; i < to_h; i++)
            for(j=to_w - 1; j >= 0; j--)
                HASH_TO_CODES;
        break;

    case 2:
        #pragma omp parallel for private(i, j)
        for(i=to_h - 1; i >= 0; i--)
            for(j=0; j < to_w; j++)
                HASH_TO_CODES;
        break;

    case 3:
        #pragma omp parallel for private(i, j)
        for(i=to_h - 1; i >= 0; i--)
            for(j=to_w - 1; j >= 0; j--)
                HASH_TO_CODES;
        break;
    }


    printf("nb_collisions: %d / %d = %f%%\n", nb_collisions, to_w * to_h,
           nb_collisions /(float)(to_w * to_h) * 100);

    printf("hash_space_used: %d / %d = %f%%\n", to_w * to_h - nb_collisions,
           (int) hash_table_size,
           (to_w * to_h - nb_collisions)/(float)hash_table_size * 100);

    // from_codes hashing (+ updating matches)
    switch(iteration) {
    case 0:
        #pragma omp parallel for private(i, j)
        for(i=0; i < from_h; i++)
            for(j=0; j < from_w; j++)
                if(mask_threshold == -1 || mask[i][j] > mask_threshold)
                    HASH_FROM_CODES;
        break;

    case 1:
        #pragma omp parallel for private(i, j)
        for(i=0; i < from_h; i++)
            for(j=from_w - 1; j >= 0; j--)
                if(mask_threshold == -1 || mask[i][j] > mask_threshold)
                    HASH_FROM_CODES;
        break;

    case 2:
        #pragma omp parallel for private(i, j)
        for(i=from_h - 1; i >= 0; i--)
            for(j=0; j < from_w; j++)
                if(mask_threshold == -1 || mask[i][j] > mask_threshold)
                    HASH_FROM_CODES;
        break;

    case 3:
        #pragma omp parallel for private(i, j)
        for(i=from_h - 1; i >= 0; i--)
            for(j=from_w - 1; j >= 0; j--)
                if(mask_threshold == -1 || mask[i][j] > mask_threshold)
                    HASH_FROM_CODES;
        break;
    }

    /* // Dump hash-table profile */
    /* int idx = 0, count = 0, empty_run = hash_table[X][i] == -1; */
    /* int nb_empty_spaces = 0; */

    /* for(i=1; i<=hash_table_size; i++) { */
    /*     if(i == hash_table_size || ( */
    /*            // Whenever the run type changes */
    /*            empty_run ? (hash_table[X][i] != -1) : (hash_table[X][i] == -1) */
    /*            )) { */
    /*         printf("%s %08d -> %04d\n", empty_run ? "[ ]" : "[x]", idx, count); */
    /*         // Reset run settings */
    /*         empty_run = hash_table[X][i] == -1; */
    /*         count = 0; */
    /*         idx = i; */
    /*     } else { */
    /*         count++; */
    /*     } */

    /*     nb_empty_spaces += (hash_table[X][i] == -1); */
    /* } */

    /* printf("nb_empty_spaces=%d\n", nb_empty_spaces); */

    free(hash_table[0]);
    free(hash_table[1]);
    free(bits_used);

    printf("nb_new_matches: %d\n", nb_new_matches);
    printf("nb_better_matches: %d\n", nb_better_matches);
}

int main(int argc, char** argv) {

    int i, j, k, l, n;

    int nthreads = 4, max_iterations = 30,
        disable_heuristics = 0, proj_lut = 0, dump_all_images = 0;

    int use_default_out_format = 1;

    int quadratic = 0;

    /*
      ~~~ TODO ~~~

      Stop the iterations after the condition is met :
      |delta(nb_better_matches + nb_new_matches)| < stop_threshold
      For N=stop_threshold_nb_passes iterations
     */
    float stop_threshold = 0.1;
    int stop_threshold_nb_passes = 3;

    char* ref_format = "leo_%03d.png";
    char* cam_format = "%03d.png";
    char* out_format = "lutCam%02d.png";

    char filename[FNAME_MAX_LEN]; // Generic filename buffer

    // Args parsing
    ARGBEGIN
        ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('p')
        proj_lut = 1;
        mask_threshold = -1;

    ARG_CASE('R')
        ref_format = ARGS;

    ARG_CASE('C')
        cam_format = ARGS;

    ARG_CASE('O')
        out_format = ARGS;
        use_default_out_format = 0;

    ARG_CASE('i')
        max_iterations = ARGI;

    ARG_CASE('d')
        disable_heuristics = 1;

    ARG_CASE('v')
        nb_values = ARGI;

    ARG_CASE('m')
        if(proj_lut)
            fprintf(stderr, "*** WARNING : -m has no effect when -p is enabled\n");
        else
            mask_threshold = ARGF;

    ARG_CASE('q')
        quadratic = 1;

    LARG_CASE("dump-all-images")
        dump_all_images = 1;

    WRONG_ARG
        printf("usage: %s [-t nb_threads=%d] [-p gen_proj_lut]\n"
               "\t[-R ref_format=\"%s\"] [-C cam_format=\"%s\"]\n"
               "\t[-O out_format=\"%s\"]\n"
               "\t[-i max_iterations=%d] [-d (disable heuristics)]\n"
               /* "\t[-b nb_boxes=%d]" */
               "\t[-v nb_values=%d]\n"
               /* "\t[-r quantizing_range=%0.2f]\n" */
               /* "\t[-r hash_table_ratio=%f]\n" */
               "\t[-m mask_threshold=%0.2f] [--dump-all-images]\n\n"
               "\t[-s stop_threshold=%0.2f] [--stop_threshold_nb_passes=%d] (TODO)\n"
               /* "\tquantizing_range should be between 0 and 1\n" */,
               argv0, nthreads, ref_format, cam_format, out_format, max_iterations,
               /* nb_boxes, */ nb_values, /* quantizing_range, */ /* hash_table_ratio,  */
               mask_threshold, stop_threshold, stop_threshold_nb_passes
            );
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


    int nb_pixels_max = fmax(to_w * to_h, from_w * from_h);

    if(nb_values == -1) {
        nb_values = (int) ceil(log2f(nb_pixels_max));
    }

    // Generating projector LUT => read captured images as projected
    // images
    if(proj_lut) {
        char* tmp_format = ref_format;
        ref_format = cam_format;
        cam_format = tmp_format;

        if(use_default_out_format)
            out_format = "lutProj%02d.png";

        int tmp;
        tmp = from_w;
        from_w = to_w;
        to_w = tmp;

        tmp = from_h;
        from_h = to_h;
        to_h = tmp;
    }

    float*** matches = malloc_f32cube(3, from_w, from_h); // matches[ x, y, distance ][h=i][w=j]

    #pragma omp parallel for private(i, j)
    for(k=0; k<3; k++)
        for(i=0; i<from_h; i++)
            for(j=0; j<from_w; j++)
                matches[k][i][j] = -1.0;

    nb_bits = nb_patterns - 1;

    if(quadratic)
        nb_bits = nb_patterns * (nb_patterns - 1) / 2;

    mpz_t** cam_codes = malloc_f32matrixul(from_w, from_h, nb_bits);
    mpz_t** ref_codes = malloc_f32matrixul(to_w, to_h, nb_bits);

    int w, h;

    // --- Load binary codes ---

    printf("Loading cam binary codes\n");

    // Camera binary codes
    for(k=0, n=0; k<nb_patterns - 1; k++) {

        // Load first cam image
        float** previous_image;

        sprintf(filename, cam_format, k);
        previous_image = load_gray(filename, &w, &h);

        float** current_image;

        // Load next images
        for(l=k + 1; l<nb_patterns; l++, n++) {

            sprintf(filename, cam_format, l);
            current_image = load_gray(filename, &w, &h);

            for(i=0; i<from_h; i++)
                for(j=0; j<from_w; j++) {

                    int diff = (current_image[i][j] - previous_image[i][j]) < 0;

                    if(diff) {
                        mpz_setbit(cam_codes[i][j], n);
                    }
                }

            free_f32matrix(current_image);

            if(!quadratic) {
                n++;
                break;
            }
        }

        // Free last image
        free_f32matrix(previous_image);
    }

    // Ref binary codes
    printf("Loading ref binary codes\n");
    for(k=0, n=0; k<nb_patterns - 1; k++) {

        // Load first cam image
        float** previous_image;

        sprintf(filename, ref_format, k);
        previous_image = load_gray(filename, &w, &h);

        float** current_image;

        // Load next images
        for(l=k + 1; l<nb_patterns; l++, n++) {

            sprintf(filename, ref_format, l);
            current_image = load_gray(filename, &w, &h);

            for(i=0; i<from_h; i++)
                for(j=0; j<from_w; j++) {

                    int diff = (current_image[i][j] - previous_image[i][j]) < 0;

                    if(diff)
                        mpz_setbit(ref_codes[i][j], n);
                }

            free_f32matrix(current_image);

            if(!quadratic) {
                n++;
                break;
            }
        }

        // Free last image
        free_f32matrix(previous_image);
    }

    if(mask_threshold != -1) {
        mask = load_mask(proj_lut ? ref_format : cam_format, nb_patterns, from_w, from_h);
    }

    /* mpz_t from_code; */
    /* mpz_t to_code; */

    time_since_start = time(NULL);

    for(l=0; l<max_iterations; l++) {
        printf("----- Iteration %02d -----\n", l);
        lsh(matches, cam_codes, ref_codes);

        // Launch heuristics
        disable_heuristics = 1;
        if(!disable_heuristics && l % 5 == 0 && l != 0) {
            printf("----- Running backward & forward matching -----\n");

            #pragma omp parallel for private(i, j, k)
            for(i=0; i<from_h; i++)
                for(j=0; j<from_w; j++) {

                    if(matches[DIST][i][j] == -1)
                        continue;

                    /* from_code = cam_codes[i][j]; */
                    /* to_code = ref_codes[i][j]; */

                    int x = matches[X][i][j];
                    int y = matches[Y][i][j];

                    /* if(mask_threshold == -1 || mask[i][j] > mask_threshold) { */
                    /*     forward_matching(matches, cam_codes, ref_codes, */
                    /*                      from_code, j, i, x, y); */
                    /* } */

                    /* backward_matching(matches, cam_codes, ref_codes, */
                    /*                   ref_codes[i][j], j, i, x, y); */
                }
        }

        // Save the iteration
        if(dump_all_images || l == max_iterations - 1 || l % 5 == 0) {
            /* if(proj_lut) */
            /*     sprintf(filename, "lutProj%02d.png", l); */
            /* else */
            sprintf(filename, out_format, l);

            save_color_map(filename, matches,
                           from_w, from_h, to_w, to_h, nb_bits);
        }

        printf("delta time = %ld\n", time(NULL) - time_since_start);
    }

    return EXIT_SUCCESS;
}
