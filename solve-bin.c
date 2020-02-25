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

void hash_to_codes(mpz_t** to_codes, int* hash_table[2], int* bits_used,
                   int i, int j, int* nb_collisions) {
    size_t hash = 0;

    for(int k=0; k < nb_values; k++) {
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

int spatial_propagation2(float*** matches, mpz_t** from_codes, mpz_t** to_codes) {

    int nb_better_matches = 0;

    for(int iteration=0; iteration<2; iteration++) {

        int delta = iteration * 2 - 1; // -1 or +1

        int i_start = 0, i_end = from_h - 1, i_inc = +1;
        int j_start = 0, j_end = from_w - 1, j_inc = +1;

        if(iteration == 0) {
            i_start = from_h - 1;
            i_end = 0;
            i_inc = -1;

            j_start = from_w - 1;
            j_end = 0;
            j_inc = -1;
        }

        for(int i=i_start; i != i_end; i += i_inc)
            for(int j=j_start; j != j_end; j += j_inc) {

                // Skip thresholded pixels
                if(mask_threshold != -1 && mask[i][j] <= mask_threshold)
                    continue;

                // Process current
                int match_x = matches[X][i][j];

                // Don't try to propagate unassigned pixels
                if(match_x == -1)
                    continue;

                int match_y = matches[Y][i][j];
                // float current_cost = matches[DIST][i][j];

                // Don't propagate outside the images
                if(!(j + delta < from_w && match_x + delta < to_w &&
                     j + delta >= 0 && match_x + delta >= 0 &&
                     i + delta < from_h && match_y + delta < to_h &&
                     i + delta >= 0 && match_y + delta >= 0)) {
                    continue;
                }

                // Propagate horizontally
                if(mask_threshold == -1 || mask[i][j + delta] > mask_threshold) {

                    // Cost of copying Current <-> Neighbor
                    int copy_cost = mpz_hamdist(
                        from_codes[i][j],
                        to_codes[match_y][match_x + delta]);

                    // Cost of propagating Neighbor <-> Neighbor
                    int propagation_cost = mpz_hamdist(
                        from_codes[i][j + delta],
                        to_codes[match_y][match_x + delta]);

                    int to_x = match_x + delta;
                    int to_y = match_y;
                    int to_cost = propagation_cost;

                    // Choose best of copy/propagation
                    if(copy_cost < propagation_cost) {
                        to_x = match_x;
                        to_cost = copy_cost;
                    }

                    // Current match cost
                    float prev_neighbour_cost = matches[DIST][i][j + delta];

                    // If the neighbor is unmatched or is already
                    // matched with a higher cost, replace it
                    if(prev_neighbour_cost == -1.0 || to_cost < prev_neighbour_cost) {
                        nb_better_matches++;

                        matches[X][i][j + delta] = to_x;
                        matches[Y][i][j + delta] = to_y;
                        matches[DIST][i][j + delta] = to_cost;
                    }

                    // FIXME : probablement pas nécessaire, ça va se faire à une étape ultérieure anyway
                    // OU... Possiblement utile, puisqu'on a déjà les données en mémoire, ça évite un calcul?
                    /* // Possibly update actual match */
                    /* if(copy_cost < current_cost) { */
                    /*     matches[X][i][j] = match_x + delta; */
                    /*     matches[Y][i][j] = match_y; */
                    /*     matches[DIST][i][j] = copy_cost; */
                    /*     current_cost = copy_cost; */
                    /* } */
                }

                // Propagate vertically
                if(mask_threshold == -1 || mask[i + delta][j] > mask_threshold) {

                    // Match by copying Current <-> Neighbor
                    int copy_cost = mpz_hamdist(
                        from_codes[i][j],
                        to_codes[match_y + delta][match_x]);

                    // Match by propagating Neighbor <-> Neighbor
                    int propagation_cost = mpz_hamdist(
                        from_codes[i + delta][j],
                        to_codes[match_y + delta][match_x]);

                    int to_x = match_x;
                    int to_y = match_y + delta;
                    int to_cost = propagation_cost;

                    // Choose best of copy/propagation
                    if(copy_cost < propagation_cost) {
                        to_y = match_y;
                        to_cost = copy_cost;
                    }

                    // Current match cost
                    float prev_neighbour_cost = matches[DIST][i + delta][j];

                    // If the neighbor is unmatched or is already
                    // matched with a higher cost, replace it
                    if(prev_neighbour_cost == -1.0 || to_cost < prev_neighbour_cost) {
                        nb_better_matches++;

                        matches[X][i + delta][j] = to_x;
                        matches[Y][i + delta][j] = to_y;
                        matches[DIST][i + delta][j] = to_cost;
                    }

                    // FIXME : probablement pas nécessaire, ça va se faire à une étape ultérieure anyway
                    // OU... Possiblement utile, puisqu'on a déjà les données en mémoire, ça évite un calcul?
                    /* // Possibly update actual match */
                    /* if(copy_cost < current_cost) { */
                    /*     matches[X][i][j] = match_x; */
                    /*     matches[Y][i][j] = match_y + delta; */
                    /*     matches[DIST][i][j] = copy_cost; */
                    /* } */
                }
            }
    }
    printf("%d\n", nb_better_matches);
    return nb_better_matches;
}

int spatial_propagation(float*** matches, mpz_t** from_codes, mpz_t** to_codes) {

    int nb_better_matches = 0;

    for(int iteration=0; iteration<2; iteration++) {

        int delta = iteration * 2 - 1; // -1 or +1

        int i_start = 0, i_end = from_h - 1, i_inc = +1;
        int j_start = 0, j_end = from_w - 1, j_inc = +1;

        if(iteration == 0) {
            i_start = from_h - 1;
            i_end = 0;
            i_inc = -1;

            j_start = from_w - 1;
            j_end = 0;
            j_inc = -1;
        }

        for(int i=i_start; i != i_end; i += i_inc)
            for(int j=j_start; j != j_end; j += j_inc) {

                // Skip thresholded pixels
                if(mask_threshold != -1 && mask[i][j] <= mask_threshold)
                    continue;

                // Process current
                int match_x = matches[X][i][j];

                // Don't try to propagate unassigned pixels
                if(match_x == -1)
                    continue;

                int match_y = matches[Y][i][j];
                // float current_cost = matches[DIST][i][j];

                // Don't propagate outside the images
                if(!(j + delta < from_w && match_x + delta < to_w &&
                     j + delta >= 0 && match_x + delta >= 0 &&
                     i + delta < from_h && match_y + delta < to_h &&
                     i + delta >= 0 && match_y + delta >= 0)) {
                    continue;
                }

                // Propagate horizontally
                if(mask_threshold == -1 || mask[i][j + delta] > mask_threshold) {

                    // Cost of copying Neighbor <-> Current
                    int copy_cost = mpz_hamdist(
                        from_codes[i][j + delta],
                        to_codes[match_y][match_x]);

                    // Cost of propagating Neighbor <-> Neighbor
                    int propagation_cost = mpz_hamdist(
                        from_codes[i][j + delta],
                        to_codes[match_y][match_x + delta]);

                    int to_x = match_x + delta;
                    int to_y = match_y;
                    int to_cost = propagation_cost;

                    // Choose best of copy/propagation
                    if(copy_cost < propagation_cost) {
                        to_x = match_x;
                        to_cost = copy_cost;
                    }

                    // Current match cost
                    float prev_neighbour_cost = matches[DIST][i][j + delta];

                    // If the neighbor is unmatched or is already
                    // matched with a higher cost, replace it
                    if(prev_neighbour_cost == -1.0 || to_cost < prev_neighbour_cost) {
                        nb_better_matches++;

                        matches[X][i][j + delta] = to_x;
                        matches[Y][i][j + delta] = to_y;
                        matches[DIST][i][j + delta] = to_cost;
                    }

                    // FIXME : probablement pas nécessaire, ça va se faire à une étape ultérieure anyway
                    // OU... Possiblement utile, puisqu'on a déjà les données en mémoire, ça évite un calcul?
                    /* // Possibly update actual match */
                    /* if(copy_cost < current_cost) { */
                    /*     matches[X][i][j] = match_x + delta; */
                    /*     matches[Y][i][j] = match_y; */
                    /*     matches[DIST][i][j] = copy_cost; */
                    /*     current_cost = copy_cost; */
                    /* } */
                }

                // Propagate vertically
                if(mask_threshold == -1 || mask[i + delta][j] > mask_threshold) {

                    // Match by copying Neighbor <-> Current
                    int copy_cost = mpz_hamdist(
                        from_codes[i + delta][j],
                        to_codes[match_y][match_x]);

                    // Match by propagating Neighbor <-> Neighbor
                    int propagation_cost = mpz_hamdist(
                        from_codes[i + delta][j],
                        to_codes[match_y + delta][match_x]);

                    int to_x = match_x;
                    int to_y = match_y + delta;
                    int to_cost = propagation_cost;

                    // Choose best of copy/propagation
                    if(copy_cost < propagation_cost) {
                        to_y = match_y;
                        to_cost = copy_cost;
                    }

                    // Current match cost
                    float prev_neighbour_cost = matches[DIST][i + delta][j];

                    // If the neighbor is unmatched or is already
                    // matched with a higher cost, replace it
                    if(prev_neighbour_cost == -1.0 || to_cost < prev_neighbour_cost) {
                        nb_better_matches++;

                        matches[X][i + delta][j] = to_x;
                        matches[Y][i + delta][j] = to_y;
                        matches[DIST][i + delta][j] = to_cost;
                    }

                    // FIXME : probablement pas nécessaire, ça va se faire à une étape ultérieure anyway
                    // OU... Possiblement utile, puisqu'on a déjà les données en mémoire, ça évite un calcul?
                    /* // Possibly update actual match */
                    /* if(copy_cost < current_cost) { */
                    /*     matches[X][i][j] = match_x; */
                    /*     matches[Y][i][j] = match_y + delta; */
                    /*     matches[DIST][i][j] = copy_cost; */
                    /* } */
                }
            }
    }
    printf("%d\n", nb_better_matches);
    return nb_better_matches;
}

int cont_spatial_propagation(float*** matches, float*** from_codes, float*** to_codes) {

    int nb_better_matches = 0;

    for(int iteration=0; iteration<2; iteration++) {

        int delta = iteration * 2 - 1; // -1 or +1

        int i_start = 0, i_end = from_h - 1, i_inc = +1;
        int j_start = 0, j_end = from_w - 1, j_inc = +1;

        if(iteration == 0) {
            i_start = from_h - 1;
            i_end = 0;
            i_inc = -1;

            j_start = from_w - 1;
            j_end = 0;
            j_inc = -1;
        }

        for(int i=i_start; i != i_end; i += i_inc)
            for(int j=j_start; j != j_end; j += j_inc) {

                // Skip thresholded pixels
                if(mask_threshold != -1 && mask[i][j] <= mask_threshold)
                    continue;

                // Process current
                int match_x = matches[X][i][j];

                // Don't try to propagate unassigned pixels
                if(match_x == -1)
                    continue;

                int match_y = matches[Y][i][j];
                // float current_cost = matches[DIST][i][j];

                // Don't propagate outside the images
                if(!(j + delta < from_w && match_x + delta < to_w &&
                     j + delta >= 0 && match_x + delta >= 0 &&
                     i + delta < from_h && match_y + delta < to_h &&
                     i + delta >= 0 && match_y + delta >= 0)) {
                    continue;
                }

                // Propagate horizontally
                if(mask_threshold == -1 || mask[i][j + delta] > mask_threshold) {

                    // Match by copying Neighbor <-> Current
                    float copy_cost = pixel_cost(
                        from_codes[i][j + delta],
                        to_codes[match_y][match_x],
                        nb_bits);

                    // Match by propagating Neighbor <-> Neighbor
                    float propagation_cost = pixel_cost(
                        from_codes[i][j + delta],
                        to_codes[match_y][match_x + delta],
                        nb_bits);

                    int to_x = match_x + delta;
                    int to_y = match_y;
                    float to_cost = propagation_cost;

                    // Choose best of copy/propagation
                    if(copy_cost < propagation_cost) {
                        to_x = match_x;
                        to_cost = copy_cost;
                    }

                    int neighbor_match_x = (int) matches[X][i][j + delta];
                    int neighbor_match_y = (int) matches[Y][i][j + delta];

                    // Current match cost pour le voisin
                    /* float prev_neighbour_cost = matches[DIST][i][j + delta]; // TODO : wrong */
                    // Old neighbor match cost
                    float prev_neighbour_cost = -1;

                    if(neighbor_match_x != -1) {
                        prev_neighbour_cost = pixel_cost(
                            from_codes[i][j + delta],
                            to_codes[neighbor_match_y][neighbor_match_x],
                            nb_bits);
                    }

                    // If the neighbor is unmatched or is already
                    // matched with a higher cost, replace it
                    if(neighbor_match_x == -1 || to_cost < prev_neighbour_cost) {
                        nb_better_matches++;

                        matches[X][i][j + delta] = to_x;
                        matches[Y][i][j + delta] = to_y;
                        matches[DIST][i][j + delta] = nb_bits - 1; /* to_cost; */ // FIXME
                    }

                    // FIXME : probablement pas nécessaire, ça va se faire à une étape ultérieure anyway
                    // OU... Possiblement utile, puisqu'on a déjà les données en mémoire, ça évite un calcul?
                    /* // Possibly update actual match */
                    /* if(copy_cost < current_cost) { */
                    /*     matches[X][i][j] = match_x + delta; */
                    /*     matches[Y][i][j] = match_y; */
                    /*     matches[DIST][i][j] = copy_cost; */
                    /*     current_cost = copy_cost; */
                    /* } */
                }

                // Propagate vertically
                if(mask_threshold == -1 || mask[i + delta][j] > mask_threshold) {

                    // Match by copying Neighbor <-> Current
                    float copy_cost = pixel_cost(
                        from_codes[i + delta][j],
                        to_codes[match_y][match_x],
                        nb_bits);

                    // Match by propagating Neighbor <-> Neighbor
                    float propagation_cost = pixel_cost(
                        from_codes[i + delta][j],
                        to_codes[match_y + delta][match_x],
                        nb_bits);

                    int to_x = match_x;
                    int to_y = match_y + delta;
                    float to_cost = propagation_cost;

                    // Choose best of copy/propagation
                    if(copy_cost < propagation_cost) {
                        to_y = match_y;
                        to_cost = copy_cost;
                    }

                    int neighbor_match_x = (int) matches[X][i + delta][j];
                    int neighbor_match_y = (int) matches[Y][i + delta][j];

                    // Current match cost pour le voisin
                    /* float prev_neighbour_cost = matches[DIST][i][j + delta]; // TODO : wrong */
                    // Old neighbor match cost
                    float prev_neighbour_cost = -1;

                    if(neighbor_match_y != -1) {
                        prev_neighbour_cost = pixel_cost(
                            from_codes[i + delta][j],
                            to_codes[neighbor_match_y][neighbor_match_x],
                            nb_bits);
                    }

                    // If the neighbor is unmatched or is already
                    // matched with a higher cost, replace it
                    if(neighbor_match_y == -1 || to_cost < prev_neighbour_cost) {

                        nb_better_matches++;

                        matches[X][i + delta][j] = to_x;
                        matches[Y][i + delta][j] = to_y;
                        matches[DIST][i + delta][j] = nb_bits - 1; /* to_cost; */ // FIXME
                    }

                    // FIXME : probablement pas nécessaire, ça va se faire à une étape ultérieure anyway
                    // OU... Possiblement utile, puisqu'on a déjà les données en mémoire, ça évite un calcul?
                    /* // Possibly update actual match */
                    /* if(copy_cost < current_cost) { */
                    /*     matches[X][i][j] = match_x; */
                    /*     matches[Y][i][j] = match_y + delta; */
                    /*     matches[DIST][i][j] = copy_cost; */
                    /* } */
                }
            }
    }
    printf("Continuous spatial better matches : %d\n", nb_better_matches);
    return nb_better_matches;
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

        int distance = mpz_hamdist(from_codes[i][j], to_codes[y][x]);

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

    free(hash_table[0]);
    free(hash_table[1]);
    free(bits_used);

    printf("nb_new_matches: %d\n", nb_new_matches);
    printf("nb_better_matches: %d\n", nb_better_matches);
}

int main(int argc, char** argv) {

    int i, j, k, l;

    int nthreads = 4, max_iterations = 30,
        disable_heuristics = 0, proj_lut = 0, dump_all_images = 0;

    int max_heuristic_iterations = 5, use_cont_heuristic = 0;

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

    LARG_CASE("use-cont-heuristic")
        use_cont_heuristic = 1;

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

    mpz_t** cam_codes;
    mpz_t** ref_codes;

    if(quadratic) {
        nb_bits = nb_patterns * (nb_patterns - 1) / 2;

        printf("Loading cam binary codes\n");
        cam_codes = load_binary_codesq(cam_format, from_w, from_h, nb_patterns);

        printf("Loading ref binary codes\n");
        ref_codes = load_binary_codesq(ref_format, to_w, to_h, nb_patterns);
    } else {
        printf("Loading cam binary codes\n");
        cam_codes = load_binary_codesl(cam_format, from_w, from_h, nb_patterns);

        printf("Loading ref binary codes\n");
        ref_codes = load_binary_codesl(ref_format, to_w, to_h, nb_patterns);
    }

    if(mask_threshold != -1) {
        mask = load_mask(proj_lut ? ref_format : cam_format, nb_patterns, from_w, from_h);
    }

    time_since_start = time(NULL);

    for(l=0; l<max_iterations; l++) {
        printf("----- Iteration %02d -----\n", l);
        lsh(matches, cam_codes, ref_codes);

        // Save the iteration
        if(dump_all_images) {
            sprintf(filename, out_format, l);

            save_color_map(filename, matches,
                           from_w, from_h, to_w, to_h, nb_bits);
        }

        printf("delta time = %ld\n", time(NULL) - time_since_start);
    }

    if(!disable_heuristics) {

        save_color_map("lutCam-no-heuristic.png", matches,
                       from_w, from_h, to_w, to_h, nb_bits);

        // Spatial propagation heuristic with binary codes
        for(int iteration=0; iteration<max_heuristic_iterations; iteration++) {
            int nb_better = 1;

            nb_better = spatial_propagation(matches, cam_codes, ref_codes);

            // Save the iteration
            if(dump_all_images || nb_better == 0 || iteration == max_heuristic_iterations - 1) {
                sprintf(filename, "lutCam%02d-bin-flood.png", iteration);

                save_color_map(filename, matches,
                               from_w, from_h, to_w, to_h, nb_bits);
            }

            printf("delta time = %ld\n", time(NULL) - time_since_start);

            if(nb_better == 0)
                break;
        }

        // TODO : second batch of heuristics with float code instead
        free_mpz_matrix(cam_codes);
        free_mpz_matrix(ref_codes);

        if(use_cont_heuristic) {
            assert(quadratic);

            printf("Loading cam continuous codes\n");
            float*** cont_cam_codes = load_continuous_codesq(cam_format, from_w, from_h, nb_patterns);
            printf("Loading ref continuous codes\n");
            float*** cont_ref_codes = load_continuous_codesq(ref_format, to_w, to_h, nb_patterns);

            for(int iteration=0; iteration<max_heuristic_iterations; iteration++) {
                int nb_better = 1;

                nb_better = cont_spatial_propagation(matches, cont_cam_codes, cont_ref_codes);

                // Save the iteration
                if(dump_all_images || nb_better == 0 || iteration == max_heuristic_iterations - 1) {
                    sprintf(filename, "lutCam%02d-cont-flood.png", iteration);

                    save_color_map(filename, matches,
                                   from_w, from_h, to_w, to_h, nb_bits);
                }

                printf("delta time = %ld\n", time(NULL) - time_since_start);

                if(nb_better == 0)
                    break;
            }
        }
    }

    save_color_map(
        proj_lut ? "lutProjPixel.png" : "lutCamPixel.png",
        matches, from_w, from_h, to_w, to_h, nb_bits);

    return EXIT_SUCCESS;
}
