#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "args.h"
#include "helpers.c"

// Global variables
int nb_boxes = 10, nb_values = 8, nb_patterns,
    from_w, from_h, to_w, to_h,
    hash_table_size;
float hash_divider;


#define HASH_TO_CODES()   hash_to_codes(to_codes, hash_table, phases_used, i, j, &nb_collisions)

#define HASH_FROM_CODES() hash_from_codes(matches, from_codes, to_codes, \
                                          hash_table, phases_used, i, j, \
                                          use_heuristics, &nb_new_matches, &nb_better_matches)

int hash_fct(float* pixel_code) {
    int hash = 0;

    for(int i=0; i<nb_values; i++) {
        hash *= nb_boxes;
        hash += fmod(pixel_code[i] + PI, 2 * PI) / (2 * PI) * nb_boxes;
    }

    return floor(hash / hash_divider);
}

void hash_to_codes(float*** to_codes, int* hash_table[2], int* phases_used,
                   int i, int j, int* nb_collisions) {
    float* pixel_code = malloc(sizeof(float) * nb_values);

    for(int k=0; k < nb_values; k++) {
        int phase = phases_used[k];

        pixel_code[k] = to_codes[phase][i][j];
    }

    int hash = hash_fct(pixel_code);

    // Premier arrivé prend la place
    if(hash_table[X][hash] == -1) {
        hash_table[X][hash] = j;
        hash_table[Y][hash] = i;
    } else {
        (*nb_collisions)++;
    }

    free(pixel_code);
}

/**
 * Attempts to improve an existing match by comparing the neighborhood of
 * the matched *to* pixel with the matched *from* pixel
 */
void forward_matching(float*** matches, float*** from_codes, float*** to_codes,
                      float* from_code, int from_x, int from_y, int to_x, int to_y) {

    float* to_code = malloc(sizeof(float) * nb_patterns);

    for(int i=fmax(to_y - 1, 0); i < fmin(to_y + 1, to_h - 1); i++)
        for(int j=fmax(to_x - 1, 0); j < fmin(to_x + 1, to_w - 1); j++) {

            for(int k=0; k < nb_patterns; k++) {
                to_code[k] = to_codes[k][i][j];
            }

            float distance = distance_modulo_pi(from_code, to_code, nb_patterns);

            // Si la nouvelle distance est plus petite, on update le match
            if(distance < matches[DIST][from_y][from_x]) {
                matches[X][from_y][from_x]    = j;
                matches[Y][from_y][from_x]    = i;
                matches[DIST][from_y][from_x] = distance;
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

    float* from_code = malloc(sizeof(float) * nb_patterns);

    for(int i=fmax(from_y - 1, 0); i < fmin(from_y + 1, from_h - 1); i++)
        for(int j=fmax(from_x - 1, 0); j < fmin(from_x + 1, from_w - 1); j++) {

            if(i == from_y && j == from_x)
                continue;

            for(int k=0; k < nb_patterns; k++) {
                from_code[k] = from_codes[k][i][j];
            }

            float distance = distance_modulo_pi(from_code, to_code, nb_patterns);

            // Si la nouvelle distance est plus petite, on update le match
            if(distance == -1.0 || distance < matches[DIST][i][j]) {
                matches[X][i][j]    = to_x;
                matches[Y][i][j]    = to_y;
                matches[DIST][i][i] = distance;
            }
        }

    free(from_code);
}

void hash_from_codes(float*** matches, float*** from_codes, float*** to_codes,
                     int* hash_table[2], int* phases_used, int i, int j,
                     char use_heuristics,
                     int* nb_new_matches, int* nb_better_matches) {

    float* pixel_code = malloc(sizeof(float) * nb_values);
    int k;

    for(k=0; k < nb_values; k++) {
        int phase = phases_used[k];

        pixel_code[k] = from_codes[phase][i][j];
    }

    int hash = hash_fct(pixel_code);

    // Collision = match
    if(hash_table[X][hash] != -1) {
        int x = hash_table[X][hash];
        int y = hash_table[Y][hash];

        float* from_code = malloc(sizeof(float) * nb_patterns);
        float* to_code = malloc(sizeof(float) * nb_patterns);

        for(k=0; k < nb_patterns; k++) {
            from_code[k] = from_codes[k][i][j];
            to_code[k] = to_codes[k][y][x];
        }

        float distance = distance_modulo_pi(from_code, to_code, nb_patterns);

        // Si la nouvelle distance est plus petite, on update le match
        if(matches[DIST][i][j] == -1.0 || distance < matches[DIST][i][j]) {
            if(matches[DIST][i][j] == -1.0)
                (*nb_new_matches)++;
            else
                (*nb_better_matches)++;

            matches[X][i][j] = x;
            matches[Y][i][j] = y;
            matches[DIST][i][j] = distance;

            /* if(use_heuristics) { */
            /*     forward_matching(matches, from_codes, to_codes, */
            /*                      from_code, j, i, x, y); */

            /*     backward_matching(matches, from_codes, to_codes, */
            /*                       to_code, j, i, x, y); */
            /* } */
        }

        free(from_code);
        free(to_code);
    }

    free(pixel_code);
}

void lsh(float*** matches, float*** from_codes, float*** to_codes, char use_heuristics) {

    int i, j;
    //int nb_values = 3; // (inexact->integer (ceiling (/ nb-patterns 2))))
    //int nb_boxes = 10; // (inexact->integer (ceiling (logb (* w h) nbr-values))))
    int nb_collisions = 0, nb_new_matches = 0, nb_better_matches = 0;

    int* hash_table[2];

    for(i=0; i<2; i++) {
        hash_table[i] = malloc(hash_table_size * sizeof(int));

        for(j=0; j<hash_table_size; j++)
            hash_table[i][j] = -1;
    }

    int* phases_used = random_phases(nb_values, nb_patterns);

    int iteration = rand() % 4;

    // Hash des to_codes
    switch(iteration) {
    case 0:
        #pragma omp parallel for private(i, j)
        for(i=0; i < to_h; i++)
            for(j=0; j < to_w; j++)
                HASH_TO_CODES();
        break;

    case 1:
        #pragma omp parallel for private(i, j)
        for(i=0; i < to_h; i++)
            for(j=to_w - 1; j >= 0; j--)
                HASH_TO_CODES();
        break;

    case 2:
        #pragma omp parallel for private(i, j)
        for(i=to_h - 1; i >= 0; i--)
            for(j=0; j < to_w; j++)
                HASH_TO_CODES();
        break;

    case 3:
        #pragma omp parallel for private(i, j)
        for(i=to_h - 1; i >= 0; i--)
            for(j=to_w - 1; j >= 0; j--)
                HASH_TO_CODES();
        break;
    }


    printf("nb_collisions: %d / %d = %f%%\n", nb_collisions, to_w * to_h,
           nb_collisions /(float)(to_w * to_h) * 100);

    printf("hash_space_used: %d / %d = %f%%\n", to_w * to_h - nb_collisions,
           hash_table_size, (to_w * to_h - nb_collisions)/(float)hash_table_size * 100);

    // from_codes hashing (+ updating matches)
    switch(iteration) {
    case 0:
        #pragma omp parallel for private(i, j)
        for(i=0; i < from_h; i++)
            for(j=0; j < from_w; j++) {
                HASH_FROM_CODES();
            }
        break;

    case 1:
        #pragma omp parallel for private(i, j)
        for(i=0; i < from_h; i++)
            for(j=from_w - 1; j >= 0; j--) {
                HASH_FROM_CODES();
            }
        break;

    case 2:
        #pragma omp parallel for private(i, j)
        for(i=from_h - 1; i >= 0; i--)
            for(j=0; j < from_w; j++) {
                HASH_FROM_CODES();
            }
        break;

    case 3:
        #pragma omp parallel for private(i, j)
        for(i=from_h - 1; i >= 0; i--)
            for(j=from_w - 1; j >= 0; j--) {
                HASH_FROM_CODES();
            }
        break;
    }

    free(hash_table[0]);
    free(hash_table[1]);
    free(phases_used);

    printf("nb_new_matches: %d\n", nb_new_matches);
    printf("nb_better_matches: %d\n", nb_better_matches);
}

int main(int argc, char** argv) {

    int i, j, k, l;

    int nthreads = 4, nb_iterations = 30, nb_waves, nb_shifts,
        disable_heuristics = 0, proj_lut = 0, use_quadratic_codes = 0;

    float hash_table_ratio = 2.0;

    char* ref_format = "leo_%d_%d_%03d_%02d.pgm";
    char* cam_format = "%03d.pgm";
    char* ref_phase_format = "phase_ref_%d_%d_%03d.pgm";
    char* cam_phase_format = "phase_cam_%d_%d_%03d.pgm";
    char filename[50];

    // Args parsing
    ARGBEGIN
        ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('p')
        proj_lut = 1;

    ARG_CASE('c')
        cam_format = ARGS;

    ARG_CASE('i')
        nb_iterations = ARGI;

    ARG_CASE('d')
        disable_heuristics = 1;

    ARG_CASE('v')
        nb_values = ARGI;

    ARG_CASE('b')
        nb_boxes = ARGI;

    ARG_CASE('r')
        hash_table_ratio = ARGF;

    ARG_CASE('q')

        use_quadratic_codes = 1;

    WRONG_ARG
        printf("usage: %s [-t nb_threads=%d] [-p gen_proj_lut]\n"
               "\t[-c cam_format=\"%s\"]\n"
               "\t[-i nb_iterations=%d] [-d (disable heuristics)]\n"
               "\t[-b nb_boxes=%d] [-v nb_values=%d]\n"
               "\t[-r hash_table_ratio=%f] [-q use_quadratic_codes]\n",
               argv0, nthreads, cam_format, nb_iterations,
               nb_boxes, nb_values, hash_table_ratio);
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

    fscanf(info, "%d %d %d %d %d",
           &to_w, &to_h, &nb_waves, &nb_patterns, &nb_shifts);

    fclose(info);

    hash_divider = floor((powf(nb_boxes, nb_values) / (to_w * to_h))
                         / hash_table_ratio);
    hash_table_size = (int) ceil(powf(nb_boxes, nb_values) / hash_divider);

    printf("hash_table_size=%d (%f times bigger than %d)\n",
           hash_table_size, hash_table_size /(float)(to_w * to_h),
           to_w * to_h);

    // Lecture d'une image pour trouver le from_w, from_h
    sprintf(filename, cam_format, 0);
    free_f32matrix(load_pgm(filename, &from_w, &from_h));


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

    float*** matches = malloc_f32cube(3, from_w, from_h); // matches[ x, y, distance ][h=i][w=j]

    #pragma omp parallel for private(i, j)
    for(k=0; k<3; k++)
        for(i=0; i<from_h; i++)
            for(j=0; j<from_w; j++)
                matches[k][i][j] = -1.0;

    float*** ref_codes;
    float*** cam_codes;

    cam_codes = load_codes(cam_phase_format, cam_format, !proj_lut, nb_patterns, nb_shifts, from_w, from_h);
    ref_codes = load_codes(ref_phase_format, ref_format, proj_lut, nb_patterns, nb_shifts, to_w, to_h);

    if(use_quadratic_codes) {
        float*** tmp_codes = cam_codes;
        cam_codes = quadratic_codes(cam_codes, nb_patterns, from_w, from_h);
        free_f32cube(tmp_codes, nb_patterns);

        tmp_codes = ref_codes;
        ref_codes = quadratic_codes(ref_codes, nb_patterns, to_w, to_h);
        free_f32cube(tmp_codes, nb_patterns);

        nb_patterns = (nb_patterns + 1) * nb_patterns / 2;
    }

    float* from_code = malloc(sizeof(float) * nb_patterns);
    float* to_code = malloc(sizeof(float) * nb_patterns);

    for(l=0; l<nb_iterations; l++) {
        printf("----- Iteration %02d -----\n", l);
        lsh(matches, cam_codes, ref_codes,
            // heuristics every 5 turn
            !disable_heuristics && l % 5 == 0 && l != 0);

        // Launch heuristics
        if(!disable_heuristics && l % 5 == 0 && l != 0) {
            printf("----- Running backward & forward matching -----\n");

            #pragma omp parallel for private(i, j, k)
            for(i=0; i<from_h; i++)
                for(j=0; j<from_w; j++) {

                    if(matches[DIST][i][j] == -1)
                        continue;

                    for(k=0; k<nb_patterns; k++) {
                        from_code[k] = cam_codes[k][i][j];
                        to_code[k] = ref_codes[k][i][j];
                    }

                    int x = matches[X][i][j];
                    int y = matches[Y][i][j];

                    forward_matching(matches, cam_codes, ref_codes,
                                     from_code, j, i, x, y);

                    backward_matching(matches, cam_codes, ref_codes,
                                      to_code, j, i, x, y);
                }
        }

        // Save the iteration
        if(proj_lut)
            sprintf(filename, "lutProj%02d.ppm", l);
        else
            sprintf(filename, "lutCam%02d.ppm", l);

        save_color_map(filename, matches,
                       from_w, from_h, to_w, to_h, nb_patterns * PI/2.0);
    }

    return EXIT_SUCCESS;
}
