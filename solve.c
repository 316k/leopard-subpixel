#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "helpers.c"

#define X 0
#define Y 1
#define DIST 2

int hash_fct(float* pixel_code, int nb_values, int nb_boxes) {
    int hash = 0;

    for(int i=0; i<nb_values; i++) {
        hash *= nb_boxes;
        hash += fmod(pixel_code[i] + PI, 2 * PI) / (2 * PI) * nb_boxes;
    }

    return hash;
}

void lsh(float*** matches, float*** from_codes, float*** to_codes, int nb_patterns,
         int from_w, int from_h, int to_w, int to_h) {

    int i, j, k, l;
    int nb_values = 8; // (inexact->integer (ceiling (/ nb-patterns 2))))
    int nb_boxes = 9; // (inexact->integer (ceiling (logb (* w h) nbr-values))))
    int nb_collisions = 0, nb_new_matches = 0, nb_better_matches = 0;
    
    int* hash_table[2];
    int hash_table_size = powf(nb_boxes, nb_values);

    printf("hash_table_size %d is %f times bigger than ~%d\n", hash_table_size, hash_table_size /(float)(to_w * to_h), to_w * to_h);

    for(i=0; i<2; i++) {
        hash_table[i] = malloc(hash_table_size * sizeof(int));
        
        for(j=0; j<hash_table_size; j++)
            hash_table[i][j] = -1;
    }

    int* phases_used = random_phases(nb_values, nb_patterns);

    // Random iteration
    int to_i_start, to_i_end, to_j_start, to_j_end,
        from_i_start, from_i_end, from_j_start,
        from_j_end, i_increment, j_increment;

    int iteration = rand() % 4;

    // i flipped ?
    if(iteration % 2) {
        to_i_start   = 0;
        to_i_end     = to_h;
        from_i_start = 0;
        from_i_end   = from_h;
        i_increment  = 1;
    } else {
        to_i_start   = to_h - 1;
        to_i_end     = -1;
        from_i_start = from_h - 1;
        from_i_end   = -1;
        i_increment  = -1;
    }

    // j flipped ?
    if(iteration >= 2) {
        to_j_start   = 0;
        to_j_end     = to_w;
        from_j_start = 0;
        from_j_end   = from_w;
        j_increment  = 1;
    } else {
        to_j_start   = to_w - 1;
        to_j_end     = -1;
        from_j_start = from_w - 1;
        from_j_end   = -1;
        j_increment  = -1;
    }
    
    // Hash des to_codes
    for(i=to_i_start; i != to_i_end; i += i_increment)
        for(j=to_j_start; j != to_j_end; j += j_increment) {
            float* pixel_code = malloc(sizeof(float) * nb_patterns);

            for(k=0; k < nb_values; k++) {
                int phase = phases_used[k];
                
                assert(phase < nb_patterns);
                
                pixel_code[k] = to_codes[phase][i][j];
            }
            
            int hash = hash_fct(pixel_code, nb_values, nb_boxes);
            
            // Premier arrivé prend la place
            if(hash_table[X][hash] == -1) {
                hash_table[X][hash] = j;
                hash_table[Y][hash] = i;
            } else {
                nb_collisions++;
            }
            
            free(pixel_code);
        }

    printf("nb_collisions: %d / %d = %f%%\n", nb_collisions, to_w * to_h, nb_collisions /(float)(to_w * to_h) * 100);
    printf("hash_space_used: %d / %d = %f%%\n", to_w * to_h - nb_collisions,
           hash_table_size, (to_w * to_h - nb_collisions)/(float)hash_table_size * 100);

    // Hash des to_codes
    for(i=from_i_start; i != from_i_end; i += i_increment)
        for(j=from_j_start; j != from_j_end; j += j_increment) {
            float* pixel_code = malloc(sizeof(float) * nb_patterns);
            
            for(k=0; k < nb_values; k++) {
                int phase = phases_used[k];
                
                pixel_code[k] = from_codes[phase][i][j];
            }
            
            int hash = hash_fct(pixel_code, nb_values, nb_boxes);
            
            // Collision = match
            if(hash_table[X][hash] != -1) {
                int x = hash_table[X][hash];
                int y = hash_table[Y][hash];

                float* to_code = malloc(sizeof(float) * nb_patterns);

                for(k=0; k < nb_values; k++) {
                    int phase = phases_used[k];
                    
                    to_code[k] = to_codes[phase][y][x];
                }
                
                float distance = distance_modulo_pi(pixel_code, to_code, nb_values);
                
                // Si la nouvelle distance est plus petite, on update le match
                if(matches[DIST][i][j] == -1.0 || distance < matches[DIST][i][j]) {
                    if(matches[DIST][i][j] == -1.0)
                        nb_new_matches++;
                    else
                        nb_better_matches++;
                    
                    matches[X][i][j] = x;
                    matches[Y][i][j] = y;
                    matches[DIST][i][j] = distance;
                }
                
                free(to_code);
            }
            
            free(pixel_code);
        }

    free(hash_table[0]);
    free(hash_table[1]);
    free(phases_used);
    
    printf("nb_new_matches: %d\n", nb_new_matches);
    printf("nb_better_matches: %d\n", nb_better_matches);
}

void save_color_map(char* filename, float*** matches, int from_w, int from_h,
                    int to_w, int to_h, float max_distance) {
    float*** gray_levels = malloc_f32cube(3, from_w, from_h);
    
    // Color map
    for(int i=0; i<from_h; i++)
        for(int j=0; j<from_w; j++) {
            if(matches[X][i][j] == -1.0) {
                gray_levels[X][i][j] = 0.0;
                gray_levels[Y][i][j] = 0.0;
                gray_levels[DIST][i][j] = 65535.0;
            } else {
                gray_levels[X][i][j] = matches[X][i][j] * 65535.0/(float)to_w;
                gray_levels[Y][i][j] = matches[Y][i][j] * 65535.0/(float)to_h;
                gray_levels[DIST][i][j] = fmin(matches[DIST][i][j] * 65535.0/(float)(max_distance), 65535.0);
            }
        }
    
    save_ppm16(filename, gray_levels[X], gray_levels[Y], gray_levels[DIST], from_w, from_h);
    
    free_f32cube(gray_levels, 3);
}

int main(char argc, char** argv) {

    int i, j, k, foo, shift;
    FILE* info = fopen("sines.txt", "r");

    int nthreads = 4, nb_iterations = 30, from_w, from_h, to_w, to_h, nb_waves, nb_shifts, nb_patterns;

    char* ref_format = "leo_%d_%d_%03d_%02d.pgm";
    char* cam_format = "%03d.pgm";
    char* ref_phase_format = "phase_ref_%d_%d_%03d.pgm";
    char* cam_phase_format = "phase_cam_%d_%d_%03d.pgm";
    char filename[50];
    
    // Args parsing
    for(i=1; i < argc; i++) {
        if(strcmp(argv[i], "-t") == 0) {
            nthreads = atoi(argv[i + 1]); i++;
        } else if(strcmp(argv[i], "-c") == 0) {
            cam_format = argv[i + 1]; i++;
        } else if(strcmp(argv[i], "-i") == 0) {
            nb_iterations = atoi(argv[i + 1]); i++;
        } else {
            printf("usage: %s [-t nb_threads=%d] [-c cam_format=\"%s\"] [-i nb_iterations=%d]\n",
                   argv[0], nthreads, cam_format, nb_iterations);
            exit(0);
        }
    }
    
    omp_set_num_threads(nthreads);

    srand(time(NULL));
    
    fscanf(info, "%d %d %d %d %d", &to_w, &to_h, &nb_waves, &nb_patterns, &nb_shifts);

    fclose(info);

    // Lecture d'une image pour trouver le from_w, from_h
    sprintf(filename, cam_format, 0);
    free_f32matrix(load_pgm(filename, &from_w, &from_h));
    
    float*** matches = malloc_f32cube(3, from_w, from_h); // matches[ x, y, distance ][h=i][w=j]

    #pragma omp parallel for private(i, j)
    for(k=0; k<3; k++)
        for(i=0; i<from_h; i++)
            for(j=0; j<from_w; j++)
                matches[k][i][j] = -1.0; 

    float*** ref_codes;
    float*** cam_codes;

    
    cam_codes = load_codes(cam_phase_format, cam_format, 1, nb_patterns, nb_shifts, from_w, from_h);
    
    ref_codes = load_codes(ref_phase_format, ref_format, 0, nb_patterns, nb_shifts, to_w, to_h);
    
    for(i=0; i<nb_iterations; i++) {
        printf("----- Iteration %02d -----\n", i);
        lsh(matches, cam_codes, ref_codes, nb_patterns, from_w, from_h, to_w, to_h);

        sprintf(filename, "matches-%02d.pgm", i);
        save_color_map(filename, matches, from_w, from_h, to_w, to_h, nb_patterns * PI/2.0);
    }

    return EXIT_SUCCESS;
}
