#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helpers.c"

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
float subpixel_value(float u, float v, float a, float b, float c, float d) {
    
    return (1 - v) * ((1 - u) * a + u * b) + v * ((1 - u) * c + u * d);
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

    int nthreads = 4, i, j, k, from_w, from_h, to_w, to_h, foo, nb_shifts, nb_patterns;
    char* ref_format = "leo_%d_%d_%03d_%02d.pgm";
    char* cam_format = "%03d.pgm";
    
    FILE* info = fopen("sines.txt", "r");

    // Check file size to avoid problems if sines.txt is empty
    fseek(info, 0, SEEK_END);
    if(!ftell(info)) {
        printf("error: empty sines.txt\n");
        exit(-1);
    }
    fseek(info, 0, SEEK_SET);

    // Args parsing
    for(i=1; i < argc - 1; i++) {
        if(strcmp(argv[i], "-t") == 0) {
            nthreads = atoi(argv[i + 1]); i++;
        } else {
        usage:
            printf("usage: %s [-t nb_threads=%d] filename\n",
                   argv[0], nthreads);
            exit(1);
        }
    }
    if(i != argc - 1) goto usage;

    omp_set_num_threads(nthreads);
    
    srand(time(NULL));
    
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
                matches[X][i][j] = floor(matches[X][i][j] / 65535.0 * to_w);
                matches[Y][i][j] = floor(matches[Y][i][j] / 65535.0 * to_h);
                matches[DIST][i][j] = matches[DIST][i][j] / 65535.0 * (nb_patterns * PI / 2.0);
            }
        }
    
    float*** subpixel = malloc_f32cube(3, from_w, from_h);

    float*** cam_codes = load_codes(cam_phase_format, cam_format, 1, nb_patterns, nb_shifts, from_w, from_h);
    float*** ref_codes = load_codes(ref_phase_format, ref_format, 0, nb_patterns, nb_shifts, to_w, to_h);
    
    #pragma omp parallel for private(i, j, k)
    for(i=0; i<from_h; i++) {
        printf("%d\n", i);
        for(j=0; j<from_w; j++) {
            
            int x = matches[X][i][j];
            int y = matches[Y][i][j];

            // Undefined matches stay undefined
            if(x < 0) {
                subpixel[X][i][j] = subpixel[Y][i][j] = subpixel[DIST][i][j] = -1;
                continue;
            }
            
            float* match = malloc(sizeof(float) * nb_patterns);

            for(k=0; k<nb_patterns; k++) {
                match[k] = cam_codes[k][i][j];
            }
            
            float** costs = malloc_f32matrix(10, 10);

            float quadrant_best[4];
            float decalage_x[4], decalage_y[4];
            
            for(k=0; k<4; k++)
                quadrant_best[k] = INFINITY;
            
            // Minimize subpixel value
            
            // Up-right, Up-left, Down-right, down-left
            int pos_x[] = {+1, -1, +1, -1};
            int pos_y[] = {-1, -1, +1, +1};

            int q; // current quadrant
            
            for(q=0; q<4; q++) {
                
                // Assure que les bornes ne sont pas dépassées
                if((pos_x[q] == -1 && x == 0) ||
                   (pos_y[q] == -1 && y == 0) ||
                   (pos_x[q] == +1 && x == to_w - 1) ||
                   (pos_y[q] == +1 && y == to_h - 1))
                    continue;
                
                for(int u=0; u<10; u++) {
                    for(int v=0; v<10; v++) {

                        int nb_costs = 0;
                        
                        for(k=0; k<nb_patterns; k++) {
                            
                            float m = match[k],
                                a = ref_codes[k][y][x],
                                b = ref_codes[k][y][x + pos_x[q]],
                                c = ref_codes[k][y + pos_y[q]][x],
                                d = ref_codes[k][y + pos_y[q]][x + pos_x[q]];
                            
                            // la phase matchée doit être un sous-pixel possible dans le quadrant
                            // pour que le coût ait un sens
                            // TODO : Considérer le warp dans [-PI,PI[
                            if((m > a && m > b && m > c && m > d) ||
                               (m < a && m < b && m < c && m < d)) {
                                costs[v][u] += PI;
                                continue;
                            }
                            
                            costs[v][u] += fabsl(match[k] - subpixel_value(u / (float)10, v / (float)10,
                                                                           a,b,c,d));
                            nb_costs++;
                        }
                    }
                }
                
                quadrant_best[q] = f32matrix_min(costs, &decalage_x[q], &decalage_y[q], 10, 10);
                decalage_x[q] = pos_x[q] * decalage_x[q] / 2.0 + 5;
                decalage_y[q] = pos_y[q] * decalage_y[q] / 2.0 + 5;
                    
                for(int u=0; u<10; u++)
                    for(int v=0; v<10; v++)
                        costs[u][v] = 0.0;
            }
                        
            free(match);
            free_f32matrix(costs);
            
            // Trouve le meilleur match de sous-pixel
            float min = INFINITY;
            int index = -1;
            for(k=0; k<4; k++) {
                // printf("%d: %f (%f %f)\n", k, quadrant_best[k], decalage_x[k], decalage_y[k]);
                if(quadrant_best[k] < min) {
                    min = quadrant_best[k];
                    index = k;
                }
            }
            // printf("Choice : %d\n", index);

            if(index != -1) {
                subpixel[X][i][j] = decalage_x[index] / 10.0;
                subpixel[Y][i][j] = decalage_y[index] / 10.0;
            } else {
                subpixel[X][i][j] = 0.5;
                subpixel[Y][i][j] = 0.5;
            }
            
            // Keep distance information
            subpixel[DIST][i][j] = matches[DIST][i][j];
            
            for(k=0; k<2; k++)
                subpixel[k][i][j] += matches[k][i][j];
        }
    }

    FILE* vals;
    
    for(k=0; k<2; k++) {
        float min = 100000;
        float dec = 0;
        float max = -100000;

        if(k == 0)
            vals = fopen("subpixel-x-vals", "w");
        else
            vals = fopen("subpixel-y-vals", "w");
        
        for(i=1; i<from_h-1; i++)
            for(j=1; j<from_w-1; j++) {
                dec += subpixel[k][i][j] - matches[k][i][j];
                min = fmin(min, subpixel[k][i][j] - matches[k][i][j]);
                max = fmax(max, subpixel[k][i][j] - matches[k][i][j]);
                fprintf(vals, "%f\n", subpixel[k][i][j] - matches[k][i][j]);
            }
        
        dec /= (from_h - 2) * (from_w - 2);
        fprintf(stderr, "avg=%f min=%f max=%f\n", dec, min, max);
        fclose(vals);
    }
    
    save_color_map("subpixel.ppm", subpixel, from_w, from_h, to_w, to_h, nb_patterns * PI/2.0);
    
    return 0;
}