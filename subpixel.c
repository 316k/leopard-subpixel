#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helpers.c"

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, l, from_w, from_h, to_w, to_h, foo, nb_shifts, nb_patterns;
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
                matches[X][i][j] = round(matches[X][i][j] / 65535.0 * to_w);
                matches[Y][i][j] = round(matches[Y][i][j] / 65535.0 * to_h);
                matches[DIST][i][j] = matches[DIST][i][j] / 65535.0 * (nb_patterns * PI / 2.0);
            }
        }

    
    float*** subpixel = malloc_f32cube(3, from_w, from_h);

    float*** cam_codes = load_codes(cam_phase_format, cam_format, 1, nb_patterns, nb_shifts, from_w, from_h);
    float*** ref_codes = load_codes(ref_phase_format, ref_format, 0, nb_patterns, nb_shifts, to_w, to_h);

    float*** supixeled = malloc_f32cube(2, from_w, from_h);

    // #pragma omp parallel for private(i, j, k, l)
    for(i=0; i<from_h; i++)
        for(j=0; j<from_w; j++) {
            
            int nx = 0;
            int ny = 0;
            
            int x = matches[X][i][j];
            int y = matches[Y][i][j];

            // Undefined matches stay undefined
            if(x == -1) {
                subpixel[X][i][j] = subpixel[Y][i][j] = subpixel[DIST][i][j] = -1;
                continue;
            }

            float val, match, original, neighbour;
            
            for(k=0; k<nb_patterns; k++) {
                
                match = cam_codes[k][i][j];
                original = ref_codes[k][i][j];

                // Left
                if(x != 0 && ((neighbour = ref_codes[k][y][x - 1]) || 1) &&
                   ((neighbour < match && match < original) ||
                    (neighbour > match && match > original))) { // borders are problematic

                    // printf("l %f %f %f %f\n", original, match, neighbour, 0.5 - ((match - original) / (neighbour - original))/2.0);
                    subpixel[X][i][j] += 0.5 - ((match - original) / (neighbour - original));
                    nx++;
                    
                } else if(x != to_w - 1 && ((neighbour = ref_codes[k][y][x + 1]) || 1) && // Right
                          ((neighbour < match && match < original) ||
                           (neighbour > match && match > original))) {
                    
                    // printf("r %f %f %f %f\n", original, match, neighbour, 0.5 + ((match - original) / (neighbour - original))/2.0);
                    subpixel[X][i][j] += 0.5 + ((match - original) / (neighbour - original));
                    nx++;
                }

                // Down
                if(y != to_h - 1 && ((neighbour = ref_codes[k][y + 1][x]) || 1) &&
                   ((neighbour < match && match < original) ||
                    (neighbour > match && match > original))) {

                    // printf("d %f %f %f %f\n", original, match, neighbour, 0.5 + ((match - original) / (neighbour - original))/2.0);
                    subpixel[Y][i][j] += 0.5 + ((match - original) / (neighbour - original));
                    ny++;
                    
                } else if(y != 0 && ((neighbour = ref_codes[k][y - 1][x]) || 1) && // Up
                          ((neighbour < match && match < original) ||
                           (neighbour > match && match > original))) {
                    
                    //printf("u %f %f %f %f\n", original, match, neighbour, 0.5 - ((match - original) / (neighbour - original))/2.0);
                    subpixel[Y][i][j] += 0.5 - ((match - original) / (neighbour - original));
                    ny++;
                }
                
                // TODO diago                
            }
            
            printf("%d %d => %d %d\n", i, j, y, x);
            // Use averages
            if(nx)
                subpixel[X][i][j] /= nx;
            
            if(ny)
                subpixel[Y][i][j] /= ny;
            
            // Keep distance information
            subpixel[DIST][i][j] = matches[DIST][i][j];
            
            for(k=0; k<2; k++)
                subpixel[k][i][j] += matches[k][i][j];
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
