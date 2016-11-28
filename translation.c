#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helpers.c"

float get(float** matrix, int i, int j, int w, int h) {

    if(i < 0 || i >= h || j < 0 || j >= w) {
        return 0;
    }

    return matrix[i][j];
    
}

int main(char argc, char** argv) {

    int nthreads = 4, i, j, k, foo, shift, w, h;

    float decalage = 0.0;

    // Args parsing
    for(i=1; i < argc - 1; i++) {
        if(strcmp(argv[i], "-t") == 0) {
            nthreads = atoi(argv[i + 1]); i++;
        } else if(strcmp(argv[i], "-d") == 0) {
            decalage = atof(argv[i + 1]); i++;
            if(!(decalage <= 1 && decalage >= 0)) goto usage;
        } else {
        usage:
            printf("usage: %s [-t nb_threads=%d] [-d decalage=%f] filename\n",
                   argv[0], nthreads, decalage);
            exit(1);
        }
    }

    if(i != argc - 1) goto usage;
    
    float** in = load_pgm(argv[argc - 1], &w, &h);

    float** out = malloc_f32matrix(w, h);

    for(i=0; i < h; i++)
        for(j=0; j<w; j++) {
            out[i][j] = (get(in, i + 1, j, w, h) * decalage + (1 - decalage) * get(in, i, j, w, h));
        }
    
    save_pgm(NULL, out, w, h, 8);
    
    return 0;
}
