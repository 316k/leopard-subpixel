/*
  Remap grayscale to another random grayscale
  Mapping will probably not be bijective
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../args.h"
#include "../helpers.c"

int main(char argc, char** argv) {

    int nthreads = 4;

    // Args parsing
    ARGBEGIN
    ARG_CASE('s')
        srand(time(NULL));

    WRONG_ARG
        usage:
        printf("usage: %s [-s random seed] gray.png out.png\n", argv0);
        exit(1);
    ARGEND

    if(argc != 2) goto usage;

    int size, w, h;

    FILE *f = fopen(argv[0], "r");

    require_file(f, argv[0]);

    read_image_header(f, &w, &h, &size);
    fclose(f);

    int depth = size == 65535 ? 16 : 8;

    float** in = load_gray(argv[0], &w, &h);

    int nbvalues = (int) powl(2, depth);

    int* map = malloc(sizeof(int) * nbvalues);

    for(int v=0; v<nbvalues; v++) {
        map[v] = rand() % nbvalues;
    }

    #pragma omp parallel for
    for(int i=0; i < h; i++) {
        for(int j=0; j < w; j++) {
            int val = (int) in[i][j];
            int remapped = (int) map[val];

            in[i][j] = remapped;
        }
    }

    save_gray_png(argv[1], in, w, h, depth);

    return 0;
}
