/*
  This program is used with traditional phase-shifting

  Generates 1D (vertical/horizontal) sine images
*/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <math.h>

#include "args.h"
#include "helpers.c"

int format_contains_percentd(char* fmt) {
    int i=0;

    int saw_percent = 0, saw_percent_end = 0, ok = 0;

    while(fmt[i] != '\0') {

        if(!saw_percent && fmt[i] == '%')
            ok = saw_percent = 1;
        else if(saw_percent && !saw_percent_end && fmt[i] == 'd') {
            saw_percent_end = 1;
        } else if(saw_percent && !saw_percent_end && !('0' <= fmt[i] && fmt[i] <= '9')) {
            ok = 0;
        }

        i++;
    }

    return ok && saw_percent_end;
}

int main(int argc, char** argv) {
    int i, j;
    int nthreads = 4, w = 1920, h = 1080, nb_shifts = 1;
    float period = 1000, freq = 1/period, vertical = 0;
    char* output_format = NULL;

    ARGBEGIN
        ARG_CASE('t')
        nthreads = ARGI;

    ARG_CASE('w')
        w = ARGI;

    ARG_CASE('h')
        h = ARGI;

    ARG_CASE('f')
        freq = ARGF;
        period = 1.0/freq;

    ARG_CASE('p')
        period = ARGF;
        freq = 1.0/period;

    ARG_CASE('s')
        nb_shifts = ARGI;

    ARG_CASE('v')
        vertical = 1;

    ARG_CASE('o')
        output_format = ARGS;

    WRONG_ARG
        printf("usage: %s [-t nb_threads=%d] [-w width=%d] [-h height=%d]\n"
               "\t\t[-f frequency=%f] [-p period=%f] [-s nb_shifts=%d] [-v vertical sine]\n"
               "\t\t[-o output_fname]\n"
               "\n"
               "\nif output_format contains %%0*d, it will be used as the shift number\n"
               "\tDefault output_format is :\n"
               "sine_(vertical ? 'v' : 'h')_(w)_(h)_(period)_(nb_shift).pgm\n",
               argv0, nthreads, w, h, freq, period, nb_shifts);
    exit(0);

    ARGEND
    omp_set_num_threads(nthreads);

    // Init image matrix
    float** image = malloc_f32matrix(w, h);

    for(int shift=0; shift < nb_shifts; shift++) {
        memset(image[0], 0, sizeof(float) * h * w);

        // Cos addition
        #pragma omp parallel for private(i, j)
        for(i=0; i < h; i++) {
            for(j=0; j < w; j++) {

                float fx = 0, fy = 0;
                float phase = shift * 2.0 * PI / (float) nb_shifts;

                if(vertical)
                    fy = 2 * PI * freq;
                else
                    fx = 2 * PI * freq;

                image[i][j] = (int) (127.5 + 127.5 * -cosf(fy * i + fx * j + phase + PI));
            }
        }

        char filename[FNAME_MAX_LEN];

        if(output_format == NULL) {
            // Save int periods correctly
            if((int) period == period) {
                sprintf(filename, "sine_%c_%d_%d_%d_%02d.png", (vertical ? 'v' : 'h'), w, h, (int) period, shift);
            } else {
                sprintf(filename, "sine_%c_%d_%d_%f_%02d.png", (vertical ? 'v' : 'h'), w, h, period, shift);
            }

            save_gray_png(filename, image, w, h, 8);

        } else {

            if(format_contains_percentd(output_format)) {

                sprintf(filename, output_format, shift);

            } else {

                sprintf(filename, output_format);

            }

            save_gray_png(filename, image, w, h, 8);
        }
    }

    return EXIT_SUCCESS;
}
