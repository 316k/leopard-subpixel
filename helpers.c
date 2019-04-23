#include <unistd.h>
#include <errno.h>
#include "lodepng.c"

extern int errno;

#define X 0
#define Y 1
#define DIST 2

#define FNAME_MAX_LEN 50

const float PI = 2 * atan2(1, 0);

void require_file(FILE* f, char* fname) {
    if(f == NULL) {
        fprintf(stderr, "%s: %s\n", fname, strerror(errno));
        exit(errno);
    }
}

// --------------- Images ---------------
float** malloc_f32matrix(int w, int h) {
    float** image = malloc(sizeof(float*) * h);
    image[0] = (float*) calloc(h * w, sizeof(float));

    for(int i=0; i < h; i++) {
        image[i] = (*image + w * i);
    }

    return image;
}

void free_f32matrix(float** image) {
    free(image[0]);
    free(image);
}

float*** malloc_f32cube(int w, int h, int k) {
    float*** cube = malloc(sizeof(float**) * w);

    for(int i=0; i<w; i++) {
        cube[i] = malloc_f32matrix(h, k);
    }

    return cube;
}

void free_f32cube(float*** cube, int k) {

    for(int i=0; i < k; i++) {
        free_f32matrix(cube[i]);
    }

    free(cube);
}

// TODO : use this function for more portability
char** image_fnames(char* img_format, char numbered_imgs, int nb_patterns, int nb_shifts, int w, int h) {

    int nb_imgs = nb_patterns * nb_shifts;

    char** fnames = malloc(sizeof(char*) * nb_imgs);
    fnames[0] = (char*) calloc(nb_imgs * FNAME_MAX_LEN, sizeof(char));

    for(int i=0; i < nb_patterns; i++) {
        for(int shift=0; shift < nb_shifts; shift++) {

            int img_nb = i * nb_shifts + shift;

            fnames[img_nb] = (*fnames + FNAME_MAX_LEN * img_nb);

            if(numbered_imgs)
                sprintf(fnames[img_nb], img_format, img_nb);
            else
                sprintf(fnames[img_nb], img_format, w, h, i, shift);
        }
    }

    return fnames;
}

void read_image_header(FILE* f, int *w, int *h, int *size) {
    char foo[100];

    // Magic number (P5, P6)
    fgets(foo, 100, f);

    // "w h" or "#" (comment)
    fgets(foo, 100, f);

    while(foo[0] == '#')
        fgets(foo, 100, f);

    sscanf(foo, "%d %d", w, h);

    fgets(foo, 100, f);

    sscanf(foo, "%d", size);
}

float** load_pgm(char* name, int *w, int *h) {
    int i, j, v, size;
    float** mat;
    FILE *f = fopen(name, "r");

    require_file(f, name);

    read_image_header(f, w, h, &size);

    size = size == 65535 ? 2 : 1;

    // Read
    mat = malloc_f32matrix(*w, *h);

    for(i=0; i < *h; i++)
        for(j=0; j < *w; j++) {
            if(size == 2) {
                // Higher
                v = getc(f) << 8;
                // Lower
                v += getc(f);
            } else {
                v = getc(f);
            }

            mat[i][j] = v;
        }

    fclose(f);

    return mat;
}

float*** load_ppm(char* name, int *w, int *h) {

    int i, j, k, v, size;
    float*** mat = malloc(sizeof(float**) * 3); // R,G,B channels

    FILE *f = fopen(name, "r");

    require_file(f, name);

    read_image_header(f, w, h, &size);

    size = size == 65535 ? 2 : 1;

    // Read
    for(i=0; i<3; i++)
        mat[i] = malloc_f32matrix(*w, *h);

    for(i=0; i < *h; i++)
        for(j=0; j < *w; j++) {
            for(k=0; k<3; k++) {
                if(size == 2) {
                    // Higher
                    v = getc(f) << 8;
                    // Lower
                    v += getc(f);
                } else {
                    v = getc(f);
                }

                mat[k][i][j] = v;
            }
        }

    fclose(f);

    return mat;
}

void save_pgm(char* filename, float** image, int w, int h, int depth) {
    int i, j, v;

    int max = depth == 16 ? 65535 : 255;

    FILE* out;
    if(filename != NULL)
        out = fopen(filename, "w+");
    else
        out = stdout;

    fprintf(out, "P5\n#\n%d %d\n%d\n", w, h, max);
    for(i=0; i < h; i++) {
        for(j=0; j<w; j++) {

            v = image[i][j];

            if(v > max || v < 0) {
                fprintf(stderr, "*** WARNING : (%d, %d) = %f\n", i, j, image[i][j]);
                v = fmin(fmax(v, 0), max);
            }

            if(depth == 16) {
                // Higher
                fputc((v >> 8), out);
                // Lower
                fputc((v & 255), out);
            } else {
                fputc(v, out);
            }
        }
    }

    fclose(out);
}

void save_ppm(char* filename, float** channels[3], int w, int h, int depth) {
    int i, j, k, v;
    int max = depth == 16 ? 65535 : 255;

    FILE* out;

    if(filename != NULL)
        out = fopen(filename, "w+");
    else
        out = stdout;

    fprintf(out, "P6\n#\n%d %d\n%d\n", w, h, max);
    for(i=0; i < h; i++) {
        for(j=0; j<w; j++) {

            for(k=0; k<3; k++) {
                v = channels[k][i][j];

                if(v > max || v < 0) {
                    fprintf(stderr, "*** WARNING : (%d, %d) = %f [channel=%d]\n", i, j, channels[k][i][j], k);
                    v = fmin(fmax(v, 0), max);
                }

                if(depth == 16) {
                    // Higher
                    fputc((v >> 8), out);
                    // Lower
                    fputc((v & 255), out);
                } else {
                    fputc(v, out);
                }
            }
        }
    }

    fclose(out);
}

unsigned char* load_png(char* fname, int *w, int *h, int *depth, LodePNGColorType expected_type) {

    unsigned int error;
    unsigned char* image;
    unsigned int width, height;
    unsigned char* png = 0;
    size_t pngsize;

    LodePNGState state;

    lodepng_state_init(&state);

    state.info_raw.colortype = expected_type;

    error = lodepng_load_file(&png, &pngsize, fname);

    lodepng_inspect(&width, &height, &state, png, pngsize);

    *depth = (int) state.info_png.color.bitdepth;
    state.info_raw.bitdepth = state.info_png.color.bitdepth;

    if(!error)
        error = lodepng_decode(&image, &width, &height, &state, png, pngsize);

    if(state.info_png.color.colortype != expected_type) {
        fprintf(stderr, "%s: color type is not %s\n", fname, expected_type == LCT_RGB ? "RGB" : "Grayscale");
        exit(-1);
    }

    if(error) {
        fprintf(stderr, "error %u: %s\n", error, lodepng_error_text(error));
        exit(error);
    }

    *w = (int) width;
    *h = (int) height;

    free(png);
    lodepng_state_cleanup(&state);

    return image;
}

float** load_png_gray(char* fname, int *w, int *h) {

    int depth;

    unsigned char* image = load_png(fname, w, h, &depth, LCT_GREY);

    float** array = malloc_f32matrix(*w, *h);

    unsigned int n = 0;
    for(int i=0; i<*h; i++)
        for(int j=0; j<*w; j++) {
            int val = image[n++];

            if(depth == 16)
                val = (val << 8) + image[n++];

            array[i][j] = val;
        }

    free(image);

    return array;
}

float*** load_png_rgb(char* fname, int *w, int *h) {

    int depth;

    unsigned char* image = load_png(fname, w, h, &depth, LCT_RGB);

    float*** array = malloc_f32cube(3, *w, *h);
    unsigned int n = 0;
    for(int i=0; i<*h; i++)
        for(int j=0; j<*w; j++) {
            for(int k=0; k<3; k++) {
                int val = image[n++];

                if(depth == 16)
                    val = (val << 8) + image[n++];

                array[k][i][j] = val;
            }
        }

    free(image);

    return array;
}

char is_png(char* fname) {

    // Ensures file exists
    FILE* f = fopen(fname, "r");
    require_file(f, fname);
    fclose(f);

    int end = strlen(fname) - 4;

    if(end < 0)
        end = 0;

    return strcmp(fname + end, ".png") == 0;
}

float*** load_color(char* fname, int *w, int *h) {

    if(is_png(fname))
        return load_png_rgb(fname, w, h);

    return load_ppm(fname, w, h);
}

float** load_gray(char* fname, int *w, int *h) {

    if(is_png(fname))
        return load_png_gray(fname, w, h);

    return load_pgm(fname, w, h);
}

void save_raw_png(char* filename, unsigned char* image, int w, int h, int depth, LodePNGColorType colortype) {

    unsigned int width = w;
    unsigned int height = h;

    unsigned int error;
    unsigned char* png;
    size_t pngsize;
    LodePNGState state;

    lodepng_state_init(&state);

    state.encoder.auto_convert = 0;

    state.info_raw.colortype = colortype;
    state.info_raw.bitdepth = depth;

    state.info_png.color.colortype = colortype;
    state.info_png.color.bitdepth = depth;

    error = lodepng_encode(&png, &pngsize, image, width, height, &state);

    if(!error)
        error = lodepng_save_file(png, pngsize, filename);

    if(error) {
        fprintf(stderr, "error %u: %s\n", error, lodepng_error_text(error));
    }

    lodepng_state_cleanup(&state);
    free(png);
}

void save_gray_png(char* filename, float** image, int w, int h, int depth) {

    unsigned char* array = malloc(w * h * (depth == 16 ? 2 : 1));

    unsigned int n = 0;

    int max = depth == 16 ? 65535 : 255;

    for(int i=0; i<h; i++)
        for(int j=0; j<w; j++) {
                int val = image[i][j];

                if(val > max || val < 0) {
                    fprintf(stderr, "*** WARNING : (%d, %d) = %f\n", i, j, image[i][j]);
                    val = fmin(fmax(val, 0), max);
                }

                if(depth == 16) {
                    array[n++] = val >> 8;
                }

                array[n++] = val & 255;
        }

    save_raw_png(filename, array, w, h, depth, LCT_GREY);

    free(array);
}

void save_color_png(char* filename, float*** image, int w, int h, int depth) {

    unsigned char* array = malloc(w * h * 3 * (depth == 16 ? 2 : 1));

    unsigned int n = 0;

    int max = depth == 16 ? 65535 : 255;

    for(int i=0; i<h; i++)
        for(int j=0; j<w; j++) {
            for(int k=0; k<3; k++) {
                int val = image[k][i][j];

                if(val > max || val < 0) {
                    fprintf(stderr, "*** WARNING : (%d, %d) = %f [channel=%d]\n", i, j, image[k][i][j], k);
                    val = fmin(fmax(val, 0), max);
                }

                if(depth == 16) {
                    array[n++] = val >> 8;
                }

                array[n++] = val & 255;
            }
        }

    save_raw_png(filename, array, w, h, depth, LCT_RGB);

    free(array);
}

// --------------- Misc calculations ---------------

int* random_phases(int size, int total) {

    int i, j, swap;
    int* list = malloc(sizeof(int) * total);

    for(i=0; i<total; i++)
        list[i] = i;


    // Shuffle
    for(i=0; i<total; i++) {
        j = rand() % total;
        swap = list[i];
        list[i] = list[j];
        list[j] = swap;
    }

    int* out = malloc(sizeof(int) * size);

    for(i=0; i<size; i++) {
        out[i] = list[i];
    }

    free(list);

    return out;
}

float solve_phase_term(float* intensities, int nb_shifts, float (*trigo)(float)) {
    float total = 0;

    for(int i=0; i<nb_shifts; i++) {
        total += trigo(2 * i * PI /(float) nb_shifts) * intensities[i];
    }

    return total;
}

float solve_phase(float* intensities, int nb_shifts) {
    return atan2(solve_phase_term(intensities, nb_shifts, sinf),
                 solve_phase_term(intensities, nb_shifts, cosf));
}

float modulo_sub_pi(float x, float y) {
    float diff = fabs(fmod(PI + x, PI) - fmod(PI + y, PI));

    return fminf(diff, PI - diff);
}

float distance_modulo_pi(float* a, float* b, int len) {
    float total = 0;

    for(int i=0; i<len; i++)
        total += modulo_sub_pi(a[i], b[i]);

    return total;
}

// Mathematica : Table[ ToString[x] -> Round[Erfc[x]/2*256], {x, -2, 2, 0.01}]
float gray_scale_erfc(float value) {
    int n = roundf(value * 100);

    if(n <= -179)
        return 255;

    // Don't worry, it's all generated...
    switch(n) {
    //-1.78
    case -178:
    case -177:
    case -176:
    case -175:
    case -174:
    case -173:
    case -172:
    case -171:
    case -170:
    case -169:
    case -168:
    case -167:
    case -166: return 254;
    case -165:
    case -164:
    case -163:
    case -162:
    case -161:
    case -160:
    case -159:
    case -158:
    case -157: return 253;
    case -156:
    case -155:
    case -154:
    case -153:
    case -152:
    case -151:
    case -150:
    case -149: return 252;
    case -148:
    case -147:
    case -146:
    case -145:
    case -144: return 251;
    case -143:
    case -142:
    case -141:
    case -140:
    case -139: return 250;
    case -138:
    case -137:
    case -136:
    case -135:
    case -134: return 249;
    case -133:
    case -132:
    case -131:
    case -130: return 248;
    case -129:
    case -128:
    case -127: return 247;
    case -126:
    case -125:
    case -124:
    case -123: return 246;
    case -122: return 245;
    case -121: return 245;
    case -120: return 245;
    case -119: return 244;
    case -118: return 244;
    case -117: return 243;
    case -116: return 243;
    case -115: return 243;
    case -114: return 242;
    case -113: return 242;
    case -112: return 242;
    case -111: return 241;
    case -110: return 241;
    case -109: return 240;
    case -108: return 240;
    case -107: return 239;
    case -106: return 239;
    case -105: return 238;
    case -104: return 238;
    case -103: return 237;
    case -102: return 237;
    case -101: return 236;
    case -100: return 236;
    case -99: return 235;
    case -98: return 235;
    case -97: return 234;
    case -96: return 234;
    case -95: return 233;
    case -94: return 232;
    case -93: return 232;
    case -92: return 231;
    case -91: return 231;
    case -90: return 230;
    case -89: return 229;
    case -88: return 229;
    case -87: return 228;
    case -86: return 227;
    case -85: return 227;
    case -84: return 226;
    case -83: return 225;
    case -82: return 224;
    case -81: return 224;
    case -80: return 223;
    case -79: return 222;
    case -78: return 221;
    case -77: return 221;
    case -76: return 220;
    case -75: return 219;
    case -74: return 218;
    case -73: return 217;
    case -72: return 217;
    case -71: return 216;
    case -70: return 215;
    case -69: return 214;
    case -68: return 213;
    case -67: return 212;
    case -66: return 211;
    case -65: return 210;
    case -64: return 209;
    case -63: return 208;
    case -62: return 207;
    case -61: return 206;
    case -60: return 205;
    case -59: return 204;
    case -58: return 203;
    case -57: return 202;
    case -56: return 201;
    case -55: return 200;
    case -54: return 199;
    case -53: return 198;
    case -52: return 197;
    case -51: return 196;
    case -50: return 195;
    case -49: return 193;
    case -48: return 192;
    case -47: return 191;
    case -46: return 190;
    case -45: return 189;
    case -44: return 188;
    case -43: return 186;
    case -42: return 185;
    case -41: return 184;
    case -40: return 183;
    case -39: return 182;
    case -38: return 180;
    case -37: return 179;
    case -36: return 178;
    case -35: return 177;
    case -34: return 175;
    case -33: return 174;
    case -32: return 173;
    case -31: return 171;
    case -30: return 170;
    case -29: return 169;
    case -28: return 167;
    case -27: return 166;
    case -26: return 165;
    case -25: return 163;
    case -24: return 162;
    case -23: return 161;
    case -22: return 159;
    case -21: return 158;
    case -20: return 157;
    case -19: return 155;
    case -18: return 154;
    case -17: return 152;
    case -16: return 151;
    case -15: return 150;
    case -14: return 148;
    case -13: return 147;
    case -12: return 145;
    case -11: return 144;
    case -10: return 142;
    case -9: return 141;
    case -8: return 140;
    case -7: return 138;
    case -6: return 137;
    case -5: return 135;
    case -4: return 134;
    case -3: return 132;
    case -2: return 131;
    case -1: return 129;
    case 0: return 128;
    case 1: return 127;
    case 2: return 125;
    case 3: return 124;
    case 4: return 122;
    case 5: return 121;
    case 6: return 119;
    case 7: return 118;
    case 8: return 116;
    case 9: return 115;
    case 10: return 114;
    case 11: return 112;
    case 12: return 111;
    case 13: return 109;
    case 14: return 108;
    case 15: return 106;
    case 16: return 105;
    case 17: return 104;
    case 18: return 102;
    case 19: return 101;
    case 20: return 99;
    case 21: return 98;
    case 22: return 97;
    case 23: return 95;
    case 24: return 94;
    case 25: return 93;
    case 26: return 91;
    case 27: return 90;
    case 28: return 89;
    case 29: return 87;
    case 30: return 86;
    case 31: return 85;
    case 32: return 83;
    case 33: return 82;
    case 34: return 81;
    case 35: return 79;
    case 36: return 78;
    case 37: return 77;
    case 38: return 76;
    case 39: return 74;
    case 40: return 73;
    case 41: return 72;
    case 42: return 71;
    case 43: return 70;
    case 44: return 68;
    case 45: return 67;
    case 46: return 66;
    case 47: return 65;
    case 48: return 64;
    case 49: return 63;
    case 50: return 61;
    case 51: return 60;
    case 52: return 59;
    case 53: return 58;
    case 54: return 57;
    case 55: return 56;
    case 56: return 55;
    case 57: return 54;
    case 58: return 53;
    case 59: return 52;
    case 60: return 51;
    case 61: return 50;
    case 62: return 49;
    case 63: return 48;
    case 64: return 47;
    case 65: return 46;
    case 66: return 45;
    case 67: return 44;
    case 68: return 43;
    case 69: return 42;
    case 70: return 41;
    case 71: return 40;
    case 72: return 39;
    case 73: return 39;
    case 74: return 38;
    case 75: return 37;
    case 76: return 36;
    case 77: return 35;
    case 78: return 35;
    case 79: return 34;
    case 80: return 33;
    case 81: return 32;
    case 82: return 32;
    case 83: return 31;
    case 84: return 30;
    case 85: return 29;
    case 86: return 29;
    case 87: return 28;
    case 88: return 27;
    case 89: return 27;
    case 90: return 26;
    case 91: return 25;
    case 92: return 25;
    case 93: return 24;
    case 94: return 24;
    case 95: return 23;
    case 96: return 22;
    case 97: return 22;
    case 98: return 21;
    case 99: return 21;
    case 100: return 20;
    case 101: return 20;
    case 102: return 19;
    case 103: return 19;
    case 104: return 18;
    case 105: return 18;
    case 106: return 17;
    case 107: return 17;
    case 108: return 16;
    case 109: return 16;
    case 110: return 15;
    case 111: return 15;
    case 112: return 14;
    case 113: return 14;
    case 114: return 14;
    case 115: return 13;
    case 116: return 13;
    case 117: return 13;
    case 118: return 12;
    case 119: return 12;
    case 120:
    case 121:
    case 122: return 11;
    case 123:
    case 124:
    case 125:
    case 126: return 10;
    case 127:
    case 128:
    case 129: return 9;
    case 130:
    case 131:
    case 132:
    case 133: return 8;
    case 134:
    case 135:
    case 136:
    case 137:
    case 138: return 7;
    case 139:
    case 140:
    case 141:
    case 142:
    case 143: return 6;
    case 144:
    case 145:
    case 146:
    case 147:
    case 148: return 5;
    case 149:
    case 150:
    case 151:
    case 152:
    case 153:
    case 154:
    case 155:
    case 156: return 4;
    case 157:
    case 158:
    case 159:
    case 160:
    case 161:
    case 162:
    case 163:
    case 164:
    case 165: return 3;
    case 166:
    case 167:
    case 168:
    case 169:
    case 170:
    case 171:
    case 172:
    case 173:
    case 174:
    case 175:
    case 176:
    case 177:
    case 178: return 2;
    default: return 1;
    }
}


// --------------- Pixel codes computing ---------------
float*** load_codes(char* phase_format, char* img_format, char numbered_imgs,
                    int nb_patterns, int nb_shifts, int w, int h) {
    float*** codes = malloc_f32cube(nb_patterns, w, h);
    char filename[FNAME_MAX_LEN];

    // Compute phases
    #pragma omp parallel for private(filename)
    for(int k=0; k<nb_patterns; k++) {
        sprintf(filename, phase_format, w, h, k);

        // If phase has already been computed, load it from the image file
        if(access(filename, F_OK) != -1) {
            codes[k] = load_gray(filename, &w, &h);

            for(int i=0; i < h; i++) {
                for(int j=0; j < w; j++) {
                    codes[k][i][j] = codes[k][i][j] / 65535.0 * 2 * PI - PI;
                }
            }
        } else {

            float*** intensities = malloc(sizeof(float**) * nb_shifts);
            float** image = malloc_f32matrix(w, h);

            #pragma omp parallel for
            for(int shift=0; shift < nb_shifts; shift++) {

                // FIXME : new image format
                if(numbered_imgs)
                    sprintf(filename, img_format, k * nb_shifts + shift);
                else
                    sprintf(filename, img_format, w, h, k, shift);

                intensities[shift] = load_gray(filename, &w, &h);

            }

            #pragma omp parallel for
            for(int i=0; i < h; i++) {
                for(int j=0; j < w; j++) {
                    float* intensities_ij = malloc(sizeof(float) * nb_shifts);

                    for(int shift=0; shift < nb_shifts; shift++) {
                        intensities_ij[shift] = intensities[shift][i][j];
                    }

                    codes[k][i][j] = solve_phase(intensities_ij, nb_shifts);
                    image[i][j] = (codes[k][i][j]  + PI) / (2 * PI) * 65535.0;

                    free(intensities_ij);
                }
            }

            sprintf(filename, phase_format, w, h, k);
            save_gray_png(filename, image, w, h, 16);

            free_f32matrix(image);
            free_f32cube(intensities, nb_shifts);
        }
    }

    return codes;
}

float** load_mask(char* img_format, int nb_patterns, int w, int h) {

    char* mask_fname = "maskCam.png";

    if(access(mask_fname, F_OK) != -1) {
        return load_gray(mask_fname, &w, &h);
    }

    float** mask = malloc_f32matrix(w, h);
    float** mask_max = malloc_f32matrix(w, h);
    float** mask_min = malloc_f32matrix(w, h);

    #pragma omp parallel for
    for(int i=0; i < h; i++)
        for(int j=0; j < w; j++)
            mask_min[i][j] = 255;

    // Compute mask
    #pragma omp parallel for
    for(int k=0; k<nb_patterns; k++) {

        char filename[FNAME_MAX_LEN];
        sprintf(filename, img_format, k);

        float** image = load_gray(filename, &w, &h);

        #pragma omp parallel for
        for(int i=0; i < h; i++) {
            for(int j=0; j < w; j++) {
                mask_min[i][j] = fmin(mask_min[i][j], image[i][j]);
                mask_max[i][j] = fmax(mask_max[i][j], image[i][j]);
            }
        }

        free_f32matrix(image);
    }

    #pragma omp parallel for
    for(int i=0; i < h; i++) {
        for(int j=0; j < w; j++) {
            mask[i][j] = (int) mask_max[i][j] - mask_min[i][j];
        }
    }

    save_gray_png(mask_fname, mask, w, h, 8);
    save_gray_png("camMin.png", mask_min, w, h, 8);
    save_gray_png("camMax.png", mask_max, w, h, 8);

    return mask;
}

float*** quadratic_codes(float*** orig_codes, int nb_patterns, int w, int h, int* new_nb_patterns) {

    int nb_extended = (nb_patterns + 1) * nb_patterns / 2;

    // Limit
    if(nb_extended > 300)
        nb_extended = 300;

    *new_nb_patterns = nb_extended;

    float*** codes = malloc_f32cube(nb_extended, w, h);

    #pragma omp parallel for
    for(int k=0; k < nb_patterns; k++)
        for(int i=0; i < h; i++)
            for(int j=0; j < w; j++)
                codes[k][i][j] = orig_codes[k][i][j];

    int k=nb_patterns;

    for(int i_ref=0; i_ref < nb_patterns; i_ref++) {
        for(int j_ref=i_ref + 1; j_ref < nb_patterns; j_ref++) {

            // Difference between pixel phases : [-2pi, 2pi] /2 => [-pi, pi]
            #pragma omp parallel for
            for(int i=0; i < h; i++)
                for(int j=0; j < w; j++)
                    codes[k][i][j] = (orig_codes[i_ref][i][j] - orig_codes[j_ref][i][j]) / 2.0;

            k++;

            if(k == nb_extended)
                goto end_copy;
        }
    }

end_copy:

    return codes;
}

/*
  Save a map (16-bits grayscale) of the phase for a set of shifted leopards
*/
void save_phase(float*** intensities, char* filename, int nb_shifts, int w, int h) {
    float** image = malloc_f32matrix(w, h);

    // Precompute phases
    #pragma omp parallel
    {
        float* intensities_ij = malloc(sizeof(float) * nb_shifts);

        for(int i=0; i < h; i++) {
            for(int j=0; j < w; j++) {

                for(int shift=0; shift < nb_shifts; shift++) {
                    intensities_ij[shift] = intensities[shift][i][j];
                }

                image[i][j] = (solve_phase(intensities_ij, nb_shifts) + PI) / (2 * PI) * 65535;
            }
        }

        free(intensities_ij);
    }

    save_gray_png(filename, image, w, h, 16);

    free_f32matrix(image);
}

/**
 * Remap an array of x/y/error to an array of 16-bit integers and save it
 */
void save_color_map(char* filename, float*** matches, int from_w, int from_h,
                    int to_w, int to_h, float max_distance) {
    float*** channels = malloc_f32cube(3, from_w, from_h);

    // Color map
    for(int i=0; i<from_h; i++)
        for(int j=0; j<from_w; j++) {
            if(matches[X][i][j] == -1.0) {
                channels[X][i][j] = 65535.0;
                channels[Y][i][j] = 65535.0;
                channels[DIST][i][j] = 65535.0;
            } else {
                channels[X][i][j] = matches[X][i][j] * 65535.0/to_w;
                channels[Y][i][j] = matches[Y][i][j] * 65535.0/to_h;
                channels[DIST][i][j] = fmin(matches[DIST][i][j] * 65535.0/max_distance, 65535.0);
            }
        }

    save_color_png(filename, channels, from_w, from_h, 16);

    free_f32cube(channels, 3);
}
