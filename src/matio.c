/*
   Helper methods to print 1D, 2D and 3D matrices to a file.

   Copyright 2015 Ian Kavanagh

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "matio.h"

void mprint(FILE *restrict stream, const int rank, const size_t *n, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x)) {
    size_t m = 1;

    fprintf(stream, "%d\n", rank);
    for (int i = 0; i < rank; ++i) {
        m *= n[i];
        fprintf(stream, "%d ", (int) n[i]);
    }
    fprintf(stream, "\n");

    for (size_t i = 0; i < m; ++i) {
        (*print)(stream, x + i * size);
    }
}

void mprintf_nm(FILE *restrict stream, const size_t n, const size_t m, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x)) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            // (col_index * row_size + row_index)
            size_t idx = i * m + j;

            (*print)(stream, x + idx * size);
        }
    }
}

void mprintf_3d(FILE *restrict stream, const size_t n, const size_t m, const size_t o, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x), const char *label) {
    for (size_t i = 0; i < n; ++i) {
        size_t idx = i*m*o;

        fprintf(stream, "%s(%d, :, :) = [\n", label, (int) i);
        mprintf_nm(stream, m, o, x + idx * size, size, (*print));
        fprintf(stream, "]\n");
    }
}

void mprintf_2d(FILE *restrict stream, const size_t n, const size_t m, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x), const char *label) {
    fprintf(stream, "%s = [\n", label);
    mprintf_nm(stream, n, m, x, size, (*print));
    fprintf(stream, "]\n");
}

void mprintf_1d(FILE *restrict stream, const size_t n, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x), const char *label) {
    mprintf_2d(stream, n, 1, x, size, (*print), label);
}

void mprintf(FILE *restrict stream, const int rank, const size_t *n, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x), const char *label) {
    switch (rank) {
        case 3:
            mprintf_3d(stream, n[0], n[1], n[2], x, size, (*print), label);
            break;
        case 2:
            mprintf_2d(stream, n[0], n[1], x, size, (*print), label);
            break;
        case 1:
            mprintf_1d(stream, n[0], x, size, (*print), label);
            break;
        default:
            fprintf(stderr, "Error: rank not one of (1, 2, 3)!\n");
            fprintf(stderr, "\tDefaulting to rank = 1.");
            size_t m = n[0];
            for (int i = 1; i < rank; ++i) {
                m += n[i];
            }
            mprintf_1d(stream, m, x, size, (*print), label);
    }
}

inline void printi(FILE *restrict stream, const void *x) {
    const int *i = x;

    fprintf(stream, "%d\n", (*i));
}

inline void printm(FILE *restrict stream, const void *x) {
    const material_t *m = x;

    fprintf(stream, "%.4f %.4f %.4f\n", (*m).epsilon_r, (*m).mu_r, (*m).sigma);
}

inline void printz(FILE *restrict stream, const void *x) {
    const double complex *z = x;

    double re = creal(*z), im = cimag(*z);
    fprintf(stream, "%6.4f%c%6.4fi\n", re, im > 0 ? '+' : '-', fabs(im));
}