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

#include <math.h>

void imprint(FILE *restrict stream, const int rank, const int *n, const int *restrict x) {
    // TODO: Lump imprint, mmprint and zmprint together so only 1 method is needed to
    //       print different types of matrices.
    int m = 1;

    fprintf(stream, "%d\n", rank);
    for (int i = 0; i < rank; ++i) {
        m *= n[i];
        fprintf(stream, "%d ", n[i]);
    }
    fprintf(stream, "\n");

    for (int i = 0; i < m; ++i) {
        fprintf(stream, "%d\n", x[i]);
    }
}

void mmprint(FILE *restrict stream, const int rank, const int *n, const material_t *restrict x) {
    int m = 1;

    fprintf(stream, "%d\n", rank);
    for (int i = 0; i < rank; ++i) {
        m *= n[i];
        fprintf(stream, "%d ", n[i]);
    }
    fprintf(stream, "\n");

    for (int i = 0; i < m; ++i) {
        material_t material = x[i];

        fprintf(stream, "%.4f %.4f %.4f\n", material.epsilon_r, material.mu_r, material.sigma);
    }
}

void zmprint(FILE *restrict stream, const int rank, const int *n, const double complex *restrict x) {
    int m = 1;

    fprintf(stream, "%d\n", rank);
    for (int i = 0; i < rank; ++i) {
        m *= n[i];
        fprintf(stream, "%d ", n[i]);
    }
    fprintf(stream, "\n");

    for (int i = 0; i < m; ++i) {
        double re = creal(x[i]), im = cimag(x[i]);

        fprintf(stream, "%6.4f%c%6.4fi\n", re, im > 0 ? '+' : '-', fabs(im));
    }
}

void zfprint_nm(FILE *restrict stream, const int n, const int m, const double complex *restrict x) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            // (col_index * row_size + row_index)
            int idx = i * m + j;
            double re = creal(x[idx]), im = cimag(x[idx]);

            fprintf(stream, "%6.4f %c %6.4fi ", re, im > 0 ? '+' : '-', fabs(im));
        }
        fprintf(stream, "\n");
    }
}

void zfprint_3d(FILE *restrict stream, const int n, const int m, const int o, const double complex *restrict x, const char *label) {
    for (int i = 0; i < n; ++i) {
        fprintf(stream, "%s(%d, :, :) = [\n", label, i);
        zfprint_nm(stream, m, o, x + i*m*o);
        fprintf(stream, "]\n");
    }
}

void zfprint_2d(FILE *restrict stream, const int n, const int m, const double complex *restrict x, const char *label) {
    fprintf(stream, "%s = [\n", label);
    zfprint_nm(stream, n, m, x);
    fprintf(stream, "]\n");
}

void zfprint_1d(FILE *restrict stream, const int n, const double complex *restrict x, const char *label) {
    zfprint_2d(stream, n, 1, x, label);
}

void zfprint(FILE *restrict stream, const int rank, const int *n, const double complex *restrict x, const char *label) {
    switch(rank) {
        case 3:
            zfprint_3d(stream, n[0], n[1], n[2], x, label);
            break;
        case 2:
            zfprint_2d(stream, n[0], n[1], x, label);
            break;
        case 1:
            zfprint_1d(stream, n[0], x, label);
        default:
            fprintf(stderr, "Error: rank not one of (1, 2, 3)!\n");
            fprintf(stderr, "\tDefaulting to rank = 1.");
            int m = n[0];
            for (int i = 1; i < rank; ++i) {
                m += n[i];
            }
            zfprint_1d(stream, m, x, label);
    }
}

void zprint_3d(const int n, const int m, const int o, const double complex *restrict x, const char *label) {
    zfprint_3d(stdout, n, m, o, x, label);
}

void zprint_2d(const int n, const int m, const double complex *restrict x, const char *label) {
    zfprint_2d(stdout, n, m, x, label);
}

void zprint_1d(const int n, const double complex *restrict x, const char *label) {
    zfprint_1d(stdout, n, x, label);
}

void zprint(const int rank, const int *n, const double complex *restrict x, const char *label) {
    zfprint(stdout, rank, n, x, label);
}
