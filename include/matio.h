/*
   Helper methods to print 1D, 2D and 3D matrices of varying types to a file.

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

#pragma once

#include <complex.h>
#include <stdio.h>

#include "physics.h"

/**
 * Purpose
 * =======
 *
 * Print an n-dimensional integer matrix to stream.
 *
 * The first line contains the rank of the matrix. The second line contains the
 * length of each dimension.
 *
 * The matrix output begins on the third line with one entry per line.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * rank     (input) INTEGER.
 *          The rank of the matrix to be printed.
 * 
 * n        (input) INTEGER array.
 *          The dimensions of the matrix.
 * 
 * x        (input) INTEGER array, dimension (n[0], n[1], ...).
 *          The matrix to be printed.
 *
 * =============================================================================
 */
void imprint(FILE *restrict stream, const int rank, const int *n, const int *restrict x);

/**
 * Purpose
 * =======
 *
 * Print an n-dimensional material matrix to stream.
 *
 * The first line contains the rank of the matrix. The second line contains the
 * length of each dimension.
 *
 * The matrix output begins on the third line with one entry per line.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * rank     (input) INTEGER.
 *          The rank of the matrix to be printed.
 * 
 * n        (input) INTEGER array.
 *          The dimensions of the matrix.
 * 
 * x        (input) MATERIAL_T array, dimension (n[0], n[1], ...).
 *          The matrix to be printed.
 *
 * =============================================================================
 */
void mmprint(FILE *restrict stream, const int rank, const int *n, const material_t *restrict x);

/**
 * Purpose
 * =======
 *
 * Print an n-dimensional complex matrix to stream.
 *
 * The first line contains the rank of the matrix. The second line contains the
 * length of each dimension.
 *
 * The matrix output begins on the third line with one entry per line.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * rank     (input) INTEGER.
 *          The rank of the matrix to be printed.
 * 
 * n        (input) INTEGER array.
 *          The dimensions of the matrix.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension (n[0], n[1], ...).
 *          The matrix to be printed.
 *
 * =============================================================================
 */
void zmprint(FILE *restrict stream, const int rank, const int *n, const double complex *restrict x);

/**
 * Purpose
 * =======
 *
 * Print a three-dimensional complex matrix to stream.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * n        (input) INTEGER.
 *          The number of rows of the matrix.
 * 
 * m        (input) INTEGER.
 *          The number of columns of the matrix.
 * 
 * o        (input) INTEGER.
 *          The depth of the matrix.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension (n, m, o).
 *          The matrix to be printed.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void zfprint_3d(FILE *restrict stream, const int n, const int m, const int o, const double complex *restrict x, const char *label);

/**
 * Purpose
 * =======
 *
 * Print a two-dimensional complex matrix to stream.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * n        (input) INTEGER.
 *          The number of rows of the matrix.
 * 
 * m        (input) INTEGER.
 *          The number of columns of the matrix.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension (n, m).
 *          The matrix to be printed.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void zfprint_2d(FILE *restrict stream, const int n, const int m, const double complex *restrict x, const char *label);

/**
 * Purpose
 * =======
 *
 * Print a one-dimensional complex matrix (vector) to stream.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * n        (input) INTEGER.
 *          Length of the vector.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension (n, 1).
 *          The matrix to be printed.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void zfprint_1d(FILE *restrict stream, const int n, const double complex *restrict x, const char *label);

/**
 * Purpose
 * =======
 *
 * Prints an n-dimensional complex matrix to stream.
 *
 * NOTE: Only 1, 2 and 3 dimensions have been implemented!
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * rank     (input) INTEGER.
 *          The rank of the matrix to be printed.
 * 
 * n        (input) INTEGER array.
 *          The dimensions of the matrix.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension
 *          (n[0], n[1], ...).
 *          The matrix to be printed.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void zfprint(FILE *restrict stream, const int rank, const int *n, const double complex* restrict x, const char *label);

/**
 * Purpose
 * =======
 *
 * Prints an n-dimensional complex matrix to stdout. Only
 * 1, 2 and 3 dimensions have been implemented.
 *
 * Arguments
 * =========
 * 
 * rank     (input) INTEGER.
 *          The rank of the matrix to be printed.
 * 
 * n        (input) INTEGER array.
 *          The dimensions of the matrix.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension
 *          (n[0], n[1], ...).
 *          The matrix to be printed.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void zprint(const int rank, const int *n, const double complex *restrict x, const char *label);

/**
 * Purpose
 * =======
 *
 * Print a three-dimensional complex matrix to stdout.
 *
 * Arguments
 * =========
 * 
 * n        (input) INTEGER.
 *          The number of rows of the matrix.
 * 
 * m        (input) INTEGER.
 *          The number of columns of the matrix.
 * 
 * o        (input) INTEGER.
 *          The depth of the matrix.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension (n, m, o).
 *          The matrix to be printed.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void zprint_3d(const int n, const int m, const int o, const double complex *restrict x, const char *label);

/**
 * Purpose
 * =======
 *
 * Print a two-dimensional complex matrix to stdout.
 *
 * Arguments
 * =========
 * 
 * n        (input) INTEGER.
 *          The number of rows of the matrix.
 * 
 * m        (input) INTEGER.
 *          The number of columns of the matrix.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension (n, m).
 *          The matrix to be printed.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void zprint_2d(const int n, const int m, const double complex *restrict x, const char *label);

/**
 * Purpose
 * =======
 *
 * Print a one-dimensional complex matrix (vector) to stdout.
 *
 * Arguments
 * =========
 * 
 * n        (input) INTEGER.
 *          The length of the vector.
 * 
 * x        (input) DOUBLE COMPLEX array, dimension (n, 1).
 *          The matrix to be printed.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void zprint_1d(const int n, const double complex *restrict x, const char *label);
