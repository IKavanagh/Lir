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
 * Print an n-dimensional matrix to stream.
 *
 * The first line contains the rank of the matrix. The second line contains the
 * length of each dimension.
 *
 * The matrix output begins on the third line with one entry per line.
 *
 * A suitable print function should be supplied to print 1 element of the array
 * per line.
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
 * n        (input) SIZE_T array.
 *          The dimensions of the matrix.
 * 
 * x        (input) VOID array, dimension (n[0], n[1], ...).
 *          The matrix to be printed.
 * 
 * size     (input) SIZE_T array.
 *          The size of 1 element of the array.
 *
 * print    (external subroutine)
 *          An external user supplied subroutine to print 1 element of the array
 *          per line.
 *
 * =============================================================================
 */
void mprint(FILE *restrict stream, const size_t rank, const size_t *n, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x));

/**
 * Purpose
 * =======
 *
 * Print a three-dimensional matrix to stream.
 *
 * A suitable print function should be supplied to print 1 element of the array
 * per line.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * n        (input) SIZE_T.
 *          The number of rows of the matrix.
 * 
 * m        (input) SIZE_T.
 *          The number of columns of the matrix.
 * 
 * o        (input) SIZE_T.
 *          The depth of the matrix.
 * 
 * x        (input) VOID array, dimension (n, m, o).
 *          The matrix to be printed.
 * 
 * size     (input) SIZE_T array.
 *          The size of 1 element of the array.
 *
 * print    (external subroutine)
 *          An external user supplied subroutine to print 1 element of the array
 *          per line.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void mprintf_3d(FILE *restrict stream, const size_t o, const size_t m, const size_t n, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x), const char *label);

/**
 * Purpose
 * =======
 *
 * Print a two-dimensional to stream.
 *
 * A suitable print function should be supplied to print 1 element of the array
 * per line.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * n        (input) SIZE_T.
 *          The number of rows of the matrix.
 * 
 * m        (input) SIZE_T.
 *          The number of columns of the matrix.
 * 
 * x        (input) VOID array, dimension (n, m).
 *          The matrix to be printed.
 * 
 * size     (input) SIZE_T array.
 *          The size of 1 element of the array.
 *
 * print    (external subroutine)
 *          An external user supplied subroutine to print 1 element of the array
 *          per line.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void mprintf_2d(FILE *restrict stream, const size_t m, const size_t n, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x), const char *label);

/**
 * Purpose
 * =======
 *
 * Print a one-dimensional complex matrix (vector) to stream.
 *
 * A suitable print function should be supplied to print 1 element of the array
 * per line.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the matrix is to be printed to.
 * 
 * n        (input) SIZE_T.
 *          Length of the vector.
 * 
 * x        (input) VOID array, dimension (n, 1).
 *          The matrix to be printed.
 * 
 * size     (input) SIZE_T array.
 *          The size of 1 element of the array.
 *
 * print    (external subroutine)
 *          An external user supplied subroutine to print 1 element of the array
 *          per line.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void mprintf_1d(FILE *restrict stream, const size_t n, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x), const char *label);

/**
 * Purpose
 * =======
 *
 * Prints an n-dimensional to stream.
 *
 * A suitable print function should be supplied to print 1 element of the array
 * per line.
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
 * n        (input) SIZE_T array.
 *          The dimensions of the matrix.
 * 
 * x        (input) VOID array, dimension (n[0], n[1], ...).
 *          The matrix to be printed.
 * 
 * size     (input) SIZE_T array.
 *          The size of 1 element of the array.
 *
 * print    (external subroutine)
 *          An external user supplied subroutine to print 1 element of the array
 *          per line.
 *
 * label    (input) STRING.
 *          A name to give the matrix when printing it.
 *
 * =============================================================================
 */
void mprintf(FILE *restrict stream, const size_t rank, const size_t *n, const void *restrict x, const size_t size, void (*print)(FILE *restrict stream, const void *x), const char *label);

/**
 * Purpose
 * =======
 *
 * Prints the element pointed to by x which is assumed to be an integer on 1
 * line in the stream.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the element is to be printed to.
 * 
 * x        (input) VOID pointer.
 *          A void pointer to the element that should be printed. This is cast
 *          to an int within the function.
 *
 * =============================================================================
 */
void printi(FILE *restrict stream, const void *x);

/**
 * Purpose
 * =======
 *
 * Prints the element pointed to by x which is assumed to be a double on 1
 * line in the stream.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the element is to be printed to.
 * 
 * x        (input) VOID pointer.
 *          A void pointer to the element that should be printed. This is cast
 *          to a double within the function.
 *
 * =============================================================================
 */
void printd(FILE *restrict stream, const void *x);

/**
 * Purpose
 * =======
 *
 * Prints the element pointed to by x which is assumed to be a material on 1
 * line in the stream.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the element is to be printed to.
 * 
 * x        (input) VOID pointer.
 *          A void pointer to the element that should be printed. This is cast
 *          to a material within the function.
 *
 * =============================================================================
 */
void printm(FILE *restrict stream, const void *x);


/**
 * Purpose
 * =======
 *
 * Prints the element pointed to by x which is assumed to be a complex double on
 * 1 line in the stream.
 *
 * Arguments
 * =========
 *
 * stream   (input) FILE pointer.
 *          A file pointer for the file that the element is to be printed to.
 * 
 * x        (input) VOID pointer.
 *          A void pointer to the element that should be printed. This is cast
 *          to a complex double within the function.
 *
 * =============================================================================
 */
void printz(FILE *restrict stream, const void *x);
