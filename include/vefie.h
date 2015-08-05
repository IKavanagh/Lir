/*
   Methods related to creating and solving the VEFIE. 

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
#include <fftw3.h>

/**
 * Purpose
 * =======
 *
 * The byte alignment used with the Intel MKL library to ensure similarity
 * between runs.
 *
 * =============================================================================
 */
extern int alignment;

/**
 * Purpose
 * =======
 *
 * The matrices used in the VEFIE when computing an indoor propagation problem.
 *
 * Matrices
 *                D: Material contrast matrix
 *                G: Greens function
 *                V: Incident field
 *
 * =============================================================================
 */
extern double complex *D, *G, *V;

/**
 * Purpose
 * =======
 *
 * Total electric field. Initial guess on input to bicgstab().
 *
 * =============================================================================
 */
extern double complex *E;

/**
 * Purpose
 * =======
 *
 * Defines the positions at which the field points are solved for.
 *
 * =============================================================================
 */
extern double complex *p;

/**
 * Purpose
 * =======
 *
 * A reduced operator to apply to the BiCGSTAB to focus it's solution.
 *
 * =============================================================================
 */
extern double complex *rfo;

/**
 * Purpose
 * =======
 *
 * The block size to be used in computing block preconditioners. Specifically
 * for computing block preconditioners.
 *
 * =============================================================================
 */
extern int block_size;

/**
 * Purpose
 * =======
 *
 * Initialises the VEFIE by creating the vectors, V, D and G. V is the incident
 * electric field, D is the material contrast vector and G is FFT of the
 * block-Toeplitz Green's function.
 *
 * V and D are of length nm, whereas G is a matrix of size 2n x 2m.
 *
 * Defines p to be the positions of the points in shape.
 *
 * Initialises the remaining variables needed for the solution of the VEFIE
 * such as those required in computing the FFT.
 *
 * Arguments
 * =========
 *
 * f        (input) DOUBLE.
 *          The frequency the VEFIE is being run for. Used to determine the
 *          number of discretised points needed.
 * 
 * antenna  (input) DOUBLE COMPLEX.
 *          The location of the antenna within the problem.
 *          The real part denotes the x-coordinate and the imaginary part the
 *          y-coordinate.
 * 
 * inc      (external subroutine)
 *          An external user supplied subroutine compute the incident field at
 *          the point p a distance R from the antenna with the free space wave
 *          number k0.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  -1: Unable to allocate memory.
 *
 *          =  -2: An error occurred whilst planning the FFT.
 *
 * =============================================================================
 */
int init_vefie(const double f, const double complex antenna, double complex (*inc)(const double k0, const double complex antenna, const double complex point));

/**
 * Purpose
 * =======
 *
 * Initialises the dimensions for the FFT, allocates memory for the input and
 * output arrays and plans the forward and backwards FFTs.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  -1: Unable to allocate memory.
 *
 *          =  -2: An error occurred whilst planning the FFT.
 *
 * =============================================================================
 */
int init_fftw(void);

/**
 * Purpose
 * =======
 *
 * Computes the FFT of a block Toeplitz matrix when given the first row or
 * column of that matrix.
 *
 * Returns the FFT of the block Toeplitz matrix in x.
 *
 * Arguments
 * =========
 * 
 * x        (input/output) DOUBLE COMPLEX array, dimension
 *          (nm, 1).
 *          On input the first row or column of the block Toeplitz matrix which
 *          is to be extended to size 2n x 2m.
 *          On output the FFT of the extended block Toeplitz matrix of
 *          dimensions 2n x 2m.
 *
 * =============================================================================
 */
void toeplitz(double complex** restrict x);

/**
 * Purpose
 * =======
 *
 * The subroutine which performs the matrix vector product
 *
 *               y := alpha*A*x + beta*y,
 *
 * in the bicgstab() for a two-dimensional problem.
 *
 * For the 2D VEFIE
 *
 *               A = I + GD
 *
 * Arguments
 * =========
 *
 * alpha    (input) DOUBLE COMPLEX.
 *          Constant to multiply the matrix vector product of A*x by.
 *
 * x        (input) DOUBLE COMPLEX array, dimension(N, 1).
 *          The input vector which is to be multiplied by the system matrix A.
 *
 * beta     (input) DOUBLE COMPLEX.
 *          Constant to multiply the vector y by before adding it to the product of alpha*A*x.
 *
 * y        (input/output) DOUBLE COMPLEX array, dimension(N, 1).
 *          On input a vector which is to be added to the product of alpha*A*x.
 *          On output contains the result
 *
 *               y := alpha*A*x + beta*y
 *
 * =============================================================================
 */
void matvec(const double complex *restrict alpha, const double complex *restrict x, const double complex *restrict beta, double complex *restrict y);

/**
 * Purpose
 * =======
 *
 * Dynamically allocates memory for all of the matrices used in the VEFIE.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int allocate_matrices(void);

/**
 * Purpose
 * =======
 *
 * Dynamically allocates memory for a matrix and optionally initialises it to 0.
 *
 * Arguments
 * =========
 *
 * N        (input) SIZE_T.
 *          The number of elements to dynamically allocate in the matrix.
 * 
 * init     (input) INTEGER.
 *          Whether the allocated matrix should be initialised to 0 or not.
 *          Setting init to 1 will initialise the allocated matrix to 0.
 * 
 * X        (output) DOUBLE COMPLEX array.
 *          The dynamically allocated array with N elements.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int matalloc(const size_t N, const int init, double complex **restrict X);

/**
 * Purpose
 * =======
 *
 * Cleanup any dynamically allocated memory.
 *
 * =============================================================================
 */
void vefie_cleanup(void);
