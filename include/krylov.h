/*
   Solve the standard matrix equation Ax = b using an iterative solver
   from the Krylov family of iterative solvers.

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
 * Solves the linear system Ax = b where A = I + GD using the BiConjugate
 * Gradient Stabilized (BiCGSTAB) iterative method with preconditioning.
 *
 * Convergence test: (norm(b - Ax) / norm(b)) < tol.
 *
 * Arguments
 * =========
 *
 * n        (input) INTEGER.
 *          Specifies the dimension of vectors.
 *
 * b        (input) DOUBLE COMPLEX array, dimension (n, 1).
 *          Right hand side vector b.
 * 
 * x        (input/output) DOUBLE COMPLEX array, dimension (n, 1).
 *          On input the initial guess. On output the final solution of Ax = b.
 * 
 * work     (input) DOUBLE COMPLEX array, dimension (ldw, 7).
 *          Workspace for resid, direction vector, etc.
 *          Vectors r and s share the same workspace.
 *
 * ldw      (input) INTEGER.
 *          Leading dimension of workspace.
 *          ldw >= max(1, n)
 * 
 * iter     (input/output) INTEGER.
 *          On input the maximum number of iterations to be performed. On output
 *          the number of iterations performed by the solver.
 *
 * resid    (input/output) INTEGER.
 *          On input the measure for the convergence test. On output the final
 *          value of norm(b - Ax) / norm(b).
 * 
 * matvec   (external subroutine)
 *          An external user supplied subroutine to perform the matrix-vector
 *          product
 *
 *               y := alpha*A*x + beta*y,
 *
 *          where alpha and beta are scalars, x and y are vectors, and A is a
 *          matrix. Vector x must remain unchanged. The solution is over-written
 *          on vector y.
 * 
 *          The call is:
 * 
 *              matvec(alpha, x, beta, y)
 *
 * pre      (external subroutine)
 *          An external user supplied subroutine to apply the preconditioner for
 *          the system.
 *
 *               M*x = b,
 *
 *          where x and b are vectors, and M a matrix. Vector x must remain
 *          unchanged. The solution is over-written on vector b.
 *
 *          The call is:
 *
 *             psolve(x, b)
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Converged successfully.
 *
 *          >  0: Convergence to tolerance not achieved.
 *
 *          <  0: Illegal input parameter or parameter breakdown occurred.
 *
 *                Illegal parameter:
 *
 *                   -1: matrix dimension n < 0
 *                   -2: ldw < n
 *                   -3: Maximum number of iterations iter <= 0.
 *
 *                Breakdown: rho or omega are approximately equal to 0.
 *
 *                  -10: rho ~= 0: rho and rtld have become orthogonal.
 *                  -11: omega ~= 0: s and t have become orthogonal relative to
 *                                   t'*t.
 *
 * =============================================================================
 */
int bicgstab(const int n, const double complex *restrict b, double complex *restrict x, double complex *restrict work, const int ldw, int *restrict iter, double *restrict resid, void (*matvec)(const double complex *restrict alpha, const double complex *restrict x, const double complex *restrict beta, double complex *restrict y), void (*pre)(double complex *restrict x, const double complex *restrict b));

/**
 * Purpose
 * =======
 *
 * Solves the linear system Ax = b where A = I + GD using the BiConjugate
 * Gradient Stabilized (BiCGSTAB) iterative method with preconditioning and an
 * additional reduced operator. The reduced operator should be used to ignore
 * particular values in the BiCGSTAB iteration so only a reduced number of more
 * important unknown values are solved within the iteration and the remaining
 * values are computed at the end. 
 *
 * Convergence test: (norm(b - Ax) / norm(b)) < tol.
 *
 * Arguments
 * =========
 *
 * n        (input) INTEGER.
 *          Specifies the dimension of vectors.
 *
 * b        (input) DOUBLE COMPLEX array, dimension (n, 1).
 *          Right hand side vector b.
 * 
 * x        (input/output) DOUBLE COMPLEX array, dimension (n, 1).
 *          On input the initial guess. On output the final solution of Ax = b.
 * 
 * work     (input) DOUBLE COMPLEX array, dimension (ldw, 8).
 *          Workspace for resid, direction vector, etc.
 *          Vectors r and s share the same workspace.
 *
 * ldw      (input) INTEGER.
 *          Leading dimension of workspace.
 *          ldw >= max(1, n)
 * 
 * iter     (input/output) INTEGER.
 *          On input the maximum number of iterations to be performed. On output
 *          the number of iterations performed by the solver.
 *
 * resid    (input/output) INTEGER.
 *          On input the measure for the convergence test. On output the final
 *          value of norm(b - Ax) / norm(b).
 * 
 * matvec   (external subroutine)
 *          An external user supplied subroutine to perform the matrix-vector
 *          product
 *
 *               y := alpha*A*x + beta*y,
 *
 *          where alpha and beta are scalars, x and y are vectors, and A is a
 *          matrix. Vector x must remain unchanged. The solution is over-written
 *          on vector y.
 * 
 *          The call is:
 * 
 *              matvec(alpha, x, beta, y)
 *
 * pre      (external subroutine)
 *          An external user supplied subroutine to apply the preconditioner for
 *          the system.
 *
 *               M*x = b,
 *
 *          where x and b are vectors, and M a matrix. Vector x must remain
 *          unchanged. The solution is over-written on vector b.
 *
 *          The call is:
 *
 *             psolve(x, b)
 * 
 * rfo      (input) DOUBLE COMPLEX array, dimension (n, 1).
 *          Reduced operator used to ignore particular unknown values in the
 *          BiCGSTAB iteration. If no reduced operator is to be used NULL can
 *          be passed in instead.
 *          The corresponding element in rfo for values to be ignored should be
 *          set to 0, whilst it should be set to 1 for all values to be used in
 *          the iteration process.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Converged successfully.
 *
 *          >  0: Convergence to tolerance not achieved.
 *
 *          <  0: Illegal input parameter or parameter breakdown occurred.
 *
 *                Illegal parameter:
 *
 *                   -1: matrix dimension n < 0
 *                   -2: ldw < n
 *                   -3: Maximum number of iterations iter <= 0.
 *
 *                Breakdown: rho or omega are approximately equal to 0.
 *
 *                  -10: rho ~= 0: rho and rtld have become orthogonal.
 *                  -11: omega ~= 0: s and t have become orthogonal relative to
 *                                   t'*t.
 *
 * =============================================================================
 */
int rbicgstab(const int n, const double complex *restrict b, double complex *restrict x, double complex *restrict work, const int ldw, int *restrict iter, double *restrict resid, void (*matvec)(const double complex *restrict alpha, const double complex *restrict x, const double complex *restrict beta, double complex *restrict y), void (*pre)(double complex *restrict x, const double complex *restrict b), const double complex *restrict rfo);

/**
 * Purpose
 * =======
 *
 * Dummy routine which can be used in the BiCGSTAB if a preconditioner is not to
 * be used.
 *
 * Performs a copy of the form
 *
 *      x := b
 *
 * Arguments
 * =========
 * 
 * x        (input/output) DOUBLE COMPLEX array, dimension (n, 1).
 *          On output the preconditioned vector.
 *
 * b        (input) DOUBLE COMPLEX array, dimension (n, 1).
 *          The right hand side vector in Ax = b which is to be preconditioned
 *          by M.
 *
 * =============================================================================
 */
void no_pre(double complex *restrict x, const double complex *restrict b);

/**
 * Purpose
 * =======
 *
 * Prints out a status message about the bicgstab() result.
 *
 * Arguments
 * =========
 *
 * info     (input) INTEGER.
 *          The result code of bicgstab().
 *
 * iter     (input) INTEGER.
 *          The number of iterations as returned by bicgstab().
 *
 * t        (input) DOUBLE.
 *          The length of time taken to by the bicgstab().
 *
 * =============================================================================
 */
void iprint(int info, int iter, double t);
