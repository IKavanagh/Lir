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

#include "krylov.h"

#include <fftw3.h>
#include <mkl.h>

/**
 * Purpose
 * =======
 *
 * Size of the problem being solved using the iterative solver,
 * i.e. the length of each vector in the iterative solver.
 *
 * ============================================================
 */
int N = 0;

/**
 * Purpose
 * =======
 *
 * The array increment value used for most vectors. Specifies
 * the number of elements in an array between consecutive items.
 * 
 * ============================================================
 */
static const int inc = 1;

/**
 * Purpose
 * =======
 *
 * The tolerance used when determining if parameter breakdown in
 * the bicgstab() has occurred.
 * 
 * Parameter's whose value is less than this will be deemed to
 * have reached breakdown.
 * 
 * ============================================================
 */
static const double partol = 1e-15;

/**
 * Purpose
 * =======
 *
 * Used in the bicgstab() when computing matrix vector products
 * with matvec().
 * 
 * ============================================================
 */
static const double complex one = 1.0;
static const double complex minus_one = -1.0;
static const double complex zero = 0.0;

inline int max(int a, int b) {
    if (b > a) return b;
    return a;
}

int bicgstab(const int n, const double complex *restrict b, double complex *restrict x, double complex *restrict work, const int ldw, int *restrict iter, double *restrict resid, void (*matvec)(const double complex *restrict alpha, const double complex *restrict x, const double complex *restrict beta, double complex *restrict y), void (*pre)(double complex *restrict x, const double complex *restrict b)) {
    return rbicgstab(n, b, x, work, ldw, iter, resid, (*matvec), (*pre), NULL);
}

int rbicgstab(const int n, const double complex *restrict b, double complex *restrict x, double complex *restrict work, const int ldw, int *restrict iter, double *restrict resid, void (*matvec)(const double complex *restrict alpha, const double complex *restrict x, const double complex *restrict beta, double complex *restrict y), void (*pre)(double complex *restrict x, const double complex *restrict b), const double complex *restrict rfo) {
    int r, rtld, p, v, t, phat, shat, neg, red, s;
    double complex alpha, beta, omega, rho, rho1 = 1.0, temp;

    // Provides more readable workspace indexing
    work -= ldw;

    if (n < 0) {
        return -1;
    } else if (ldw < max(1, n)) {
        return -2;
    } else if (*iter <= 0) {
        return -3;
    }

    N = n;
    
    int maxit = *iter;
    double tol = *resid;

    // Workspace columns
    r = 1;
    rtld = 2;
    p = 3;
    v = 4;
    t = 5;
    phat = 6;
    shat = 7;
    neg = 8;
    red = 9;
    s = 1;

    if (rfo) {
        // Flip reduced operator
        for (int i = 0; i < n; ++i) {
            if (creal(rfo[i]) > partol) {
                work[neg * ldw + i] = 0.0;
            } else {
                work[neg * ldw + i] = 1.0;
            }
        }
    }

    // Initial residual
    cblas_zcopy(n, &b[0], inc, &work[r * ldw], inc);
    if (cblas_dznrm2(n, &x[0], inc) > partol) {
        (*matvec)(&minus_one, &x[0], &one, &work[r * ldw]); // r = -A*x + b
        if (cblas_dznrm2(n, &work[r * ldw], inc) <= tol) {
            // Initial guess is less than tolerance
            return 0;
        }
    }
    if (rfo) {
        cblas_zgbmv(CblasRowMajor, CblasNoTrans, n, n, 0, 0, &one, rfo, inc, &work[r * ldw], inc, &zero, &work[rtld * ldw], inc); // Reduced operator
        cblas_zcopy(n, &work[rtld * ldw], inc, &work[r * ldw], inc);
    } else {
        cblas_zcopy(n, &work[r * ldw], inc, &work[rtld * ldw], inc);
    }

    double bnrm2 = cblas_dznrm2(n, &work[r * ldw], inc);
    if (bnrm2 < partol) {
        bnrm2 = 1.0;
    }

    *iter = 0;
    *resid = 1;

    while (*resid > tol && *iter < maxit) {
        ++(*iter);

        cblas_zdotc_sub(n, &work[rtld * ldw], inc, &work[r * ldw], inc, &rho); // rho = <rtld, r>
        if (cabs(rho) < partol) {
            return -10;
        }

        if (*iter > 1) {
            beta = (rho / rho1) * (alpha / omega); // beta = rho / rho1 * (alpha / omega)

            temp = -omega;
            cblas_zaxpy(n, &temp, &work[v * ldw], inc, &work[p * ldw], inc); // p = p - w*v
            cblas_zscal(n, &beta, &work[p * ldw], inc); // p = beta*(p - w*v)
            cblas_zaxpy(n, &one, &work[r * ldw], inc, &work[p * ldw], inc); // p = r + beta*(p - w*v)
        } else {
            cblas_zcopy(n, &work[r * ldw], inc, &work[p * ldw], inc);
        }

        (*pre)(&work[phat * ldw], &work[p * ldw]);
        (*matvec)(&one, &work[phat * ldw], &zero, &work[v * ldw]); // v = A*phat
        if (rfo) { // Reduced operator
            cblas_zgbmv(CblasRowMajor, CblasNoTrans, n, n, 0, 0, &one, rfo, inc, &work[v * ldw], inc, &zero, &work[red * ldw], inc);
            cblas_zcopy(n, &work[red * ldw], inc, &work[v * ldw], inc);
        }

        cblas_zdotc_sub(n, &work[rtld * ldw], inc, &work[v * ldw], inc, &alpha); // alpha = <rtld, v>
        alpha = rho / alpha; // alpha = rho / <rtld, v>

        temp = -alpha;
        cblas_zaxpy(n, &temp, &work[v * ldw], inc, &work[r * ldw], inc); // s = r - alpha*v
        if (cblas_dznrm2(n, &work[s * ldw], inc) <= tol) {
            cblas_zaxpy(n, &alpha, &work[phat * ldw], inc, &x[0], inc); // x = x + alpha*p
            *resid = cblas_dznrm2(n, &work[s * ldw], inc) / bnrm2;
            break;
        }

        (*pre)(&work[shat * ldw], &work[s * ldw]);
        (*matvec)(&one, &work[shat * ldw], &zero, &work[t * ldw]); // t = A*shat
        if (rfo) { // Reduced operator
            cblas_zgbmv(CblasRowMajor, CblasNoTrans, n, n, 0, 0, &one, rfo, inc, &work[t * ldw], inc, &zero, &work[red * ldw], inc); // Reduced operator
            cblas_zcopy(n, &work[red * ldw], inc, &work[t * ldw], inc);
        }

        cblas_zdotc_sub(n, &work[t * ldw], inc, &work[s * ldw], inc, &omega); // omega = <t, s>
        cblas_zdotc_sub(n, &work[t * ldw], inc, &work[t * ldw], inc, &temp); // temp = <t, t>
        omega /= temp; // omega = <t, s> / <t, t>

        cblas_zaxpy(n, &alpha, &work[phat * ldw], inc, &x[0], inc); // x = x + alpha*phat
        cblas_zaxpy(n, &omega, &work[shat * ldw], inc, &x[0], inc); // x = x + alpha*phat + omega*shat

        temp = -omega;
        cblas_zaxpy(n, &temp, &work[t * ldw], inc, &work[r * ldw], inc); // r = s - omega*t

        *resid = cblas_dznrm2(n, &work[r * ldw], inc) / bnrm2;

        if (*resid <= tol || *iter == maxit) {
            break;
        }

        if (rfo) {
            if (*iter % (maxit / 10) == 0) { // Check real error over reduced error
                // Compute ignored values
                cblas_zcopy(n, &b[0], inc, &work[red * ldw], inc);
                (*matvec)(&minus_one, &x[0], &one, &work[red * ldw]); // red = b - A*x
                cblas_zgbmv(CblasRowMajor, CblasNoTrans, n, n, 0, 0, &one, &work[neg * ldw], inc, &work[red * ldw], inc, &one, &x[0], inc); // x = x + (~rfo)(b - A*x)

                cblas_zcopy(n, &b[0], inc, &work[red * ldw], inc);
                (*matvec)(&minus_one, &x[0], &one, &work[red * ldw]); // red = b - A*x
                *resid = cblas_dznrm2(n, &work[red * ldw], inc) / cblas_dznrm2(n, &b[0], inc);

                if (*resid <= tol) {
                    return 0;
                }
            }
        }

        if (cabs(omega) < partol) {
            return -11;
        }

        rho1 = rho;
    }

    if (rfo) {
        // Compute ignored values
        cblas_zcopy(n, &b[0], inc, &work[r * ldw], inc);
        (*matvec)(&minus_one, &x[0], &one, &work[r * ldw]); // r = b - A*x;
        cblas_zgbmv(CblasRowMajor, CblasNoTrans, n, n, 0, 0, &one, &work[neg * ldw], inc, &work[r * ldw], inc, &one, &x[0], inc); // x = x + (~rfo)(b - A*x)
    }

    if (*iter == maxit) {
        return 1;
    }
    return 0;
}

void no_pre(double complex *restrict x, const double complex *restrict b) {
    cblas_zcopy(N, b, inc, x, inc);
}

void iprint(int info, int iter, double t) {
    switch(info) {
        case 0:
            printf("BiCGSTAB converged successfully in %.4f seconds after %d iterations.\n", t, iter);
            break;
        case -1:
            printf("Illegal parameter - n must be >= 0.\n");
            break;
        case -2:
            printf("Illegal parameter - ldw must be >= n.\n");
            break;
        case -3:
            printf("Illegal parameter - iter must be > 0.\n");
            break;
        case -10:
            printf("Breakdown after %d iteration(s) - rho and rtld have become orthogonal.\n", iter);
            break;
        case -11:
            printf("Breakdown after %d iteration(s) - s and t have become orthogonal relative to t'*t.\n", iter);
            break;
        default:
            printf("BiCGSTAB did not converge after %d iterations\n", iter);
            break;
    }
}
