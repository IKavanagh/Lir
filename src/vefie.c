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

#include "vefie.h"

#include <mkl.h>
#include <string.h>

#if defined(_OPENMP)
#include <omp.h>
#else
void omp_set_num_threads(int num_threads) { return; }
int omp_get_num_threads(void) { return 1; }
int omp_get_max_threads(void) { return 1; }
#endif

#include "matlib.h"
#include "physics.h"
#include "shape.h"

/**
 * Purpose
 * =======
 *
 * The dimensions of the FFT used to perform the matrix vector
 * product of A*x.
 *
 * ============================================================
 */
int dims[3];

/**
 * Purpose
 * =======
 *
 * Input and output arrays for the FFT.
 *
 * ============================================================
 */
double complex *in, *out;

/**
 * Purpose
 * =======
 *
 * The plans for performing the forward and backward FFT when
 * computing matrix vector product of A*x.
 *
 * ============================================================
 */
fftw_plan forward_plan, backward_plan;

double complex *D, *G, *V;
double complex *E;
double complex *p;
double complex *rfo;

int block_size = 20;

int init_vefie(const double f, const double complex antenna, double complex (*inc)(const double k0, const double complex antenna, const double complex point)) {
    // TODO: Document possible return values
    int ret_code;
    if (omp_get_max_threads() != 1 && fftw_init_threads()) {
        omp_set_num_threads(omp_get_max_threads());
        fftw_plan_with_nthreads(omp_get_max_threads());
    } else {
        fprintf(stderr, "Warning: fftw_init_threads() failed!\n");
    }

    if (f < 0) {
        return -10;
    }
    
    ret_code = allocate_matrices();
    if (ret_code != 0) {
        return ret_code;
    }
    
    ret_code = init_fftw();
    if (ret_code != 0) {
        return ret_code;
    }

    const char *unit = "Hz";

    if (3 * floor(log10(f) / 3) - 6 < 1e-15) {
        unit = "MHz";
    } else if (3 * floor(log10(f) / 3) - 9 < 1e-15) {
        unit = "GHz";
    }

    #pragma omp parallel
    #pragma omp master
    printf("Running VEFIE with %d threads for a shape of size %.2fm x %.2fm with a source radiating at %.2f%s positioned at (%.2f,%.2f).\n", omp_get_num_threads(), x, y, f / (strcmp(unit, "MHz") == 0 ? 1e6 : strcmp(unit, "GHz") == 0 ? 1e9 : 1), unit, creal(antenna), cimag(antenna));

    double dx = x / (double) n, dy = y / (double) m, a = sqrt(dx * dy / M_PI);

    double k0 = creal(kd(f, 1, 1, 0));

    #pragma omp parallel for
    for (int k = 0; k < n*m; ++k) {
        int i = k % m;
        int j = k / m;

        // Define position as centre of each cube
        p[k] = (xlim[0] + I*ylim[0]) + (((double) i + 0.5)*dx + I*((double) j + 0.5)*dy);

        V[k] = (*inc)(k0, antenna, p[k]);

        material_t mat = material[shape[j + m * i]];
        double complex km = kd(f, mat.epsilon_r, mat.mu_r, mat.sigma);

        D[k] = km*km - k0*k0;
        if (cabs(D[k]) > 1e-15) {
            rfo[k] = 1.0;
        }

        double Rmn = cabs(p[0] - p[k]);
        G[k] = (I / 4) * (((2*M_PI*a) / k0) * j1(k0*a) * h2n(0, k0*Rmn));
    }
    G[0] = (I / 4) * (((2*M_PI*a) / k0) * h2n(1, k0*a) - 4*I / (k0*k0));

    return toeplitz(&G); // Extend G to size 2N * 2M and take FFT
}

int init_fftw(void) {
    dims[0] = 2*n;
    dims[1] = 2*m;

    int N = dims[0] * dims[1];

    int ret_code = fftw_malloc_s(&in, (size_t) N * sizeof *in);
    if (ret_code != 0) {
        return ret_code;
    }

    ret_code = fftw_malloc_s(&out, (size_t) N * sizeof *out);
    if (ret_code != 0) {
        return ret_code;
    }

    forward_plan = fftw_plan_dft(2, dims, in, out, FFTW_FORWARD, FFTW_MEASURE);
    backward_plan = fftw_plan_dft(2, dims, in, out, FFTW_BACKWARD, FFTW_MEASURE);

    if (!forward_plan || !backward_plan) {
        fprintf(stderr, "Failed: Unable to plan FFT.\n");
        return -2;
    }

    return 0;
}

int toeplitz(double complex** restrict X) {
    // TODO: Fix when compiling with gcc
    //       * Produces a different result (appears to be less accurate)
    //       * Could depend on which library is being used for the FFT
    //         (Intel or FFTW3)
    //       * Other likely option, use of #pragma omp parallel for varies
    //       * between compiler
    double complex* Y = *X;

    double complex *Z;
    int ret_code = fftw_malloc_s(&Z, (size_t) (dims[0]*dims[1]) * sizeof *Z);
    if (ret_code != 0) {
        // TODO: Document possible return code
        return ret_code;
    }

    fftw_plan plan = fftw_plan_dft(2, dims, Z, Z, FFTW_FORWARD, FFTW_ESTIMATE);

    // Embed Y into a circular convolution problem of size 2Nx2M
    #pragma omp parallel for
    for (int k = 0; k < n*m; ++k) {
        int i = k / m; // Row
        int j = k % m; // Column

        Z[i * dims[1] + j] = Y[i * m + j]; // Top left
        if (j == 0) {
            Z[i * dims[1] + (j + m)] = 0.0; // Top right
        } else {
            Z[i * dims[1] + (j + m)] = Y[(i + 1) * m - j]; // Top right
        }

        if (i == 0) {
            Z[(i + m) * dims[1] + j] = 0.0; // Bottom left
            Z[(i + m) * dims[1] + (j + m)] = 0.0; // Bottom right
        } else {
            Z[(i + m) * dims[1] + j] = Y[(m - i) * m + j]; // Bottom left
            if (j == 0) {
                Z[(i + m) * dims[1] + (j + m)] = 0.0; // Bottom right
            } else {
                Z[(i + m) * dims[1] + (j + m)] = Y[(m + 1 - i) * m - j]; // Bottom right
            }
        }
    }

    fftw_execute(plan);

    fftw_free(Y);
    *X = Z;

    fftw_destroy_plan(plan);
    fftw_cleanup();

    return 0;
}

void matvec(const double complex *restrict alpha, const double complex *restrict X, const double complex *restrict beta, double complex *restrict Y) {
    // TODO: Fix when compiling with gcc
    //       * Produces a different result (appears to be less accurate)
    //       * Could depend on which library is being used for the FFT
    //         (Intel or FFTW3)
    //       * Other likely option, use of #pragma omp parallel for varies
    //       * between compiler
    #pragma omp parallel for
    for (int k = 0; k < n*m; ++k) {
        int i = k / m;
        int j = k % m;

        in[i * dims[1] + j] = D[i * m + j] * X[i * m + j];
        in[i * dims[1] + (j + m)] = 0.0;
        in[(i + m) * dims[1] + j] = 0.0;
        in[(i + m) * dims[1] + (j + m)] = 0.0;
    }

    fftw_execute(forward_plan);

    #pragma omp parallel for
    for (int k = 0; k < dims[0] * dims[1]; ++k) {
        in[k] = G[k] * out[k];
    }

    fftw_execute(backward_plan);

    #pragma omp parallel for
    for (int k = 0; k < n*m; ++k) {
        int i = k / m;
        int j = k % m;

        Y[k] = *alpha * (X[k] + (out[i * dims[1] + j] / (dims[0]*dims[1]))) + *beta*Y[k];
    }
}

int allocate_matrices(void) {
    size_t N = (size_t) (n*m);

    int ret_code = matalloc(N, 0, &D);
    if (ret_code != 0) {
        return ret_code;
    }

    ret_code = fftw_malloc_s(&G, N * sizeof *G);
    if (ret_code != 0) {
        return ret_code;
    }

    ret_code = matalloc(N, 0, &V);
    if (ret_code != 0) {
        return ret_code;
    }

    ret_code = matalloc(N, 0, &E);
    if (ret_code != 0) {
        return ret_code;
    }

    ret_code = matalloc(N, 0, &p);
    if (ret_code != 0) {
        return ret_code;
    }

    ret_code = matalloc(N, 0, &rfo);
    if (ret_code != 0) {
        return ret_code;
    }

    return 0;
}

int matalloc(const size_t N, const int init, double complex **restrict X) {
    int ret_code;
    if (init) {
        ret_code = calloc_align_s((void **) X, N, sizeof **X);
    } else {
        ret_code = malloc_align_s((void **) X, N * sizeof **X);
    }
    return ret_code;
}

void vefie_cleanup(void) {
    if (!D) free(D);
    if (!G) fftw_free(G);
    if (!V) free(V);

    if (!E) free(E);

    if (!p) free(p);

    if (!rfo) free(rfo);

    if (!in) fftw_free(in);
    if (!out) fftw_free(out);

    if (!forward_plan) fftw_destroy_plan(forward_plan);
    if (!backward_plan) fftw_destroy_plan(backward_plan);

    fftw_cleanup_threads();
    fftw_cleanup();
}
