#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MKL_Complex16
#define MKL_Complex16 double complex
#endif

#include <mkl.h>

#include "krylov.h"
#include "matio.h"
#include "monotonic_clock.h"
#include "shape.h"
#include "vefie.h"

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Error: Wrong number of arguments!\n");
        fprintf(stderr, "usage: %s f a_x a_y\n", argv[0]);
        exit(4 - argc);
    }

    int err = 0;
    double t = monotonic_clock();

    const char *filename = "input/building.txt";
    double f = atof(argv[1]);
    double complex antenna = atof(argv[2]) + I*atof(argv[3]);

    if (argc > 4) {
        filename = argv[4];
    }

    if ((err = init_shape(filename, f, 10, block_size)) ||
        (err = init_vefie(f, antenna, hertzian_dipole))) {
        exit(err);
    }

    t = monotonic_clock() - t;

    printf("Initialised VEFIE elements for a problem size of %d in %.4f seconds.\n", n*m, t);

    t = monotonic_clock();

    int iter = 1000;
    double resid = 1e-3;
    int N = n*m;

    double complex *work = mkl_malloc((size_t) N * 8 * sizeof *work, alignment);
    int info = rbicgstab(n*m, V, E, work, n*m, &iter, &resid, matvec, no_pre, rfo, 1);

    t = monotonic_clock() - t;
    
    iprint(info, iter, t);
    if (info == 0) {
        FILE *fp;
        int dimensions[2] = {1};

        fp = fopen("output/E.txt", "w");
        dimensions[0] = n*m;
        dimensions[1] = 1;
        zmprint(fp, 1, dimensions, E);
        fclose(fp);

        fp = fopen("output/position.txt", "w");
        dimensions[0] = n;
        dimensions[1] = m;
        zmprint(fp, 2, dimensions, p);
        fclose(fp);

        fp = fopen("output/shape.txt", "w");
        dimensions[0] = n;
        dimensions[1] = m;
        imprint(fp, 2, dimensions, shape);
        fclose(fp);

        fp = fopen("output/material.txt", "w");
        dimensions[0] = materials;
        dimensions[1] = 1;
        mmprint(fp, 2, dimensions, material);
        fclose(fp);

        fp = fopen("output/D.txt", "w");
        dimensions[0] = n*m;
        dimensions[1] = 1;
        zmprint(fp, 1, dimensions, D);
        fclose(fp);

        fp = fopen("output/G.txt", "w");
        dimensions[0] = n*m;
        dimensions[1] = 1;
        zmprint(fp, 1, dimensions, G);
        fclose(fp);

        fp = fopen("output/V.txt", "w");
        dimensions[0] = n*m;
        dimensions[1] = 1;
        zmprint(fp, 1, dimensions, V);
        fclose(fp);

        fp = fopen("output/rfo.txt", "w");
        dimensions[0] = n*m;
        dimensions[1] = 1;
        zmprint(fp, 1, dimensions, rfo);
        fclose(fp);
    }

    shape_cleanup();
    vefie_cleanup();

    return 0;
}