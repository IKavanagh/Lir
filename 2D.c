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
        fprintf(stderr, "usage: %s f x y [file]\n", argv[0]);
        exit(4 - argc);
    }

    int err = 0;

    const char *filename = "input/building.txt";
    double f = atof(argv[1]);
    double complex antenna = atof(argv[2]) + I*atof(argv[3]);

    if (argc > 4) {
        filename = argv[4];
    }

    double t = monotonic_clock();

    if ((err = init_shape(filename, f, 10, block_size)) ||
        (err = init_vefie(f, antenna, hertzian_dipole))) {
        exit(err);
    }

    t = monotonic_clock() - t;

    printf("Initialised VEFIE elements for a problem size of %d in %.4f seconds.\n", n*m, t);

    int iter = 1000, N = n*m, info;
    double resid = 1e-3;

    double complex *work = mkl_malloc((size_t) N * 9 * sizeof *work, alignment);

    t = monotonic_clock();
    info = rbicgstab(N, V, E, work, N, &iter, &resid, matvec, no_pre, rfo);
    t = monotonic_clock() - t;

    iprint(info, iter, t);
    if (info == 0) {
        FILE *fp;
        size_t dimensions[2] = {1};

        fp = fopen("output/E.txt", "w");
        dimensions[0] = (size_t) N;
        dimensions[1] = 1;
        mprint(fp, 2, dimensions, E, sizeof *E, printz);
        fclose(fp);

        fp = fopen("output/position.txt", "w");
        dimensions[0] = (size_t) n;
        dimensions[1] = (size_t) m;
        mprint(fp, 2, dimensions, p, sizeof *p, printz);
        fclose(fp);

        fp = fopen("output/shape.txt", "w");
        dimensions[0] = (size_t) n;
        dimensions[1] = (size_t) m;
        mprint(fp, 2, dimensions, shape, sizeof *shape, printi);
        fclose(fp);

        fp = fopen("output/material.txt", "w");
        dimensions[0] = (size_t) materials;
        dimensions[1] = 1;
        mprint(fp, 2, dimensions, material, sizeof *material, printm);
        fclose(fp);

        fp = fopen("output/D.txt", "w");
        dimensions[0] = (size_t) N;
        dimensions[1] = 1;
        mprint(fp, 2, dimensions, D, sizeof *D, printz);
        fclose(fp);

        fp = fopen("output/G.txt", "w");
        dimensions[0] = (size_t) N;
        dimensions[1] = 1;
        mprint(fp, 2, dimensions, G, sizeof *G, printz);
        fclose(fp);

        fp = fopen("output/V.txt", "w");
        dimensions[0] = (size_t) N;
        dimensions[1] = 1;
        mprint(fp, 2, dimensions, V, sizeof *V, printz);
        fclose(fp);

        fp = fopen("output/rfo.txt", "w");
        dimensions[0] = (size_t) N;
        dimensions[1] = 1;
        mprint(fp, 2, dimensions, rfo, sizeof *rfo, printz);
        fclose(fp);
    }

    shape_cleanup();
    vefie_cleanup();

    return 0;
}