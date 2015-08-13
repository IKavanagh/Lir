/*
   Defines the shape and operations relating to it for which the VEFIE
   is applied to.

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

#include "shape.h"

#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>

#include "matlib.h"

double x = 0.0, y = 0.0, z = 0.0;
int n = 0, m = 0, o = 0;

double xlim[2] = {DBL_MAX, -DBL_MAX};
double ylim[2] = {DBL_MAX, -DBL_MAX};
double zlim[2] = {DBL_MAX, -DBL_MAX};

int materials = 0;

int *shape = NULL;
material_t *material = NULL;

int init_shape(const char *restrict filename, const double f, const int disc_per_lambda, const int block_size) {
    // TODO: Add error checking
    int region_count;
    region_t *regions;

    int ret_code = read_shape(filename, &regions, &region_count);
    if (ret_code == 0) {
        ret_code = create_shape(f, disc_per_lambda, block_size, regions, region_count);
    }

    free(regions);

    return ret_code;
}

int create_shape(const double f, const int disc_per_lambda, const int block_size, const region_t *restrict regions, const int region_count) {
    double lambda0 = lambda(kd(f, 1, 1, 0));

    n = (int) floor (x / (lambda0 / disc_per_lambda));
    n += block_size - n % block_size;

    double dx = x / (double) n;

    m = (int) floor (y / (lambda0 / disc_per_lambda));
    m += block_size - m % block_size;

    double dy = y / (double) m;

    o = (int) floor (z / (lambda0 / disc_per_lambda));
    o += block_size - o % block_size;

    double dz = z / (double) o;

    if (!shape) {
        free(shape); // Guard against memory leaks
    }
    int ret_code = calloc_s((void **) &shape, (size_t) (o*m*n), sizeof *shape);
    if (ret_code != 0) {
        return ret_code;
    }

    for (int l = 0; l < region_count; ++l) {
        region_t region = regions[l];

        int x_range[2] = {INT_MAX}, y_range[2] = {INT_MAX}, z_range[2] = {INT_MAX};

        x_range[1] = 0;
        y_range[1] = 0;
        z_range[1] = 0;

        for (int i = 0; i < 4; ++i) {
            int Nx = (int) round ((region.vertex[i].x - xlim[0]) / dx);
            int My = (int) round ((region.vertex[i].y - ylim[0]) / dy);

            x_range[0] = min(x_range[0], Nx);
            x_range[1] = max(x_range[1], Nx);

            y_range[0] = min(y_range[0], My);
            y_range[1] = max(y_range[1], My);

            // z must be treated differently as its given in a base coordinate
            // and a height
            int Oz = (int) round ((region.vertex[i].z - zlim[0]) / dz);
            z_range[0] = min(z_range[0], Oz);

            Oz = (int) round ((region.height - zlim[0]) / dz);
            z_range[1] = max(z_range[1], Oz);
        }

        for (int k = z_range[0]; k < z_range[1]; ++k) {
            for (int j = y_range[0]; j < y_range[1]; ++j) {
                for (int i = x_range[0]; i < x_range[1]; ++i) {
                    int idx = (k * m + j) * n + i;
                    shape[idx] = region.material_id;
                }
            }
        }
    }

    return 0;
}

int read_shape(const char *restrict filename, region_t **restrict regions, int *restrict region_count) {
    int region_size = 4, count = 0;

    region_t *region;
    int ret_code = malloc_s((void **) &region, (size_t) region_size * sizeof *region);
    if (ret_code != 0) {
        return ret_code;
    }

    if (materials == 0) {
        if (!material) {
            free(material);
        }
        ret_code = malloc_s((void **) &material, sizeof *material);
        if (ret_code != 0) {
            return ret_code;
        }

        material[0] = free_space(); // Free space
        materials = 1;
    }

    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "%s: File or directory does not exist!", filename);
        return 1;
    }

    int c; // int required to handle EOF
    while ((c = fgetc(fp)) != EOF) {
        if (c == '#') { // Comment skip to next line
            while ((c = fgetc(fp)) != EOF && c != '\n');
            c = fgetc(fp); // First character of next line
        }

        if (c == '(') { // Beginning of coordinate list
            ret_code = ungetc(c, fp);
            if (ret_code == EOF) {
                fprintf(stderr, "ungetc(): Failed to put %c back into %s at line %d\n", c, __FILE__, __LINE__-6);
                return 2;
            }

            ret_code = read_region(fp, &region[count]);
            if (ret_code != 0) {
                return ret_code;
            }

            count++;
            if (count == region_size) {
                region_size *= 2;

                ret_code = realloc_s((void **) &region, (size_t) region_size * sizeof *region);
                if(ret_code != 0) {
                    return ret_code;
                }
            }
        }
    }

    x = xlim[1] - ylim[0];
    y = ylim[1] - ylim[0];
    z = zlim[1] - zlim[0];

    region_t *r = *regions;
    if (!r) {
        free(r);
    }
    *regions = region;
    *region_count = count;

    if (ferror(fp)) {
        fprintf(stderr, "fgetc(): Failed in %s at #%d\n", __FILE__, __LINE__-6);
        return 2;
    } else if (feof(fp)) {
        fclose(fp);
    }

    return 0;
}

int read_region(FILE *restrict fp, region_t *restrict region) {
    double height = 0;
    point_t p;
    for (int i = 0; i < 4; ++i) {
        if(fscanf(fp, "(%lf, %lf, %lf) ", &p.x, &p.y, &p.z) == EOF) {
            fprintf(stderr, "fscanf(): Failed to read coordinates from %s at line %d\n", __FILE__, __LINE__-6);
            return 2;
        }

        (*region).vertex[i] = p;

        xlim[0] = fmin(xlim[0], p.x);
        xlim[1] = fmax(xlim[1], p.x);

        ylim[0] = fmin(ylim[0], p.y);
        ylim[1] = fmax(ylim[1], p.y);
    }
    if(fscanf(fp, "%lf ", &height) == EOF) {
        fprintf(stderr, "fscanf(): Failed to read height from %s at line %d\n", __FILE__, __LINE__-6);
        return 2;
    }
    (*region).height = height;

    zlim[0] = fmin(zlim[0], p.z);
    zlim[1] = fmax(zlim[1], p.z + height);

    material_t mat;
    if(fscanf(fp, "%lf %lf %lf", &mat.epsilon_r, &mat.mu_r, &mat.sigma) == EOF) {
        fprintf(stderr, "fscanf(): Failed to read material from %s at line %d\n", __FILE__, __LINE__-6);
        return 2;
    }

    int material_id = -1;
    for (int i = 0; i < materials && material_id == -1; ++i) {
        if (matcmp(mat, material[i])) {
            material_id = i;
        }
    }

    if (material_id == -1) { // New material
        materials++;

        int ret_code = realloc_s((void **) &material, (size_t) materials * sizeof *material);
        if (ret_code != 0) {
            return ret_code;
        }

        material_id = materials - 1;
        material[materials - 1] = mat;
    }

    (*region).material_id = material_id;

    return 0;
}

void shape_cleanup(void) {
    if (!shape) free(shape);
    if (!material) free(material);
}
