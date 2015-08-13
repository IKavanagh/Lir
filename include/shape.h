/*
   Defines the shape and operations relating to it for which the electric field
   will be solved for.

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

#include <stdio.h>

#include "physics.h"

/**
 * Purpose
 * =======
 *
 * The length of each side in the x, y and z directions.
 *
 * =============================================================================
 */
extern double x, y, z;

/**
 * Purpose
 * =======
 *
 * The minimum and maximum values of the x, y and z sides.
 *
 * =============================================================================
 */
extern double xlim[2], ylim[2], zlim[2];

/**
 * Purpose
 * =======
 *
 * The number of discretisations in the x, y and z directions.
 *
 * =============================================================================
 */
extern size_t n, m, o;

/**
 * Purpose
 * =======
 *
 * A description of the shape in which the electric field will be solved for.
 *
 * =============================================================================
 */
extern int *shape;

/**
 * Purpose
 * =======
 *
 * The number of materials we have a description of in material.
 *
 * =============================================================================
 */
size_t materials;

/**
 * Purpose
 * =======
 *
 * An array of the different materials present in the shape.
 * Index 0 represents free space (air).
 *
 * =============================================================================
 */
extern material_t *material;

/**
 * Purpose
 * =======
 *
 * Reads in an entire shape description from the file filename and creates the
 * shape description based on the frequency of the problem, the number of
 * discretisations we want per wavelength and the minimum size of blocks we
 * want.
 *
 * The created shape is stored in the global variable shape.
 *
 * Arguments
 * =========
 *
 * filename (input) STRING.
 *          The name of the file the shape description is to be read from.
 *
 * f        (input) DOUBLE.
 *          The frequency we are solving the electric field for.
 *          Used to determine the number of discretised points needed.
 *
 * disc_per_lambda (input) INTEGER.
 *          The number of discretised points we require per wavelength.
 *
 * block_size (input) INTEGER.
 *          The minimum block size constraint we want to impose.
 *          The number of discretised points will be a multiple of this value.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Building created successfully.
 *
 *          =  1: File filename doesn't exist.
 *
 *          =  2: Error reading file.
 *
 *          =  -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int init_shape(const char *restrict filename, const double f, const int disc_per_lambda, const int block_size);

/**
 * Purpose
 * =======
 *
 * Creates the description of the shape based on the regions which have already
 * been determined, the frequency of the problem, the number of discretisations
 * we want per wavelength and the minimum size of blocks we want.
 * 
 * The number of discretised points along the x and y axis are determined based
 * on the frequency, f, the number of discretisations per lambda, 
 * disc_per_lambda, and the minimum block size, block_size.
 *
 * The created shape is stored in the global variable shape.
 *
 * Arguments
 * =========
 *
 * f        (input) DOUBLE.
 *          The frequency we are solving the electric field for.
 *          Used to determine the number of discretised points needed.
 *
 * disc_per_lambda (input) INTEGER.
 *          The number of discretised points we require per wavelength.
 *
 * block_size (input) INTEGER.
 *          The minimum block size constraint we want to impose.
 *          The number of discretised points will be a multiple of this value.
 * 
 * regions  (input) REGION_T array, dimension (region_count, 1).
 *          An array of regions. The shape is created using the dimensions of
 *          these regions and their material id.
 *
 * region_count (input) INTEGER.
 *          The number of regions in the regions array.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Shape created successfully.
 *
 *          =  -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int create_shape(const double f, const int disc_per_lambda, const int block_size, const region_t *restrict regions, const size_t region_count);

/**
 * Purpose
 * =======
 *
 * Reads in one region from the file pointed to by fp and stores it in region.
 * 
 * If the region is made up of a new material it stores this in the material
 * array and increases the number of materials known.
 *
 * Arguments
 * =========
 *
 * fp       (input) FILE.
 *          A pointer to the current location in the file where the region is to
 *          be read from.
 *
 * region   (output) REGION_T.
 *          Stores the region that was read in from the file.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Building created successfully.
 *
 *          =  2: Error reading file.
 *
 *          =  -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int read_region(FILE *restrict fp, region_t *restrict region);

/**
 * Purpose
 * =======
 *
 * Reads in an entire shape description from file and stores each region in the
 * regions array and the number of regions read in into region_count.
 *
 * Arguments
 * =========
 *
 * filename (input) STRING.
 *          The name of the file the shape description is to be read from.
 *
 * regions  (output) REGION_T array.
 *          Stores all of the regions that were read in from the file.
 *
 * region_count (output) INTEGER.
 *          The number of regions that were read in from the file and
 *          consequently are stored in regions.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Building created successfully.
 *
 *          =  1: File filename doesn't exist.
 *
 *          =  2: Error reading file.
 *
 *          =  -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int read_shape(const char *restrict filename, region_t **restrict regions, size_t *restrict region_count);

/**
 * Purpose
 * =======
 *
 * Cleanup any dynamically allocated memory.
 *
 * =============================================================================
 */
void shape_cleanup(void);
