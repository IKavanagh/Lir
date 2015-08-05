/*
   Helper functions to compute various mathematical and physical constants.

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
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * Purpose
 * =======
 *
 * Represents a material. Contains the relative permittivity, epsilon_r, and
 * permeability, mu_r, of the material as well as its conductivity sigma.
 *
 * =============================================================================
 */
typedef struct {
    double epsilon_r;
    double mu_r;
    double sigma;
} material_t;

/**
 * Purpose
 * =======
 *
 * Represents a point in 3D space. Contains an x, y and z coordinate.
 *
 * =============================================================================
 */
typedef struct {
    double x;
    double y;
    double z;
} point_t;

/**
 * Purpose
 * =======
 *
 * Represents a region in space bounded by 4 points with a vertical height and a
 * specific material assigned to it. Each vertex is a 3D point.
 *
 * The material is stored as an id which is used to reference a material type.
 *
 * =============================================================================
 */
typedef struct {
    point_t vertex[4];
    double height;
    int material_id;
} region_t;

/**
 * Purpose
 * =======
 *
 * Determine which of its argument is the smallest value.
 *
 * Arguments
 * =========
 *
 * a        (input) INTEGER.
 *          First argument to compare with.
 *
 * b        (input) DOUBLE.
 *          Second argument to compare with.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          The smallest of its two arguments.
 *
 * =============================================================================
 */
int min(int a, int b);

/**
 * Purpose
 * =======
 *
 * Determine which of its argument is the largest value.
 *
 * Arguments
 * =========
 *
 * a        (input) INTEGER.
 *          First argument to compare with.
 *
 * b        (input) DOUBLE.
 *          Second argument to compare with.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          The largest of its two arguments.
 *
 * =============================================================================
 */
int max(int a, int b);

/**
 * Purpose
 * =======
 *
 * Checks for equality between two material structs.
 *
 * Arguments
 * =========
 *
 * a        (input) MATERIAL.
 *          First material to compare.
 *
 * b        (input) MATERIAL.
 *          Second material to compare.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          = 0: Materials match.
 *
 *          = otherwise: Materials do not match.          
 *
 * =============================================================================
 */
int matcmp(material_t m1, material_t m2);

/**
 * Purpose
 * =======
 *
 * Defines the material properties of free space.
 *
 * Output
 * ======
 *
 *          MATERIAL_T.
 *          A material_t struct containing the appropriate values for
 *          epsilon_r, mu_r and sigma to denote a free space material.          
 *
 * =============================================================================
 */
material_t free_space(void);

/**
 * Purpose
 * =======
 *
 * Computes the Hankel function of the first kind.
 *
 * Arguments
 * =========
 * 
 * n        (input) INTEGER.
 *          Order of the Hankel function.
 * 
 * x        (input) DOUBLE.
 *          Value to compute the Hankel function for.
 * 
 * Output
 * ======
 *
 *          COMPLEX DOUBLE.
 *          The computed Hankel function.
 *
 * =============================================================================
 */
double complex h1n(int n, double x);

/**
 * Purpose
 * =======
 *
 * Computes the Hankel function of the second kind.
 *
 * Arguments
 * =========
 * 
 * n        (input) INTEGER.
 *          Order of the Hankel function.
 * 
 * x        (input) DOUBLE.
 *          Value to compute the Hankel function for.
 * 
 * Output
 * ======
 *
 *          COMPLEX DOUBLE.
 *          The computed Hankel function.
 *
 * =============================================================================
 */
double complex h2n(int n, double x);

/**
 * Purpose
 * =======
 *
 * Computes the wave number associated with an electromagnetic wave passing
 * through a medium.
 *
 * Arguments
 * =========
 * 
 * f        (input) DOUBLE.
 *          Frequency of the electromagnetic wave.
 * 
 * epsilonr (input) DOUBLE.
 *          Relative permittivity of the medium.
 * 
 * mur      (input) DOUBLE.
 *          Relative permeability of the medium.
 * 
 * sigma    (input) DOUBLE.
 *          Conductivity of the medium.
 * 
 * Output
 * ======
 *
 *          COMPLEX DOUBLE.
 *          Wave number of an electromagnetic wave passing through a defined
 *          medium.
 *
 * =============================================================================
 */
double complex k(double f, double epsilonr, double mur, double sigma);

/**
 * Purpose
 * =======
 *
 * Computes the wavelength of an electromagnetic wave passing through a medium.
 *
 * Arguments
 * =========
 * 
 * k        (input) DOUBLE COMPLEX.
 *          Wave number of the electromagnetic wave.
 * 
 * Output
 * ======
 *
 *          DOUBLE.
 *          Wavelength of the electromagnetic wave.
 *
 * =============================================================================
 */
double lambda(double complex k);

/**
 * Purpose
 * =======
 *
 * Computes the value of a hertzian dipole source at the point p relative to
 * the antenna. In 2D the Hertzian dipole is represented as a line source given
 * by
 *              H_{0}^{2}(k0*R)
 * where H_{0}^{(2)} is the zeroth order bessel function of the second kind.
 *
 * Arguments
 * =========
 * 
 * k0       (input) DOUBLE.
 *          Free space wave number.
 *
 * antenna  (input) DOUBLE COMPLEX.
 *          Position of the antenna in the problem.
 *          The real part denotes the x-coordinate and the imaginary part the
 *          y-coordinate.
 * 
 * point    (input) DOUBLE COMPLEX.
 *          x-y position in the problem the incident field is to be computed at.
 *          The real part denotes the x-coordinate and the imaginary part the
 *          y-coordinate.
 * 
 * Output
 * ======
 *
 *          DOUBLE COMPLEX.
 *          Value of a hertzian dipole source at the point p relative to the
 *          antenna.
 *
 * =============================================================================
 */
double complex hertzian_dipole(const double k0, const double complex antenna, const double complex point);

/**
 * Purpose
 * =======
 *
 * Computes the value of a plane value given by
 *              exp(-jkp)
 *
 * Arguments
 * =========
 * 
 * k0       (input) DOUBLE.
 *          Free space wave number.
 *
 * antenna  (input) DOUBLE COMPLEX.
 *          Position of the antenna in the problem.
 *          The real part denotes the x-coordinate and the imaginary part the
 *          y-coordinate.
 * 
 * point    (input) DOUBLE COMPLEX.
 *          x-y position in the problem the incident field is to be computed at.
 *          The real part denotes the x-coordinate and the imaginary part the
 *          y-coordinate.
 * 
 * Output
 * ======
 *
 *          DOUBLE COMPLEX.
 *          Value of plane wave source at the point p relative to the antenna.
 *
 * =============================================================================
 */
double complex plane_wave(const double k0, const double complex antenna, const double complex point);
