/*
   Helper functions to compute various mathematical and physical
   constants.

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

#include "physics.h"

#define _GNU_SOURCE // Required to enable special math functions (jn and yn)

inline int min(int a, int b) {
    if (b < a) return b;
    return a;
}

inline int max(int a, int b) {
    if (b > a) return b;
    return a;
}

inline int matcmp(material_t m1, material_t m2) {
    return fabs(m1.epsilon_r - m2.epsilon_r) < 1e-15 && fabs(m1.mu_r - m2.mu_r) < 1e-15 && fabs(m1.sigma - m2.sigma) < 1e-15;
}

material_t free_space(void) {
    return (material_t) { .epsilon_r = 1.0, .mu_r = 1.0, .sigma = 0.0 };
}

inline double complex h1n(int n, double x) {
    return jn(n, x) + yn(n, x)*I;
}

inline double complex h2n(int n, double x) {
    return jn(n, x) - yn(n, x)*I;
}

double complex kd(double f, double epsilonr, double mur, double sigma) {
    double omega = 2.0 * M_PI * f;
    double mu = mur * (4.0 * M_PI * 1e-7);
    double epsilon = epsilonr * 8.854e-12;

    return -I * csqrt(I * omega * mu * (sigma + I * omega * epsilon));
}

inline double lambda(double complex k) {
    return 2.0 * M_PI / creal(k);
}

inline double complex hertzian_dipole(const double k0, const double complex antenna, const double complex point) {
    return h2n(0, k0*cabs(point - antenna));
}

inline double complex plane_wave(const double k0, const double complex antenna, const double complex point) {
    return cexp(-I*k0*creal(point));
}
