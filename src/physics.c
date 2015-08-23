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

inline material_t free_space(void) {
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

    return -I * sqrt(I * omega * mu * (sigma + I * omega * epsilon));
}

inline double lambda(double complex k) {
    return 2.0 * M_PI / creal(k);
}

void hertzian_dipole_2d(const double k0, const double antenna[2], const double point[2], double complex field[3]) {
    double R = sqrt(pow(antenna[0] - point[0], 2) + pow(antenna[1] - point[1], 2));

    field[0] = 0.0;
    field[1] = 0.0;
    field[2] = h2n(0, k0*R);
}

void plane_wave_2d(const double k0, const double antenna[2], const double point[2], double complex field[3]) {
    field[0] = 0.0;
    field[1] = 0.0;
    field[2] = exp(-I*k0*(point[0] - antenna[0]));
}

void hertzian_dipole_3d(const double k0, const double f, const double antenna[3], const double point[3], double complex field[3]) {
    double mu0 = 4e-7*M_PI;
    double omega = 2*M_PI*f;

    double l = lambda(k0) / 2;
    double i = 0.01;

    double x = point[0] - antenna[0];
    double y = point[1] - antenna[1];
    double z = point[2] - antenna[2];

    double R = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    double phi = atan2(y, x);
    double theta = atan2(sqrt(pow(x, 2) + pow(y, 2)), z);

    double complex theta_field = (1 / R) * (((i*l) / (4 * M_PI)) * I * omega * mu0);

    theta_field *= exp(-I*k0*R) * sin(theta);

    field[0] = theta_field * cos(theta) * cos(phi);
    field[1] = theta_field * cos(theta) * sin(phi);
    field[2] = theta_field * sin(theta);
}

void plane_wave_3d(const double k0, const double f, const double antenna[3], const double point[3], double complex field[3]) {
    field[0] = 0.0;
    field[1] = 0.0;
    field[2] = exp(-I*k0*(point[0] - antenna[0]));
}
