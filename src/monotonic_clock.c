/*
   A monotonic clock which increases by 1 second for every real second, can be
   used to time program execution.

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

   Adapted from:
   ThomasHabets/monotonic_clock <thomas@habets.se> 2010
   https://github.com/ThomasHabets/monotonic_clock
*/

#include "monotonic_clock.h"

#if defined (__posix__)
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200112L
#endif
#elif defined(__MACH__)
#include <mach/mach_time.h>
#elif defined (__WIN32__)
#include <Windows.h>
#endif

#include <sys/time.h>
#include <time.h>

/**
 * Purpose
 * =======
 *
 * Fall back clock using either gettimeofday() or time() to
 * return a value which can be used to time program execution.
 * 
 * Output
 * ======
 *
 *          DOUBLE.
 *          A number of seconds.
 *
 * ============================================================
 */
double monotonic_clock_fallback(void) {
#ifdef __linux__
    struct timeval tv;
    if (!gettimeofday(&tv, NULL)) {
        return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
    }
#endif /* __linux__ */
    return (double) time(NULL);
}

double monotonic_clock(void) {
  // TODO: Test Windows version
#if defined (__posix__)
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return ts.tv_sec + ts.tv_nsec / 1000000000.0;
    }
    
    return monotonic_clock_fallback();
#elif defined(__MACH__)
    uint64_t t = mach_absolute_time();
    static double scaling_factor = 0;

    if (scaling_factor < 1e-15) {
        mach_timebase_info_data_t info;
        if (mach_timebase_info(&info) != 0) {
            return monotonic_clock_fallback();
        }
        scaling_factor = (double) info.numer / (double) info.denom;
    }

    return (double) t * scaling_factor / 1000000000.0;
#elif defined (__WIN32__)
    static uint64_t scaling_factor = 0;

    LARGE_INTEGER count;
    if (QueryPerformanceCounter(&count) == 0) {
        return monotonic_clock_fallback();
    }

    if (scaling_factor == 0) {
        LARGE_INTEGER frequency;
        if (QueryPerformanceFrequency(&frequency) == 0) {
            return monotonic_clock_fallback();
        }
        scaling_factor = frequency.QuadPart * 1000000000;
    }

    return (double) count.QuadPart / scaling_factor;
#endif
}
