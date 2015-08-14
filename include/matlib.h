/*
   Safer alloc() functions which check the returned pointer is valid.

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
#include <fftw3.h>
#include <stddef.h>

/**
 * Purpose
 * =======
 *
 * Safer version of calloc() which checks if the returned pointer was valid.
 *
 * Arguments
 * =========
 *
 * n        (input/output) ADDRESS of POINTER.
 *          The address of the pointer which should point to the new
 *          dynamically allocated section of memory.
 *
 * nelem    (input) SIZE_T.
 *          The number of elements to allocate.
 * 
 * elsize   (input) SIZE_T.
 *          The size of one element in memory to allocate.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Memory dynamically allocated and zeroed successfully.
 *
 *          = -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int calloc_s(void **ptr, size_t nelem, size_t elsize);

/**
 * Purpose
 * =======
 *
 * Safer version of malloc() which checks if the returned pointer was valid.
 *
 * Arguments
 * =========
 *
 * n        (input/output) ADDRESS of POINTER.
 *          The address of the pointer which should point to the new
 *          dynamically allocated section of memory.
 *
 * size     (input) SIZE_T.
 *          The size of memory to allocate.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Memory dynamically allocated and zeroed successfully.
 *
 *          = -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int malloc_s(void **ptr, size_t size);

/**
 * Purpose
 * =======
 *
 * Safer version of realloc() which checks if the returned pointer was valid.
 *
 * Arguments
 * =========
 *
 * n        (input/output) ADDRESS of POINTER.
 *          The address of the pointer which should point to the new
 *          dynamically allocated section of memory.
 *
 * size     (input) SIZE_T.
 *          The new size of memory to allocate.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Memory dynamically allocated and zeroed successfully.
 *
 *          = -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int realloc_s(void **ptr, size_t new_size);


/**
 * Purpose
 * =======
 *
 * Safer version of fftw_malloc() which checks if the returned pointer was
 * valid.
 *
 * Arguments
 * =========
 *
 * n        (input/output) ADDRESS of POINTER to FFTW_COMPLEX.
 *          The address of the pointer which should point to the new
 *          dynamically allocated section of memory.
 *          Note: The pointer should point to fftw_complex.
 *
 * size     (input) SIZE_T.
 *          The size of memory to allocate.
 *
 * Output
 * ======
 *
 *          INTEGER.
 *          =  0: Memory dynamically allocated and zeroed successfully.
 *
 *          = -1: Unable to allocate memory.
 *
 * =============================================================================
 */
int fftw_malloc_s(fftw_complex **ptr, size_t size);