#include "matlib.h"

#include <stdlib.h>
#include <stdio.h>

int aerror(void *ptr, size_t size);

int calloc_s(void **ptr, size_t nelem, size_t elsize) {
    *ptr = calloc(nelem, elsize);
    return aerror(*ptr, nelem * elsize);
}

int malloc_s(void **ptr, size_t size) {
    *ptr = malloc(size);
    return aerror(*ptr, size);
}

int realloc_s(void **ptr, size_t size) {
    *ptr = realloc(*ptr, size);
    return aerror(*ptr, size);
}

int fftw_malloc_s(fftw_complex **ptr, size_t size) {
    *ptr = fftw_malloc(size);
    return aerror((void *) *ptr, size);
}

int aerror(void *ptr, size_t size) {
    if (!ptr) {
        fprintf(stderr, "Failed: Unable to allocate %lu bytes.\n", size);
        return -1;
    }
    return 0;
}