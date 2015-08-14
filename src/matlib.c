#include "matlib.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int alignment = 64;

int aerror(void *ptr, size_t size) {
    if (!ptr) {
        fprintf(stderr, "Failed: Unable to allocate %lu bytes.\n", size);
        return -1;
    }
    return 0;
}

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

int calloc_align_s(void **ptr, size_t nelem, size_t elsize) {
    size_t size = nelem*elsize;

    void *X;
    
    int ans = posix_memalign(&X, alignment, size);
    if (ans != 0) {
        return aerror(NULL, size);
    }

    memset(X, 0, size);

    *ptr = X;
    return aerror(*ptr, size);
}

int malloc_align_s(void **ptr, size_t size) {
    int ans = posix_memalign(ptr, alignment, size);
    if (ans != 0) {
        return aerror(NULL, size);
    }
    return aerror(*ptr, size);
}
