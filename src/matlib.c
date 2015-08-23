#include "matlib.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

size_t alignment = 64;

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
    void *new_ptr = realloc(*ptr, size);
    if (!new_ptr) {
        free(*ptr);
    }
    *ptr = new_ptr;
    return aerror(*ptr, size);
}

int fftw_malloc_s(fftw_complex **ptr, size_t size) {
    *ptr = fftw_malloc(size);
    return aerror((void *) *ptr, size);
}

int calloc_align_s(void **ptr, size_t nelem, size_t elsize) {
    size_t size = nelem*elsize;
    
    int ans = posix_memalign(ptr, alignment, size);
    if (ans != 0) {
        return aerror(NULL, size);
    }

    memset(*ptr, 0, size);

    return aerror(*ptr, size);
}

int malloc_align_s(void **ptr, size_t size) {
    int ans = posix_memalign(ptr, alignment, size);
    if (ans != 0) {
        return aerror(NULL, size);
    }
    return aerror(*ptr, size);
}
