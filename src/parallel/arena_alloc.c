#include "arena_alloc.h"

#include <stdio.h>
#include <stdlib.h>

int arena_init(Arena *a, size_t cap) {

    a->arena = malloc(cap);
    a->idx = a->arena;
    a->max_idx = a->arena;
    if (!a->arena) return 1;
   
    a->max_idx = a->idx + cap;

    return 0;
}

int area_deinit(Arena *a) {
    if (!a->arena)
        return 1;

    free(a->arena);
    a->idx = NULL;
    a->max_idx = NULL;

    return 0;
}

void *arena_alloc(Arena *a, size_t n_bytes) {

    if (a->idx + n_bytes > a->max_idx)
        return NULL;
    
    void *ret = a->idx;
    a->idx += n_bytes;

    return ret;
}
