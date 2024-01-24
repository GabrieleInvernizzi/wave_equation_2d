#pragma once

#include <stddef.h>

typedef struct {
    void* arena;
    void* idx;
    void* max_idx;
} Arena;


int arena_init(Arena* a, size_t cap);
void* arena_alloc(Arena* a, size_t n_bytes);
int area_deinit(Arena* a);


