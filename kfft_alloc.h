#pragma once

#include <stdlib.h>

typedef struct {
    size_t allocated;   //area size
    void* head;         //current head address
    void* tail;         //current tail address

    void* cur;          // cursor pointer

    char area[1];
} kfft_allocator_t;

kfft_allocator_t*
kfft_allocator_create (size_t size);

void*
kfft_internal_alloc (size_t nmem);

void
kfft_internal_zmem (void* ptr, size_t size);

void
kfft_allocator_free (kfft_allocator_t* A);
