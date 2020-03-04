#pragma once

#include <stdint.h>
#include <stdlib.h>

typedef struct {
    size_t allocated; // area size
    uint8_t* head;    // current head address
    uint8_t* tail;    // current tail address

    uint8_t* cur; // cursor pointer

    uint8_t area[1];
} kfft_pool_t;

kfft_pool_t*
kfft_allocator_create(const size_t size);

kfft_pool_t*
kfft_allocator_init(void* mem, const size_t nmem);

void*
kfft_internal_alloc(kfft_pool_t* A, const size_t nmem);

size_t
kfft_allocator_empty(const kfft_pool_t* A);

void
kfft_internal_zmem(const kfft_pool_t* A, void* ptr, const size_t size);

void
kfft_allocator_clear(kfft_pool_t* A);

void
kfft_allocator_free(kfft_pool_t* A);

#define kfft_allocator_free_and_null(A)                                                            \
    do {                                                                                           \
        kfft_allocator_free(A);                                                                    \
        A = NULL;                                                                                  \
    } while (0)
