#pragma once

#include <stdint.h>
#include <stdlib.h>

typedef struct {
    size_t allocated; // area size
    uint8_t align;    // memory align in pool
    kfft_simd_t vex;  // system extensions info

    uint8_t* head; // current head address
    uint8_t* tail; // current tail address
    uint8_t* cur;  // cursor pointer

    uint8_t area[1];
} kfft_pool_t;

kfft_pool_t*
kfft_pool_create(const size_t size);

void*
kfft_pool_alloc(kfft_pool_t* A, const size_t nmem);

size_t
kfft_pool_empty(const kfft_pool_t* A);

void
kfft_pool_zmem(const kfft_pool_t* A, void* ptr, const size_t size);

void
kfft_pool_clear(kfft_pool_t* A);

void
kfft_pool_free(kfft_pool_t* A);

#define kfft_pool_free_and_null(A)                                                                 \
    do {                                                                                           \
        kfft_pool_free(A);                                                                         \
        A = NULL;                                                                                  \
    } while (0)
