#pragma once

#include <stdlib.h>
#include "kfft_guts.h"

kfft_pool_t*
kfft_allocator_create(size_t size);

kfft_pool_t*
kfft_allocator_init(void* mem, size_t nmem);

void*
kfft_internal_alloc(kfft_pool_t* A, size_t nmem);

size_t
kfft_allocator_empty(kfft_pool_t* A);

void
kfft_internal_zmem(kfft_pool_t* A, void* ptr, size_t size);

void
kfft_allocator_clear(kfft_pool_t* A);

void
kfft_allocator_free(kfft_pool_t** A);
