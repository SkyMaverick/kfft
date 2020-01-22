#pragma once

#include <stdlib.h>
#include "kfft_guts.h"

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
kfft_allocator_free(kfft_pool_t** A);
