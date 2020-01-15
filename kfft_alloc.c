#include "kfft_alloc.h"

kfft_pool_t*
kfft_allocator_init(void* mem, size_t nmem) {
    kfft_pool_t* ret = (kfft_pool_t*)mem;
    if (mem && (nmem > sizeof(kfft_pool_t))) {
        ret->allocated = nmem;

        ret->head = ret->area;
        ret->tail = (uint8_t*)ret + nmem;

        ret->cur = ret->head;
    }
    return ret;
}

kfft_pool_t*
kfft_allocator_create(size_t size) {
    size_t memneed = sizeof(kfft_pool_t) + size;
    kfft_pool_t* ret = KFFT_MALLOC(memneed);
    return kfft_allocator_init(ret, memneed);
}

void*
kfft_internal_alloc(kfft_pool_t* A, size_t nmem) {
    uint8_t* ret = NULL;

    if (A && (nmem > 0)) {
        if (A->cur + nmem <= A->tail) {
            ret = A->cur;
            A->cur += nmem;
        }
    }
    return (void*)ret;
}

void
kfft_internal_zmem(kfft_pool_t* A, void* ptr, size_t size) {
    if (A && ptr && (size > 0))
        if ((uint8_t*)ptr + size <= A->tail)
            KFFT_ZEROMEM(ptr, size);
}

void
kfft_allocator_clear(kfft_pool_t* A) {
    kfft_internal_zmem(A, (void*)(A->area), A->tail - A->head);
    A->cur = A->head;
}

void
kfft_allocator_free(kfft_pool_t** A) {
    if (*A) {
        KFFT_FREE(*A);
    }
}

size_t
kfft_allocator_empty(kfft_pool_t* A) {
    return (A) ? A->tail - A->cur : 0;
}
