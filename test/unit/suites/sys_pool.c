#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

#include "util.h"
#include "CUnit/Basic.h"

#include "kfft.h"

#define TEST_MEMORY_SIZE 100
#define TEST_CYCLE_COUNT 10

_TEST(create_free) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    CU_ASSERT_PTR_NOT_NULL(P);
    if (P) {
        CU_ASSERT_EQUAL((P->allocated), TEST_MEMORY_SIZE + sizeof(kfft_pool_t));
        CU_ASSERT_EQUAL((P->tail - P->head), TEST_MEMORY_SIZE);
        CU_ASSERT_EQUAL(P->head, P->cur);
        kfft_pool_free(P);
    }

    STDIO_ON(A);
}

_TEST(zero) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(0);
    if (P) {
        CU_ASSERT_EQUAL(P->allocated, sizeof(kfft_pool_t));
        CU_ASSERT_PTR_EQUAL(P->tail, P->head);
        CU_ASSERT_PTR_EQUAL(P->head, P->cur);
        kfft_pool_free(P);
    }

    STDIO_ON(A);
}

_TEST(zero_overflow) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(0);
    if (P) {
        void* ret = kfft_pool_alloc(P, 1);
        CU_ASSERT_PTR_NULL(ret);
        kfft_pool_free(P);
    }

    STDIO_ON(A);
}

_TEST(alloc) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    if (P) {
        void* ret = kfft_pool_alloc(P, TEST_MEMORY_SIZE / 2);
        CU_ASSERT_PTR_NOT_NULL(ret);
        kfft_pool_free(P);
    }

    STDIO_ON(A);
}
_TEST(alloc_zero) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    if (P) {
        void* ret = kfft_pool_alloc(P, 0);
        CU_ASSERT_PTR_NULL(ret);
        kfft_pool_free(P);
    }

    STDIO_ON(A);
}
_TEST(alloc_hight) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    if (P) {
        void* ret = kfft_pool_alloc(P, TEST_MEMORY_SIZE);
        CU_ASSERT_PTR_NOT_NULL(ret);
        CU_ASSERT_PTR_EQUAL(P->cur, P->tail);

        ret = kfft_pool_alloc(P, 1);
        CU_ASSERT_PTR_NULL(ret);

        kfft_pool_free(P);
    }

    STDIO_ON(A);
}
_TEST(alloc_chain) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    if (P) {
        void* ret = kfft_pool_alloc(P, TEST_MEMORY_SIZE / 5);
        CU_ASSERT_PTR_NOT_NULL(ret);
        ret = kfft_pool_alloc(P, TEST_MEMORY_SIZE / 5);
        CU_ASSERT_PTR_NOT_NULL(ret);
        ret = kfft_pool_alloc(P, TEST_MEMORY_SIZE / 5);
        CU_ASSERT_PTR_NOT_NULL(ret);
        ret = kfft_pool_alloc(P, TEST_MEMORY_SIZE / 5);
        CU_ASSERT_PTR_NOT_NULL(ret);
        kfft_pool_free(P);
    }

    STDIO_ON(A);
}
_TEST(alloc_overflow) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    if (P) {
        void* ret = kfft_pool_alloc(P, TEST_MEMORY_SIZE + 1);
        CU_ASSERT_PTR_NULL(ret);
        CU_ASSERT_EQUAL(P->head, P->cur);
        kfft_pool_free(P);
    }

    STDIO_ON(A);
}

_TEST(zmem_size) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    size_t block_size = TEST_MEMORY_SIZE / 4;
    if (P) {
        void* A = kfft_pool_alloc(P, block_size);
        memset(A, 0xff, block_size);
        void* B = kfft_pool_alloc(P, block_size);
        memset(B, 0xff, block_size);
        void* C = kfft_pool_alloc(P, block_size);
        memset(C, 0xff, block_size);

        kfft_pool_zmem(P, B, block_size);
        CU_ASSERT_EQUAL(((uint8_t*)B)[0], 0);
        CU_ASSERT_EQUAL(((uint8_t*)B)[-1], 0xff);
        CU_ASSERT_EQUAL(((uint8_t*)C)[-1], 0);
        CU_ASSERT_EQUAL(((uint8_t*)C)[0], 0xff);

        kfft_pool_free(P);
    }

    STDIO_ON(A);
}
_TEST(zmem_overflow) {
    STDIO_OFF(A);

    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    size_t block_size = TEST_MEMORY_SIZE / 4;
    if (P) {
        void* A = kfft_pool_alloc(P, block_size);
        memset(A, 0xff, block_size);
        void* B = kfft_pool_alloc(P, block_size);
        memset(B, 0xff, block_size);
        void* C = kfft_pool_alloc(P, block_size);
        memset(C, 0xff, block_size);

        kfft_pool_zmem(P, B, TEST_MEMORY_SIZE);
        CU_ASSERT_EQUAL(((uint8_t*)B)[0], 0xff);

        kfft_pool_free(P);
    }

    STDIO_ON(A);
}

_TEST(clear) {
    STDIO_OFF(A);

    srand(time(NULL));
    kfft_pool_t* P = kfft_pool_create(TEST_MEMORY_SIZE);
    size_t block_size = TEST_MEMORY_SIZE / 4;
    if (P) {
        void* A = kfft_pool_alloc(P, block_size);
        memset(A, 0xaa, block_size);
        void* B = kfft_pool_alloc(P, block_size);
        memset(B, 0xcc, block_size);
        void* C = kfft_pool_alloc(P, block_size);
        memset(C, 0xff, block_size);

        kfft_pool_clear(P);
        size_t sum = 0;
        for (unsigned i = 0; i < TEST_CYCLE_COUNT; i++)
            sum += *((uint8_t*)P->head + (rand() % (P->head - P->tail)));
        CU_ASSERT_EQUAL(sum, 0);
        CU_ASSERT_PTR_EQUAL(P->head, P->cur);

        kfft_pool_free(P);
    }

    STDIO_ON(A);
}

void
_run_suite(void) {
    // ... all startup functionality
    CU_pSuite suite = CUnitCreateSuite("plan internal allocator");
    if (suite) {
        _ADD_TEST(suite, create_free);
        _ADD_TEST(suite, zero);
        _ADD_TEST(suite, zero_overflow);

        _ADD_TEST(suite, alloc);
        _ADD_TEST(suite, alloc_zero);
        _ADD_TEST(suite, alloc_hight);
        _ADD_TEST(suite, alloc_chain);
        _ADD_TEST(suite, alloc_overflow);

        _ADD_TEST(suite, zmem_size);
        _ADD_TEST(suite, zmem_overflow);

        _ADD_TEST(suite, clear);
    }
    return;
}
