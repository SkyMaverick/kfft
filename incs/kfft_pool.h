#pragma once

#include <stdint.h>
#include <stdlib.h>

/*!
  \file
  \brief Small internal plan allocator

  Use this fast allocator for internal plan configuration structures.
  This MUST have each plan object (or group).
  \hint Estimate the amount of space needed before creating the pool.
  */
typedef struct {
    size_t allocated; ///< pool area size
    uint8_t align;    ///< memory align in pool
    kfft_simd_t vex;  ///< system VEX info

    uint8_t* head; ///< current head address
    uint8_t* tail; ///< current tail address
    uint8_t* cur;  ///< cursor current pointer

    uint8_t area[1]; ///< pool data
} kfft_pool_t;

/*!
    Create new allocator object
    \param[in] size - pool size in byte
    \result pool object pointer
  */
kfft_pool_t*
kfft_pool_create(const size_t size);

/*!
    Allocate memory area with pool allocator
    \param[in] A - pool object
    \param[in] nmem - allocation memory size (byte)
    \result memory area pointer
 */
void*
kfft_pool_alloc(kfft_pool_t* A, const size_t nmem);

/*!
    Estimates the amount of free space in the allocator
    \param[in] A - pool object pointer (maybe NULL)
    \result free space (bytes)
 */
size_t
kfft_pool_empty(const kfft_pool_t* A);

/*!
    Zero memory block in the allocator
    \param[in] A - pool object pointer
    \param[in] ptr - memory block pointer
    \param[in] size  - memory block size (bytes)
 */
void
kfft_pool_zmem(const kfft_pool_t* A, void* ptr, const size_t size);

/*!
    Renew allocator. Zero all memory and clear pointers
    \param[in] A - pool object pointer
 */
void
kfft_pool_clear(kfft_pool_t* A);
/*!
    Release allocator.
    \param[in] A - pool object pointer
 */
void
kfft_pool_free(kfft_pool_t* A);

/*!
    Release allocator and NULL variable.
    \param[in] A - pool object pointer
 */
#define kfft_pool_free_and_null(A)                                                                 \
    do {                                                                                           \
        kfft_pool_free(A);                                                                         \
        A = NULL;                                                                                  \
    } while (0)
