#pragma once

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "kfft_config.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef double kfft_scalar;

typedef struct {
    kfft_scalar r;
    kfft_scalar i;
} kfft_cpx;

enum {
    KFFT_INFO_TRACE = 1 << 0,
    KFFT_INFO_USE_SIMD = 1 << 1,
    KFFT_INFO_USE_ALLOCA = 1 << 2,
    KFFT_INFO_USE_SYSMATH = 1 << 3,
    KFFT_INFO_RADER_ALGO = 1 << 4,
    KFFT_INFO_MEMLESS_MODE = 1 << 5,
};

enum {
    KFFT_FLAG_NORMAL = 0,
    KFFT_FLAG_INVERSE = 1 << 0,
    KFFT_FLAG_RENEW = 1 << 1,
    KFFT_FLAG_GENERIC = 1 << 2,
    KFFT_FLAG_GENERIC_ONLY = 1 << 3
};

typedef struct {
    uint16_t vmajor;
    uint16_t vminor;
    uint16_t vpatch;

    uint16_t flags;
} kfft_info_t;

enum { KFFT_ERROR_ALLOCATE_FAIL, KFFT_ERROR_BUFFER_OVVERFLOW };

#include "incs/kfft_macro.h"
#include "incs/kfft_system.h"
#include "incs/kfft_alloc.h"

#include "incs/kfft_cpx.h"
#include "incs/kfft_real.h"

#define KFFT_PLAN_ALLOCATOR(X) (*((kfft_pool_t**)(X)))

KFFT_API void
kfft_info(kfft_info_t* info);
KFFT_API void
kfft_cleanup(uintptr_t mem);

#define kfft_free(X) kfft_cleanup((uintptr_t)(X))

#ifdef __cplusplus
}
#endif
