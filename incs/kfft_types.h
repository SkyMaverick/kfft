#pragma once

#include <inttypes.h>

#if defined(KFFT_OS_WINDOWS)
    #include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

#if defined(KFFT_HALF_SCALAR)
typedef float kfft_scalar;
#else
typedef double kfft_scalar;
#endif

typedef struct {
    kfft_scalar r;
    kfft_scalar i;
} kfft_cpx;

enum {
    KFFT_INFO_TRACE = 1 << 0,
    KFFT_INFO_USE_SIMD = 1 << 1,
    KFFT_INFO_USE_ALLOCA = 1 << 2,
    KFFT_INFO_USE_SYSMATH = 1 << 3,
    KFFT_INFO_USE_OPENMP = 1 << 4,
    KFFT_INFO_RADER_ALGO = 1 << 5,
    KFFT_INFO_MEMLESS_MODE = 1 << 6,
    KFFT_INFO_HALF_SCALAR = 1 << 7,
};

enum {
    KFFT_FLAG_NORMAL = 0,
    KFFT_FLAG_INVERSE = 1 << 0,
    KFFT_FLAG_RENEW = 1 << 1,
    KFFT_FLAG_GENERIC = 1 << 2,
    KFFT_FLAG_GENERIC_ONLY = 1 << 3,
};

typedef struct {
    uint16_t vmajor;
    uint16_t vminor;
    uint16_t vpatch;

    uint16_t flags;
} kfft_info_t;

enum {
    KFFT_RET_SUCCESS = 0x0000,
    KFFT_RET_ALLOC_FAIL = 0x0001,
    KFFT_RET_BUFFER_FAIL = 0x0002,
    KFFT_RET_FREE_NULL = 0x0003,
    KFFT_RET_IMPROPER_PLAN = 0x0004,
    KFFT_RET_BADARGUMENTS = 0x0005
};
typedef unsigned kfft_return_t;

typedef struct {
    uint8_t arch; // Architecture ID
    uint32_t ext; // HW extensions extensionse (with operation system correct)
    uint16_t pad; // Paddung for uint64_t type size
} kfft_simd_t;
