#pragma once

#include <stdbool.h>

#define HW_ALIGN_X8664_SSE 16
#define HW_ALIGN_X8664_AVX 32
#define HW_ALIGN_ARM_NEON 16

enum {
    HW_ARCH_UKNW = 0,
    HW_ARCH_X86 = 1u << 0,
    HW_ARCH_X64 = 1u << 1,
    HW_ARCH_ARM = 1u << 2,
};

/* x86 architecture extensions */
#define HW_MMX 1u << 0
#define HW_x64 1u << 1
#define HW_ABM 1u << 2
#define HW_RDRAND 1u << 3
#define HW_BMI1 1u << 4
#define HW_BMI2 1u << 5
#define HW_ADX 1u << 6
#define HW_PREFETCHWT1 1u << 7
#define HW_MPX 1u << 8
#define HW_SSE 1u << 9
#define HW_SSE2 1u << 10
#define HW_SSE3 1u << 11
#define HW_SSSE3 1u << 12
#define HW_SSE41 1u << 13
#define HW_SSE42 1u << 14
#define HW_SSE4a 1u << 15
#define HW_AES 1u << 16
#define HW_SHA 1u << 17
#define HW_AVX 1u << 18
#define HW_XOP 1u << 19
#define HW_FMA3 1u << 20
#define HW_FMA4 1u << 21
#define HW_AVX2 1u << 22
#define HW_AVX512_F 1u << 23
#define HW_AVX512_PF 1u << 24
#define HW_AVX512_ER 1u << 25
#define HW_AVX512_CD 1u << 26
#define HW_AVX512_VL 1u << 27
#define HW_AVX512_BW 1u << 28
#define HW_AVX512_DQ 1u << 29
#define HW_AVX512_IFMA 1u << 30
#define HW_AVX512_VBMI 1u << 31

KFFT_API kfft_simd_t
kfft_simd_analize(void);

static inline bool
kfft_simd_check(const kfft_simd_t S, uint32_t flag) {
    return (S.ext & flag);
}

static inline uint8_t
kfft_simd_align(const kfft_simd_t S) {
    switch (S.arch) {
    case HW_ARCH_X86:
    case HW_ARCH_X64: {
        return
#if defined(KFFT_SIMD_AVX2_SUPPORT)
            (kfft_simd_check(S, HW_AVX2))
                ? HW_ALIGN_X8664_AVX
                :
#endif
#if defined(KFFT_SIMD_AVX_SUPPORT)
                (kfft_simd_check(S, HW_AVX))
                    ? HW_ALIGN_X8664_AVX
                    :
#endif
#if defined(KFFT_SIMD_SSE_SUPPORT)
                    (kfft_simd_check(S, HW_SSE | HW_SSE2)) ? HW_ALIGN_X8664_SSE :
#endif
                                                           0;
    }
    case HW_ARCH_ARM: {
#if defined(KFFT_SIMD_NEON_SUPPORT)
        // TODO detect NEON
        return HW_ALIGN_ARM_NEON;
#endif
    }
    case HW_ARCH_UKNW:
        return 0;
    }
    return 0;
}

KFFT_API kfft_simd_t
kfft_simd_info(void);

#if defined(KFFT_DYNAPI_ENABLE)
typedef kfft_simd_t (*kfft_callback_simd_analize)(void);
typedef kfft_simd_t (*kfft_callback_simd_info)(void);
#endif /* KFFT_DYNAPI_ENABLE */
