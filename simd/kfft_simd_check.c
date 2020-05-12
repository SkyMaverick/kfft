#include "kfft_macro.h"
#include "kfft_config.h"

#include "kfft_simd.h"

#include <inttypes.h>
#include <stdbool.h>

#if defined(KFFT_ARCH_X86)
    #ifdef KFFT_OS_WINDOWS
        #define cpuid(info, x) __cpuidex(info, x, 0)
static inline __int64
xgetbv(unsigned int x) {
    return _xgetbv(x);
}
    #else /* WINDOWS */
        #include <cpuid.h>
        #if defined(__pic__) && defined(__i386__)
void
cpuid(int cpu_info[4], int info_type) {
    __asm__ volatile("mov %%ebx, %%edi\n"
                     "cpuid\n"
                     "xchg %%edi, %%ebx\n"
                     : "=a"(cpu_info[0]), "=D"(cpu_info[1]), "=c"(cpu_info[2]), "=d"(cpu_info[3])
                     : "a"(info_type), "c"(0));
}
        #else  /* PIC and i386 */
void
cpuid(int cpu_info[4], int info_type) {
    __asm__ volatile("cpuid\n"
                     : "=a"(cpu_info[0]), "=b"(cpu_info[1]), "=c"(cpu_info[2]), "=d"(cpu_info[3])
                     : "a"(info_type), "c"(0));
}
        #endif /* PIC and i386 */
static inline uint64_t
xgetbv(unsigned int index) {
    uint32_t eax, edx;
    __asm__ __volatile__("xgetbv" : "=a"(eax), "=d"(edx) : "c"(index));
    return ((uint64_t)edx << 32) | eax;
}
        #define _XCR_XFEATURE_ENABLED_MASK 0
    #endif /* Win32 */

static inline bool
protect_avx() {
    int info[4];
    cpuid(info, 1);

    return (((info[2] & (1 << 27)) != 0) &&                        /* XSAVE && XSTORE */
            ((info[2] & (1 << 28)) != 0))                          /* AVX support */
               ? (xgetbv(_XCR_XFEATURE_ENABLED_MASK) & 0x6) == 0x6 /* OS support check */
               : false;
}

static inline bool
protect_avx512(void) {
    return protect_avx() ? (xgetbv(_XCR_XFEATURE_ENABLED_MASK) & 0xe6) == 0xe6 : false;
}

static uint32_t
x86_exts(void) {
    uint32_t ret = 0;

    int info[4];
    cpuid(info, 0);

    int nIds = info[0];
    cpuid(info, 0x80000000);

    unsigned nExIds = info[0];

    //  Detect Features
    if (nIds >= 0x00000001) {
        cpuid(info, 0x00000001);
        ret |= (info[3] & ((int)1 << 23)) ? HW_MMX : 0;
        ret |= (info[3] & ((int)1 << 25)) ? HW_SSE : 0;
        ret |= (info[3] & ((int)1 << 26)) ? HW_SSE2 : 0;
        ret |= (info[2] & ((int)1 << 0)) ? HW_SSE3 : 0;

        ret |= (info[2] & ((int)1 << 9)) ? HW_SSSE3 : 0;
        ret |= (info[2] & ((int)1 << 19)) ? HW_SSE41 : 0;
        ret |= (info[2] & ((int)1 << 20)) ? HW_SSE42 : 0;
        ret |= (info[2] & ((int)1 << 25)) ? HW_AES : 0;

        ret |= ((info[2] & ((int)1 << 28)) && (protect_avx())) ? HW_AVX : 0;
        ret |= (info[2] & ((int)1 << 12)) ? HW_FMA3 : 0;

        ret |= (info[2] & ((int)1 << 30)) ? HW_RDRAND : 0;
    }
    if (nIds >= 0x00000007) {
        cpuid(info, 0x00000007);
        ret |= ((info[1] & ((int)1 << 5)) && (protect_avx())) ? HW_AVX2 : 0;

        ret |= (info[1] & ((int)1 << 3)) ? HW_BMI1 : 0;
        ret |= (info[1] & ((int)1 << 8)) ? HW_BMI2 : 0;
        ret |= (info[1] & ((int)1 << 19)) ? HW_ADX : 0;
        ret |= (info[1] & ((int)1 << 14)) ? HW_MPX : 0;
        ret |= (info[1] & ((int)1 << 29)) ? HW_SHA : 0;
        ret |= (info[2] & ((int)1 << 0)) ? HW_PREFETCHWT1 : 0;

        ret |= (info[1] & ((int)1 << 16)) && (protect_avx512()) ? HW_AVX512_F : 0;
        ret |= (info[1] & ((int)1 << 28)) && (protect_avx512()) ? HW_AVX512_CD : 0;
        ret |= (info[1] & ((int)1 << 26)) && (protect_avx512()) ? HW_AVX512_PF : 0;
        ret |= (info[1] & ((int)1 << 27)) && (protect_avx512()) ? HW_AVX512_ER : 0;
        ret |= (info[1] & ((int)1 << 31)) && (protect_avx512()) ? HW_AVX512_VL : 0;
        ret |= (info[1] & ((int)1 << 30)) && (protect_avx512()) ? HW_AVX512_BW : 0;
        ret |= (info[1] & ((int)1 << 17)) && (protect_avx512()) ? HW_AVX512_DQ : 0;
        ret |= (info[1] & ((int)1 << 21)) && (protect_avx512()) ? HW_AVX512_IFMA : 0;
        ret |= (info[2] & ((int)1 << 1)) && (protect_avx512()) ? HW_AVX512_VBMI : 0;
    }
    if (nExIds >= 0x80000001) {
        cpuid(info, 0x80000001);

        ret |= (info[3] & ((int)1 << 29)) ? HW_x64 : 0;
        ret |= (info[2] & ((int)1 << 5)) ? HW_ABM : 0;
        ret |= (info[2] & ((int)1 << 6)) ? HW_SSE4a : 0;
        ret |= (info[2] & ((int)1 << 16)) ? HW_FMA4 : 0;
        ret |= (info[2] & ((int)1 << 11)) ? HW_XOP : 0;
    }
    return ret;
}

static uint8_t
x86_arch(void) {
    uint8_t ret = 0;
    return ret;
}
#endif /*KFFT_ARCH_X86*/

/******************************************************************************/

#if defined(KFFT_ARCH_ARM)
// TODO
#endif /* KFFT_ARCH_ARM */

KFFT_API kfft_simd_t
kfft_simd_analize(void) {
    kfft_simd_t S = {0, 0, 0};
#if defined(KFFT_ARCH_X86)
    S.arch = x86_arch();
    S.ext = x86_exts();
#elif defined(KFFT_ARCH_ARM)
    S.arch = HW_ARCH_ARM;
#else
    S.arch = HW_ARCH_UKNW;
#endif
    return S;
}
