#include "kfa_defs.h"

static void
display_flags(kfft_info_t* info) {
    fprintf(stdout, "%s - %s\n", "Enable trace messages",
            (info->flags & KFFT_INFO_TRACE) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use SIMD instructions",
            (info->flags & KFFT_INFO_USE_SIMD) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use alloca() function",
            (info->flags & KFFT_INFO_USE_ALLOCA) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use system math functions",
            (info->flags & KFFT_INFO_USE_SYSMATH) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use OpenMP functions",
            (info->flags & KFFT_INFO_USE_OPENMP) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use Rader algoritm",
            (info->flags & KFFT_INFO_RADER_ALGO) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable lesser memory mode",
            (info->flags & KFFT_INFO_MEMLESS_MODE) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable half scalar mode",
            (info->flags & KFFT_INFO_HALF_SCALAR) ? "YES" : "NO");
}

#if defined(KFFT_USE_SIMD)
static void
display_simd(state_t* st) {

    KCB_TYPE(simd_analize) cb_func = KFFT_CALLBACK(st, simd_analize);
    if (cb_func) {
        fprintf(stdout, "\n%s:\n", "Found CPU extensions");

        kfft_simd_t s_info = kfft_simd_analize();
    #if defined(KFFT_ARCH_X86)
        //  Misc.
        // clang-format off
        if (kfft_simd_check(s_info, HW_MMX        )) fprintf (stdout, "%s ", "MMX");
        if (kfft_simd_check(s_info, HW_x64        )) fprintf (stdout, "%s ", "x64");
        if (kfft_simd_check(s_info, HW_ABM        )) fprintf (stdout, "%s ", "ABM");
        if (kfft_simd_check(s_info, HW_RDRAND     )) fprintf (stdout, "%s ", "RDRAND");
        if (kfft_simd_check(s_info, HW_BMI1       )) fprintf (stdout, "%s ", "BMI1");
        if (kfft_simd_check(s_info, HW_BMI2       )) fprintf (stdout, "%s ", "BMI2");
        if (kfft_simd_check(s_info, HW_ADX        )) fprintf (stdout, "%s ", "ADX");
        if (kfft_simd_check(s_info, HW_PREFETCHWT1)) fprintf (stdout, "%s ", "PREFETCHWT1");
        if (kfft_simd_check(s_info, HW_MPX        )) fprintf (stdout, "%s ", "MPX");
        if (kfft_simd_check(s_info, HW_SSE        )) fprintf (stdout, "%s ", "SSE");
        if (kfft_simd_check(s_info, HW_SSE2       )) fprintf (stdout, "%s ", "SSE2");
        if (kfft_simd_check(s_info, HW_SSE3       )) fprintf (stdout, "%s ", "SSE3");
        if (kfft_simd_check(s_info, HW_SSSE3      )) fprintf (stdout, "%s ", "SSSE3");
        if (kfft_simd_check(s_info, HW_SSE41      )) fprintf (stdout, "%s ", "SSE41");
        if (kfft_simd_check(s_info, HW_SSE42      )) fprintf (stdout, "%s ", "SSE42");
        if (kfft_simd_check(s_info, HW_SSE4a      )) fprintf (stdout, "%s ", "SSE4a");
        if (kfft_simd_check(s_info, HW_AES        )) fprintf (stdout, "%s ", "AES");
        if (kfft_simd_check(s_info, HW_SHA        )) fprintf (stdout, "%s ", "SHA");
        if (kfft_simd_check(s_info, HW_AVX        )) fprintf (stdout, "%s ", "AVX");
        if (kfft_simd_check(s_info, HW_XOP        )) fprintf (stdout, "%s ", "XOP");
        if (kfft_simd_check(s_info, HW_FMA3       )) fprintf (stdout, "%s ", "FMA3");
        if (kfft_simd_check(s_info, HW_FMA4       )) fprintf (stdout, "%s ", "FMA4");
        if (kfft_simd_check(s_info, HW_AVX2       )) fprintf (stdout, "%s ", "AVX2");
        if (kfft_simd_check(s_info, HW_AVX512_F   )) fprintf (stdout, "%s ", "AVX512_F");
        if (kfft_simd_check(s_info, HW_AVX512_PF  )) fprintf (stdout, "%s ", "AVX512_PF");
        if (kfft_simd_check(s_info, HW_AVX512_ER  )) fprintf (stdout, "%s ", "AVX512_ER");
        if (kfft_simd_check(s_info, HW_AVX512_CD  )) fprintf (stdout, "%s ", "AVX512_CD");
        if (kfft_simd_check(s_info, HW_AVX512_VL  )) fprintf (stdout, "%s ", "AVX512_VL");
        if (kfft_simd_check(s_info, HW_AVX512_BW  )) fprintf (stdout, "%s ", "AVX512_BW");
        if (kfft_simd_check(s_info, HW_AVX512_DQ  )) fprintf (stdout, "%s ", "AVX512_DQ");
        if (kfft_simd_check(s_info, HW_AVX512_IFMA)) fprintf (stdout, "%s ", "AVX512_IFMA");
        if (kfft_simd_check(s_info, HW_AVX512_VBMI)) fprintf (stdout, "%s ", "AVX512_VBMI");
    // clang-format on
    #elif defined(KFFT_ARCH_ARM)
            // TODO
    #endif /* ARCH_X86 */

        fprintf(stdout, "\n%s:\n", "Build with CPU extensions");
        fprintf(stdout, "%s%s - %s\n", "Enable SSE support",
    #if defined(KFFT_HAVE_SSE3)
                " (with SSE3)"
    #else
                ""
    #endif
                ,
    #if defined(KFFT_SIMD_SSE_SUPPORT)
                "YES"
    #else
                "NO"
    #endif
        );
        fprintf(stdout, "%s - %s\n", "Enable AVX support",
    #if defined(KFFT_SIMD_AVX_SUPPORT)
                "YES"
    #else
                "NO"
    #endif
        );
    } else
        fprintf(stdout, "%s\n", "Don't found kfft_simd_analize() callback");
}
#endif /* KFFT_USE_SIMD */

static void
display_info(state_t* st) {
    kfft_info_t info;
    KRNL_FUNCS(st).cb_info(&info);

    fprintf(stdout, "%s version: %d.%d.%d\n", APP_NAME, VER_MAJOR, VER_MINOR, VER_PATCH);
    fprintf(stdout, "Uses libkfft version : %d.%d.%d\n\n", info.vmajor, info.vminor, info.vpatch);

    display_flags(&info);
#if defined(KFFT_USE_SIMD)
    display_simd(st);
#endif /* KFFT_USE_SIMD */
}

static void
display_help(void) {
    fprintf(stdout, "%s", help_msg);
}
