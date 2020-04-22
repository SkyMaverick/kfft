// clang-format off
static void
display_simd (void) {
    kfft_simd_t s_info = kfft_simd_analize();
    //  Misc.
    fprintf(stdout, "%s - %s\n", "HW_MMX         ", kfft_simd_check(s_info, HW_MMX        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_x64         ", kfft_simd_check(s_info, HW_x64        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_ABM         ", kfft_simd_check(s_info, HW_ABM        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_RDRAND      ", kfft_simd_check(s_info, HW_RDRAND     ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_BMI1        ", kfft_simd_check(s_info, HW_BMI1       ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_BMI2        ", kfft_simd_check(s_info, HW_BMI2       ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_ADX         ", kfft_simd_check(s_info, HW_ADX        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_PREFETCHWT1 ", kfft_simd_check(s_info, HW_PREFETCHWT1) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_MPX         ", kfft_simd_check(s_info, HW_MPX        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_SSE         ", kfft_simd_check(s_info, HW_SSE        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_SSE2        ", kfft_simd_check(s_info, HW_SSE2       ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_SSE3        ", kfft_simd_check(s_info, HW_SSE3       ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_SSSE3       ", kfft_simd_check(s_info, HW_SSSE3      ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_SSE41       ", kfft_simd_check(s_info, HW_SSE41      ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_SSE42       ", kfft_simd_check(s_info, HW_SSE42      ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_SSE4a       ", kfft_simd_check(s_info, HW_SSE4a      ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AES         ", kfft_simd_check(s_info, HW_AES        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_SHA         ", kfft_simd_check(s_info, HW_SHA        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX         ", kfft_simd_check(s_info, HW_AVX        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_XOP         ", kfft_simd_check(s_info, HW_XOP        ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_FMA3        ", kfft_simd_check(s_info, HW_FMA3       ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_FMA4        ", kfft_simd_check(s_info, HW_FMA4       ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX2        ", kfft_simd_check(s_info, HW_AVX2       ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_F    ", kfft_simd_check(s_info, HW_AVX512_F   ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_PF   ", kfft_simd_check(s_info, HW_AVX512_PF  ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_ER   ", kfft_simd_check(s_info, HW_AVX512_ER  ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_CD   ", kfft_simd_check(s_info, HW_AVX512_CD  ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_VL   ", kfft_simd_check(s_info, HW_AVX512_VL  ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_BW   ", kfft_simd_check(s_info, HW_AVX512_BW  ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_DQ   ", kfft_simd_check(s_info, HW_AVX512_DQ  ) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_IFMA ", kfft_simd_check(s_info, HW_AVX512_IFMA) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "HW_AVX512_VBMI ", kfft_simd_check(s_info, HW_AVX512_VBMI) ? "YES" : "NO");
}                                                         
// clang-format on

static void
display_info(void) {
    kfft_info_t info;
    kfft_info(&info);

    fprintf(stdout, "%s version: %d.%d.%d\n", APP_NAME, VER_MAJOR, VER_MINOR, VER_PATCH);
    fprintf(stdout, "Uses libkfft version : %d.%d.%d\n\n", info.vmajor, info.vminor, info.vpatch);

    fprintf(stdout, "%s - %s\n", "Enable trace messages",
            (info.flags & KFFT_INFO_TRACE) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use SIMD instructions",
            (info.flags & KFFT_INFO_USE_SIMD) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use alloca() function",
            (info.flags & KFFT_INFO_USE_ALLOCA) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use system math functions",
            (info.flags & KFFT_INFO_USE_SYSMATH) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use OpenMP functions",
            (info.flags & KFFT_INFO_USE_OPENMP) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use Rader algoritm",
            (info.flags & KFFT_INFO_RADER_ALGO) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable lesser memory mode",
            (info.flags & KFFT_INFO_MEMLESS_MODE) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable half scalar mode",
            (info.flags & KFFT_INFO_HALF_SCALAR) ? "YES" : "NO");

    fprintf(stdout, "%s:\n", "Check CPU instructions");
    display_simd();

}

static void
display_help(void) {
    fprintf(stdout, "%s", help_msg);
}
