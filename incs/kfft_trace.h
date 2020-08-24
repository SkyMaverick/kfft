#pragma once

#if (defined KFFT_TRACE)
    #include <stdio.h>
    #if (defined(KFFT_USE_OPENMP) && !defined(KFFT_OMP_ISBLOCKED))
        #include <omp.h>
        #define kfft_trace(fmt, ...)                                                               \
            fprintf(stdout, "<Thr.%d> " fmt, omp_get_thread_num(), __VA_ARGS__)
        #define kfft_trace_raw(fmt, ...) fprintf(stdout, fmt, __VA_ARGS__)
    #else
        #define kfft_trace(fmt, ...) fprintf(stdout, fmt, __VA_ARGS__)
        #define kfft_trace_raw(fmt, ...) kfft_trace(fmt, __VA_ARGS__)
    #endif
#else
    #define kfft_trace(fmt, ...)     // noop
    #define kfft_trace_raw(fmt, ...) // noop
#endif
