#pragma once 

#include <stdio.h>
#include "config.h"

#if (defined KFFT_TRACE)
    #define kfft_trace(fmt, ...) fprintf(stdout, fmt, __VA_ARGS__)
    #define kfft_sztrace(msg, X)                                                                   \
        do {                                                                                       \
            register double nmem = (double)(X);                                                    \
            if (nmem < 0x0400) {                                                                   \
                kfft_trace(msg "%4.2f byte\n", nmem);                                              \
            } else if ((nmem /= 0x0400) < 0x0400) {                                                \
                kfft_trace(msg "%4.2f Kbyte\n", nmem);                                             \
            } else if ((nmem /= 0x0400) < 0x0400) {                                                \
                kfft_trace(msg "%4.2f Mbyte\n", nmem);                                             \
            } else if ((nmem /= 0x0400) < 0x0400) {                                                \
                kfft_trace(msg "%4.2f Gbyte\n", nmem);                                             \
            }                                                                                      \
        } while (0)

#else
    #define kfft_trace(fmt, ...) // noop
    #define kfft_sztrace(msg, X) // noop
#endif
