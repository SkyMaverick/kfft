#pragma once

#include "config.h"

#if (defined KFFT_TRACE)
    #include <stdio.h>
    #define kfft_trace(fmt, ...) fprintf(stdout, fmt, __VA_ARGS__)
#else
    #define kfft_trace(fmt, ...) // noop
#endif
