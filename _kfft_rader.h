#pragma once

#include <stdio.h>

#include "_kfft_guts.h"
#include "kfft_core.h"
#include "kfft.h"

static inline __fft_cfg
__fft_recursive_alloc (const __fft_cfg st,
                       int p)
{
    return __kiss_fft_config ( p, st->inverse, st->delta, st->step, st->level + 1, NULL, NULL);
}

static void kf_rader (
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const __fft_cfg st,
        int m,
        int p
        )
{
    kfft_trace ("%s:\t m: %d\tp: %d\tlvl: %d\n", "Generic FFT use Rader", m, p, st->level);
    // TODO Add Rader algo and strart recursive FFT with new buffer
    __fft_cfg st_req = __fft_recursive_alloc (st, p);
    
    // TODO Traditional FFT with new buffer
    
    __kiss_fft_free(st_req);
}


