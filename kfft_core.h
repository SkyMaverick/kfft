#pragma once

#include "_kfft_guts.h"

kfft_kplan_p
kfft_kconfig(int nfft, int inverse_fft, int level, void* mem, size_t* lenmem);

void
__kfft(kfft_kplan_p cfg, const kfft_cpx* fin, kfft_cpx* fout);

#define kfft_kfree(X)                                                                              \
    do {                                                                                           \
        free(X);                                                                                   \
        X = NULL;                                                                                  \
    } while (0)
