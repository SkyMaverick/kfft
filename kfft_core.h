#pragma once

#include "_kfft_guts.h"

__fft_cfg
__kfft_config (int nfft,
                   int inverse_fft,
                   int level,
                   void * mem,
                   size_t * lenmem );

void
__kfft(__fft_cfg cfg,const kfft_cpx *fin,kfft_cpx *fout);

#define __kfft_free(X)  \
    do {                    \
        free(X);            \
        X=NULL;             \
    }while(0)
