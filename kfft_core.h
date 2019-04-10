#pragma once

#include "_kfft_guts.h"

__fft_cfg
__kiss_fft_config (int nfft,
                   int inverse_fft,
                   int delta,
                   int step,
                   int level,
                   void * mem,
                   size_t * lenmem );

void
__kiss_fft(__fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout);
