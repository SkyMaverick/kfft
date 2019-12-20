#pragma once

#include "kfft_guts.h"

typedef kfft_rplan_t* kfft_rplan;

kfft_rplan
kfft_rconfig(int nfft, int inverse_fft, int level, void* mem, size_t* lenmem);
void
kfft_rfree(kfft_rplan* cfg);
