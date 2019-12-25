#pragma once

#include "kfft_guts.h"

void
kf_rader(kfft_cpx* Fout, const size_t fstride, const kfft_kplan_t* st, int m, int p);
