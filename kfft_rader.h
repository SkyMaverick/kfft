#pragma once

#include "kfft_guts.h"

#if defined(KFFT_RADER_ALGO)

void
kf_rader(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, uint32_t p);

#endif /* KFFT_READER_ALGO */
