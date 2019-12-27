#pragma once

#if defined(KFFT_RADER_ALGO)

    #include "kfft_guts.h"

void
kf_rader(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, uint32_t m, uint32_t p);

#endif /* KFFT_READER_ALGO */
