/* perform the butterfly for one stage of a mixed radix FFT */
#include "kfft_generic.h"
#include "kfft_math_intern.h"

kfft_return_t
FUNC_SSE(std_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                          const kfft_comp_t* st, uint32_t u, uint32_t m, uint32_t p) {
    kfft_trace_core(st->level, "[GenStd (SSE)] fstride - %zu | m - %u | p - %u\n", fstride, m, p);
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t k = u, q1, q;
    __m128 T;

    for (q1 = 0; q1 < p; ++q1, k += m)
        _mm_store_ps((float*)&Ftmp[q1], _mm_load_ps((float*)&Fout[k]));

    k = u;

    for (q1 = 0; q1 < p; ++q1, k += m) {
        uint32_t twidx = 0;
        __m128 FOK = _mm_load_ps((float*)&Ftmp[0]);
        for (q = 1; q < p; ++q) {
            twidx += fstride * k;
            if (twidx >= st->nfft)
                twidx -= st->nfft;
            kfft_cpx ctw = TWIDDLE(twidx, st);
            C_MUL_SSE(T, _mm_load_ps((float*)&Ftmp[q]), _mm_load_ps((float*)&ctw));
            C_ADD_SSE(FOK, FOK, T);
            _mm_store_ps((float*)&Fout[k], FOK);
        }
    }
    return ret;
}
