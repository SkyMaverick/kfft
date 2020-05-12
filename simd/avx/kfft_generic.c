/* perform the butterfly for one stage of a mixed radix FFT */
#include "kfft_generic.h"
#include "kfft_math_intern.h"

kfft_return_t
FUNC_AVX(std_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                          const kfft_comp_t* st, uint32_t u, uint32_t m, uint32_t p) {
    kfft_trace_core(st->level, "[GenStd (AVX)] fstride - %zu | m - %u | p - %u\n", fstride, m, p);
    kfft_return_t ret = KFFT_RET_SUCCESS;

    // TODO

    return ret;
}
