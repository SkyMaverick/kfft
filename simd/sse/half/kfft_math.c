#include "kfft_simd_compat.h"
#include "kfft_simd.h"

#include "kfft_math_intern.h"

void
FUNC_SSE(kfft_math_hadamard_cpx)(kfft_cpx* Fout, const kfft_cpx* Fin, uint32_t size) {
    // FIXME Now primitive algorithm
    __m128 Fo, Fi, T;
    for (uint32_t i = 0; i < size; i++) {
        Fo = CLOAD1212(&Fout[i]);
        /* FIXME if use _mm_load_ps debug build crash with SIGSEGV */
        Fi = _mm_set_ps((float)Fin[i].i, (float)Fin[i].r, (float)Fin[i].i, (float)Fin[i].r);
        C_MUL_SSE(T, Fo, Fi);
        _mm_storel_pi((__m64*)&Fout[i], T);
    }
}
