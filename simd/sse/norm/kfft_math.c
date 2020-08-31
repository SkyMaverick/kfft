#include "kfft_simd_compat.h"
#include "kfft_simd.h"

#include "kfft_math_intern.h"

void
FUNC_SSE(kfft_math_hadamard_cpx)(kfft_cpx* Fout, const kfft_cpx* Fin, uint32_t size) {
    // FIXME Now primitive algorithm
    //    kfft_math_adamar_cpx(Fout, Fin, size);
    KFFT_OMP( omp parallel for schedule(static) private(tmp))
    for (uint32_t i = 0; i < size; i++) {
        __m128d Fo, Fi, T;
        Fo = _mm_load_pd((double*)&Fout[i]);
        /* FIXME if use _mm_load_pd debug build crash with SIGSEGV */
        Fi = _mm_set_pd((double)Fin[i].i, (double)Fin[i].r);
        C_MUL_SSE(T, Fo, Fi);
        _mm_store_pd((double*)&Fout[i], T);
    }
}
