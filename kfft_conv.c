#include "kfft.h"
#include "kfft_math.h"

static inline void
kfft_mul(kfft_cpx* Fout, kfft_cpx* Fin, uint32_t size) {
    // FIXME Now primitive algorithm
    kfft_cpx tmp;
    for (uint32_t i = 0; i < size; i++) {
        C_MUL(tmp, Fout[i], Fin[i]);
        C_CPY(Fout[i], tmp);
    }
}

/*
    ALGORITHM
        - FFT Fout sequence with P plan
        - FFT Fin sequence with P plan
        - Per-elements multiplify to Fout
        - Inverse FFT Fout sequence
 */
static inline int
kfft_convolution(kfft_cpx* Fout, kfft_cpx* Fin, kfft_comp_t* P, kfft_comp_t* Pi) {
    kfft_eval_cpx(P, Fin, Fin);
    kfft_eval_cpx(P, Fout, Fout);
    kfft_mul(Fout, Fin, P->nfft);

    kfft_eval_cpx(Pi, Fout, Fout);
}
