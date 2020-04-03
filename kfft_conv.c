#include "kfft.h"

static inline kfft_return_t
kfft_part_convolution(kfft_cpx* Fout, kfft_cpx* Fin, kfft_comp_t* P, kfft_comp_t* Pi) {
    kfft_return_t ret = kfft_eval_cpx(P, Fout, Fout);
    if (ret == KFFT_RET_SUCCESS) {
        kfft_math_adamar_cpx(Fout, Fin, P->nfft);
        ret = kfft_eval_cpx(Pi, Fout, Fout);
    }
    return ret;
}

/*
    ALGORITHM
        - FFT Fout sequence with P plan
        - FFT Fin sequence with P plan
        - Per-elements multiplify to Fout
        - Inverse FFT Fout sequence
 */
KFFT_API kfft_return_t
kfft_convolution(kfft_cpx* Fout, kfft_cpx* Fin, kfft_comp_t* P, kfft_comp_t* Pi) {

    kfft_return_t ret = kfft_eval_cpx(P, Fin, Fin);
    return (ret == KFFT_RET_SUCCESS) ? kfft_part_convolution(Fout, Fin, P, Pi) : ret;
}
