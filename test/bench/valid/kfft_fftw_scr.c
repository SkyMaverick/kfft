#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <fftw3.h>
#include "kfft.h"

enum { RETURN_EQUAL = 0, RETURN_NONEQUAL = 1, RETURN_MMFAIL = 2, RETURN_ARGFAIL = 3 };

#if defined(TEST_HALF_SCALAR)
    #define FFTW(X) fftwf_##X
    #define fft_scalar float
#else
    #define FFTW(X) fftw_##X
    #define fft_scalar double
#endif

#define TEST_PRECISION 0.1
#define AMP_LIMIT 1000

FFTW(complex) * fftw_out;
kfft_cpx* kfft_out;

fft_scalar *fftw_in, *kfft_in;

#define TRACE_FAULT_CPX(P, F, K)                                                                   \
    fprintf(stderr,                                                                                \
            "Fault on %zu position:\nFFTW:\n\treal\t%f\n\timg\t%f\nKFFT:\n\treal\t%f\n\timg "      \
            "\t%f\nDiff:\n\treal\t%f\n\timg\t%f\n",                                                \
            (P), (F)[(P)][0], (F)[(P)][1], (K)[(P)].r, (K)[(P)].i, (F)[(P)][0] - (K)[(P)].r,       \
            (F)[(P)][1] - (K)[(P)].i)

#define TRACE_FAULT_SCALAR(P, F, K)                                                                \
    fprintf(stderr,                                                                                \
            "Fault on %zu position:\nFFTW:\n\treal\t%f\nKFFT:\n\treal\t%f\nDiff:\n\treal\t%f\n",   \
            (P), (F)[(P)], (K)[(P)], (F)[(P)] - (K)[(P)])

static unsigned
compare_spectr_inv(size_t size) {
    unsigned ret = RETURN_MMFAIL;

    FFTW(plan)
    fftw_plan = FFTW(plan_dft_c2r_1d)(size, fftw_out, fftw_in, FFTW_ESTIMATE);
    if (fftw_plan) {
        kfft_plan_sclr* kfft_plan =
            kfft_config_scalar(size, KFFT_FLAG_INVERSE | KFFT_FLAG_DISABLE_NORM, NULL, NULL);
        if (kfft_plan) {
            // processing FFTW
            FFTW(execute)(fftw_plan);
            // processing KFFT
            kfft_evali_scalar(kfft_plan, kfft_out, kfft_in);

            // COMPARE SEQUENCES
            ret = RETURN_EQUAL;

            for (size_t i = 0; i < size; i++) {
                if (KFFT_FABS(fftw_in[i] - kfft_in[i]) > TEST_PRECISION) {
                    ret = RETURN_NONEQUAL;
                    TRACE_FAULT_SCALAR(i, fftw_in, kfft_in);
                    goto bailout;
                }
            }

        bailout:
            kfft_cleanup(kfft_plan);
        }
        FFTW(destroy_plan)(fftw_plan);
    }
    return ret;
}

static unsigned
compare_spectr_fwd(size_t size) {
    unsigned ret = RETURN_MMFAIL;

    FFTW(plan)
    fftw_plan = FFTW(plan_dft_r2c_1d)(size, fftw_in, fftw_out, FFTW_ESTIMATE);
    if (fftw_plan) {
        kfft_plan_sclr* kfft_plan = kfft_config_scalar(size, KFFT_FLAG_NORMAL, NULL, NULL);
        if (kfft_plan) {
            // processing FFTW
            FFTW(execute)(fftw_plan);
            // processing KFFT
            kfft_eval_scalar(kfft_plan, kfft_in, kfft_out);

            // COMPARE SEQUENCES
            ret = RETURN_EQUAL;
            for (size_t i = 0; i < size / 2; i++) {
                if (KFFT_FABS(fftw_out[i][0] - kfft_out[i].r) > TEST_PRECISION) {
                    ret = RETURN_NONEQUAL;
                    TRACE_FAULT_CPX(i, fftw_out, kfft_out);
                    goto bailout;
                }
                if (KFFT_FABS(fftw_out[i][1] - kfft_out[i].i) > TEST_PRECISION) {
                    ret = RETURN_NONEQUAL;
                    TRACE_FAULT_CPX(i, fftw_out, kfft_out);
                    goto bailout;
                }
            }

        bailout:
            kfft_cleanup(kfft_plan);
        }
        FFTW(destroy_plan)(fftw_plan);
    }
    return ret;
}

int
main(int argc, char* argv[]) {
    unsigned ret = RETURN_NONEQUAL;
    if (argc > 1) {
        size_t size = atol(argv[1]);

        fftw_in = FFTW(malloc)(size * sizeof(fft_scalar));
        fftw_out = FFTW(malloc)(size * sizeof(FFTW(complex)));

        if (fftw_in && fftw_out) {

            kfft_in = kfft_malloc(size * sizeof(kfft_scalar));
            kfft_out = kfft_malloc(size * sizeof(kfft_cpx));

            if (kfft_in && kfft_out) {
                for (size_t i = 0; i < size; i++) {
                    kfft_in[i] = (fft_scalar)(rand() % (AMP_LIMIT + 1)) /* i+1 */;
                    fftw_in[i] = kfft_in[i];
                }
                ret = compare_spectr_fwd(size);
                if (ret == RETURN_EQUAL) {
                    memset(fftw_in, 0, size * sizeof(fft_scalar));
                    memset(kfft_in, 0, size * sizeof(fft_scalar));
                    ret = compare_spectr_inv(size);
                }

                kfft_free(kfft_in);
                kfft_free(kfft_out);
            }

            FFTW(free)(fftw_in);
            FFTW(free)(fftw_out);
        }

        return ret;
    } else {
        fprintf(stderr, "%s\n", "Need sequense size as parameter");
        return RETURN_ARGFAIL;
    }
}
