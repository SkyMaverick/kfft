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

#define TEST_PRECISION 0.001

static unsigned
compare_spectr_inv(fft_scalar* fftw_buffer, fft_scalar* kfft_buffer, size_t size) {
    unsigned ret = RETURN_MMFAIL;
    FFTW(complex)* fftw_in = (FFTW(complex)*)(fftw_buffer + size);
    kfft_cpx* kfft_in = (kfft_cpx*)(kfft_buffer + size);

    FFTW(plan) fftw_plan = FFTW(plan_dft_c2r_1d)(size, fftw_in, fftw_buffer, FFTW_ESTIMATE);
    if (fftw_plan) {
        kfft_plan_sclr* kfft_plan = kfft_config_scalar(size, KFFT_FLAG_INVERSE, NULL, NULL);
        if (kfft_plan) {
            // processing FFTW
            FFTW(execute)(fftw_plan);
            // processing KFFT
            kfft_evali_scalar(kfft_plan, kfft_in, kfft_buffer);

            // COMPARE SEQUENCES
            for (size_t i = 0; i < size; i++) {
                if (fabs(fftw_buffer[i] - kfft_buffer[i]) > TEST_PRECISION) {
                    ret = RETURN_NONEQUAL;
                    goto bailout;
                }
            }

        bailout:
            kfft_cleanup(kfft_plan);
        }
        FFTW(destroy_plan)(fftw_plan);
        ret = RETURN_EQUAL;
    }
    return ret;
}

static unsigned
compare_spectr_fwd(fft_scalar* fftw_buffer, fft_scalar* kfft_buffer, size_t size) {
    unsigned ret = RETURN_MMFAIL;
    FFTW(complex)* fftw_out = (FFTW(complex)*)(fftw_buffer + size);
    kfft_cpx* kfft_out = (kfft_cpx*)(kfft_buffer + size);

    FFTW(plan) fftw_plan = FFTW(plan_dft_r2c_1d)(size, fftw_buffer, fftw_out, FFTW_ESTIMATE);
    if (fftw_plan) {
        kfft_plan_sclr* kfft_plan = kfft_config_scalar(size, KFFT_FLAG_NORMAL, NULL, NULL);
        if (kfft_plan) {
            // processing FFTW
            FFTW(execute)(fftw_plan);
            // processing KFFT
            kfft_eval_scalar(kfft_plan, kfft_buffer, kfft_out);

            // COMPARE SEQUENCES
            for (size_t i = 0; i < size; i++) {
                if (fabs(fftw_out[i][0] - kfft_out[i].r) > TEST_PRECISION) {
                    ret = RETURN_NONEQUAL;
                    goto bailout;
                }
                if (fabs(fftw_out[i][1] - kfft_out[i].i) > TEST_PRECISION) {
                    ret = RETURN_NONEQUAL;
                    goto bailout;
                }
            }

        bailout:
            kfft_cleanup(kfft_plan);
        }
        FFTW(destroy_plan)(fftw_plan);
        ret = RETURN_EQUAL;
    }

    return ret;
}

int
main(int argc, char* argv[]) {
    unsigned ret = RETURN_NONEQUAL;
    if (argc > 1) {
        size_t size = atol(argv[1]);

        fft_scalar* fftw_buffer = FFTW(malloc)(size * (sizeof(FFTW(complex)) + sizeof(fft_scalar)));
        if (fftw_buffer) {
            fft_scalar* kfft_buffer =
                kfft_malloc(size * (sizeof(FFTW(complex)) + sizeof(fft_scalar)));
            if (kfft_buffer) {
                for (size_t i = 0; i < size; i++) {
                    kfft_buffer[i] = (fft_scalar)rand();
                    fftw_buffer[i] = kfft_buffer[i];
                }

                ret = compare_spectr_fwd(fftw_buffer, kfft_buffer, size);
                if (ret == RETURN_EQUAL)
                    ret = compare_spectr_inv(fftw_buffer, kfft_buffer, size);
                kfft_free(&kfft_buffer);
            }
            FFTW(free)(fftw_buffer);
        }
    } else {
        fprintf(stderr, "%s\n", "Need sequense size as parameter");
        return RETURN_ARGFAIL;
    }
    return ret;
}
