#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "CUnit/Basic.h"

#include "defs_cpx1.h"

#define BUFSZ(X) (X) * sizeof(kfft_cpx)
#define TEST_PRECISION 0.001

_TEST(create_std) {}
_TEST(release_std) {}
_TEST(check_size) {}
_TEST(insert) {}
_TEST(insert_underflow) {}
_TEST(insert_overflow) {}
_TEST(renew) {}

#define FFT_FORWARD_ALGO(X)                                                                        \
    do {                                                                                           \
        STDIO_OFF(helper);                                                                         \
                                                                                                   \
        kfft_plan_cpx* plan = NULL;                                                                \
        kfft_cpx* fft_out = NULL;                                                                  \
                                                                                                   \
        if ((plan = kfft_config_cpx(FFT_LENGHT##X, KFFT_FLAG_NORMAL, NULL, NULL)) != NULL) {       \
            if ((fft_out = kfft_malloc(BUFSZ(FFT_LENGHT##X))) != NULL) {                           \
                if (kfft_eval_cpx(plan, fft##X##_in, fft_out) == KFFT_RET_SUCCESS) {               \
                    for (uint32_t i = 0; i < FFT_LENGHT##X; i++)                                   \
                        if ((KFFT_FABS(fft##X##_out[i].r - fft_out[i].r) > TEST_PRECISION) ||      \
                            (KFFT_FABS(fft##X##_out[i].i - fft_out[i].i) > TEST_PRECISION)) {      \
                            CU_FAIL("Precision error");                                            \
                            goto bailout;                                                          \
                        }                                                                          \
                } else                                                                             \
                    CU_FAIL("Complex evaluation fail");                                            \
            bailout:                                                                               \
                kfft_free(fft_out); /* free buffer */                                              \
            } else                                                                                 \
                CU_FAIL("Memory allocation fail");                                                 \
            kfft_cleanup(plan); /* release plan */                                                 \
        } else                                                                                     \
            CU_FAIL("Don't create plan");                                                          \
                                                                                                   \
        STDIO_ON(helper);                                                                          \
    } while (0)

_TEST(simple_valid_2) { FFT_FORWARD_ALGO(2); }
_TEST(simple_valid_3) { FFT_FORWARD_ALGO(3); }
_TEST(simple_valid_5) { FFT_FORWARD_ALGO(5); }
_TEST(simple_valid_prime) { FFT_FORWARD_ALGO(7); }
_TEST(simple_valid_complex) { FFT_FORWARD_ALGO(330); }

#define FFT_INVERSE_ALGO(X)                                                                        \
    do {                                                                                           \
        STDIO_OFF(helper);                                                                         \
                                                                                                   \
        kfft_plan_cpx* plan = NULL;                                                                \
        kfft_cpx* fft_out = NULL;                                                                  \
                                                                                                   \
        if ((plan = kfft_config_cpx(FFT_LENGHT##X, KFFT_FLAG_INVERSE, NULL, NULL)) != NULL) {      \
            if ((fft_out = kfft_malloc(BUFSZ(FFT_LENGHT##X))) != NULL) {                           \
                if (kfft_eval_cpx(plan, fft##X##_out, fft_out) == KFFT_RET_SUCCESS) {              \
                    for (uint32_t i = 0; i < FFT_LENGHT##X; i++)                                   \
                        if ((KFFT_FABS(fft##X##_in[i].r - fft_out[i].r) > TEST_PRECISION) ||       \
                            (KFFT_FABS(fft##X##_in[i].i - fft_out[i].i) > TEST_PRECISION)) {       \
                            CU_FAIL("Precision error");                                            \
                            goto bailout;                                                          \
                        }                                                                          \
                } else                                                                             \
                    CU_FAIL("Complex evaluation fail");                                            \
            bailout:                                                                               \
                kfft_free(fft_out); /* free buffer */                                              \
            } else                                                                                 \
                CU_FAIL("Memory allocation fail");                                                 \
            kfft_cleanup(plan); /* release plan */                                                 \
        } else                                                                                     \
            CU_FAIL("Don't create plan");                                                          \
                                                                                                   \
        STDIO_ON(helper);                                                                          \
    } while (0)

_TEST(simple_inverse_valid_2) { FFT_INVERSE_ALGO(2); }
_TEST(simple_inverse_valid_3) { FFT_INVERSE_ALGO(3); }
_TEST(simple_inverse_valid_5) { FFT_INVERSE_ALGO(5); }
_TEST(simple_inverse_valid_prime) { FFT_INVERSE_ALGO(7); }
_TEST(simple_inverse_valid_complex) { FFT_INVERSE_ALGO(330); }

void
_run_suite(void) {
    // ... all startup functionality
    CU_pSuite suite = CUnitCreateSuite("1d complex API");
    if (suite) {
        _ADD_TEST(suite, create_std);
        _ADD_TEST(suite, release_std);
        _ADD_TEST(suite, check_size);
        _ADD_TEST(suite, insert);
        _ADD_TEST(suite, insert_underflow);
        _ADD_TEST(suite, insert_overflow);
        _ADD_TEST(suite, renew);

        _ADD_TEST(suite, simple_valid_2);
        _ADD_TEST(suite, simple_valid_3);
        _ADD_TEST(suite, simple_valid_5);
        _ADD_TEST(suite, simple_valid_prime);
        _ADD_TEST(suite, simple_valid_complex);
        _ADD_TEST(suite, simple_inverse_valid_2);
        _ADD_TEST(suite, simple_inverse_valid_3);
        _ADD_TEST(suite, simple_inverse_valid_5);
        _ADD_TEST(suite, simple_inverse_valid_prime);
        _ADD_TEST(suite, simple_inverse_valid_complex);
    }
    return;
}
