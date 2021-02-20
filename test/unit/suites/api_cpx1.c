#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "CUnit/Basic.h"

#include "defs_cpx1.h"

#define BUFSZ(X) (X) * sizeof(kfft_cpx)
#define TEST_PRECISION 0.001

_TEST(create_std) {
    STDIO_OFF(helper);

    kfft_plan_cpx* plan = kfft_config_cpx(1, KFFT_FLAG_NORMAL, NULL, NULL);

    CU_ASSERT_PTR_NOT_NULL(plan);
    if (plan)
        kfft_release(plan);

    CU_ASSERT_PTR_NULL(plan);

    STDIO_ON(helper);
}

_TEST(check_size_zero) {
    STDIO_OFF(helper);

    size_t check = 0;
    kfft_plan_cpx* plan = kfft_config_cpx(0, KFFT_FLAG_NORMAL, NULL, &check);
    if (plan)
        kfft_release(plan);
    CU_ASSERT_EQUAL(check, sizeof(kfft_plan_cpx));

    STDIO_ON(helper);
}

_TEST(check_size_norm) {
    STDIO_OFF(helper);

    size_t check = 0;
    kfft_plan_cpx* plan = kfft_config_cpx(101, KFFT_FLAG_NORMAL, NULL, &check);
    if (plan)
        kfft_release(plan);
    CU_ASSERT_NOT_EQUAL(KFFT_FABS(check - sizeof(kfft_plan_cpx)), 0);

    STDIO_ON(helper);
}

_TEST(insert) {
    STDIO_OFF(helper);

    size_t szBlock_1, szBlock_2;
    szBlock_1 = szBlock_2 = 0;

    kfft_config_cpx(101, KFFT_FLAG_INVERSE, NULL, &szBlock_1);
    kfft_config_cpx(78, KFFT_FLAG_NORMAL, NULL, &szBlock_2);

    if (szBlock_1 && szBlock_2) {
        kfft_pool_t* P = kfft_pool_create(szBlock_1 + szBlock_2);
        if (P) {
            kfft_plan_cpx* kpcpx1 = kfft_config_cpx(78, KFFT_FLAG_NORMAL, P, NULL);
            kfft_plan_cpx* kpcpx2 = kfft_config_cpx(101, KFFT_FLAG_INVERSE, P, NULL);

            CU_ASSERT_PTR_NOT_NULL(kpcpx1);
            CU_ASSERT_PTR_NOT_NULL(kpcpx2);

            kfft_pool_free(P);
        } else
            CU_FAIL("Allocator error");
    } else
        CU_FAIL("Calculate error");

    STDIO_ON(helper);
}

_TEST(insert_overflow) {
    STDIO_OFF(helper);

    size_t szBlock_1, szBlock_2;
    szBlock_1 = szBlock_2 = 0;

    kfft_config_cpx(101, KFFT_FLAG_INVERSE, NULL, &szBlock_1);
    kfft_config_cpx(78, KFFT_FLAG_NORMAL, NULL, &szBlock_2);

    if (szBlock_1 && szBlock_2) {
        kfft_pool_t* P = kfft_pool_create(szBlock_1 + szBlock_2);
        if (P) {
            kfft_plan_cpx* kpcpx1 = kfft_config_cpx(78, KFFT_FLAG_NORMAL, P, NULL);
            kfft_plan_cpx* kpcpx2 = kfft_config_cpx(107, KFFT_FLAG_INVERSE, P, NULL);

            CU_ASSERT_PTR_NOT_NULL(kpcpx1);
            CU_ASSERT_PTR_NULL(kpcpx2);

            kfft_pool_free(P);
        } else
            CU_FAIL("Allocator error");
    } else
        CU_FAIL("Calculate error");

    STDIO_ON(helper);
}

_TEST(renew_simple) {
    STDIO_OFF(handler);

    kfft_plan_cpx *plan_inv = NULL, *plan = NULL;
    kfft_cpx *fft_out = NULL, *fft_in = NULL;

    if ((plan = kfft_config_cpx(FFT_LENGHT330, KFFT_FLAG_NORMAL, NULL, NULL)) != NULL) {
        if ((fft_out = kfft_malloc(BUFSZ(FFT_LENGHT330))) != NULL) {
            if (kfft_eval_cpx(plan, fft330_in, fft_out) == KFFT_RET_SUCCESS) {
                for (uint32_t i = 0; i < FFT_LENGHT330; i++)
                    if ((KFFT_FABS(fft330_out[i].r - fft_out[i].r) > TEST_PRECISION) ||
                        (KFFT_FABS(fft330_out[i].i - fft_out[i].i) > TEST_PRECISION)) {
                        CU_FAIL("Precision error");
                        goto cleanup;
                    }
            } else {
                CU_FAIL("Complex forward evaluation fail");
                goto cleanup;
            }
        } else
            CU_FAIL("Memory forward allocation fail");

        // RENEW - new inverse plan
        if ((plan_inv = kfft_config_cpx(FFT_LENGHT330, KFFT_FLAG_INVERSE | KFFT_FLAG_RENEW,
                                        KFFT_PLAN_MMGR(plan), NULL)) != NULL) {
            if ((fft_in = kfft_malloc(BUFSZ(FFT_LENGHT330))) != NULL) {
                if (kfft_eval_cpx(plan, fft_out, fft_in) == KFFT_RET_SUCCESS) {
                    for (uint32_t i = 0; i < FFT_LENGHT330; i++)
                        if ((KFFT_FABS(fft_in[i].r - fft_in[i].r) > TEST_PRECISION) ||
                            (KFFT_FABS(fft_in[i].i - fft_in[i].i) > TEST_PRECISION)) {
                            CU_FAIL("Precision error");
                            goto cleanup;
                        }
                } else {
                    CU_FAIL("Complex inverse evaluation fail");
                    goto cleanup;
                }
            } else
                CU_FAIL("Memory inverse allocation fail");
        } else
            CU_FAIL("Don't create plan");
    cleanup:
        if (fft_in)
            kfft_free(fft_in);
        if (fft_out)
            kfft_free(fft_out);
        kfft_release(plan); /* release plan */
    } else
        CU_FAIL("Don't create plan");

    STDIO_ON(handler);
}

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
            kfft_release(plan); /* release plan */                                                 \
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
            kfft_release(plan); /* release plan */                                                 \
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
        _ADD_TEST(suite, create_std);      // Standart create plan
        _ADD_TEST(suite, check_size_zero); // Get size zero-lenght plan
        _ADD_TEST(suite, check_size_norm); // Get size normal situation
        _ADD_TEST(suite, insert);          // Insert more plans in pool
        _ADD_TEST(suite, insert_overflow); // Insert overflow plan in pool
        _ADD_TEST(suite, renew_simple);    // Renew for inverse plan (any theme)

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
