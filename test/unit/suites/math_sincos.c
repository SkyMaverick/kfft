#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#include "CUnit/Basic.h"
#include "util.h"

#include "kfft_custom_math.h"

#define TEST_CYCLE_COUNT 101
#define TEST_ACCURACY 0.001

_TEST(null) {
    double co, si;
    co = si = 0.0;

    int ret = kfft_sincos_double(NULL, NULL, 0);
    CU_ASSERT_EQUAL(ret, KFFT_MATH_ISNAN);
    ret = kfft_sincos_double(&co, NULL, 0);
    CU_ASSERT_EQUAL(ret, KFFT_MATH_ISNAN);
    ret = kfft_sincos_double(NULL, &si, 0);
    CU_ASSERT_EQUAL(ret, KFFT_MATH_ISNAN);
    ret = kfft_sincos_double(&co, &si, 0);
    CU_ASSERT_EQUAL(ret, KFFT_MATH_NORMAL);
}
_TEST(zero) {
    double co, si;
    co = si = 0.;

    int ret = kfft_sincos_double(&co, &si, 0);

    CU_ASSERT_EQUAL(ret, KFFT_MATH_NORMAL);
    CU_ASSERT_EQUAL(co, 1);
    CU_ASSERT_EQUAL(si, 0);
}
_TEST(normal) {
    double co, si;

    srand(time(NULL));
    for (unsigned i = 0; i < TEST_CYCLE_COUNT; i++) {
        co = si = 0.0;
        unsigned gnm = rand();

        unsigned ret = kfft_sincos_double(&co, &si, gnm);
        if (ret != KFFT_MATH_NORMAL)
            CU_FAIL_FATAL("Return not KFFT_MATH_NORMAL");
        if (((co - cos(gnm)) > TEST_ACCURACY) || ((si - sin(gnm)) > TEST_ACCURACY))
            CU_FAIL_FATAL("Return value not accuracy");
    }
}
_TEST(normal_inverse) {
    double co, si;

    srand(time(NULL));
    for (unsigned i = 0; i < TEST_CYCLE_COUNT; i++) {
        co = si = 0.0;
        int gnm = -1 * rand();

        unsigned ret = kfft_sincos_double(&co, &si, gnm);
        if (ret != KFFT_MATH_NORMAL)
            CU_FAIL_FATAL("Return not KFFT_MATH_NORMAL");
        if (((co - cos(gnm)) > TEST_ACCURACY) || ((si - sin(gnm)) > TEST_ACCURACY))
            CU_FAIL_FATAL("Return value not accuracy");
    }
}
_TEST(hight) {
    double co, si;
    co = si = 0.0;
    double gnm = UINT32_MAX;

    unsigned ret = kfft_sincos_double(&co, &si, gnm);
    if (ret != KFFT_MATH_NORMAL)
        CU_FAIL_FATAL("Return not KFFT_MATH_NORMAL");
    if (((co - cos(gnm)) > TEST_ACCURACY) || ((si - sin(gnm)) > TEST_ACCURACY))
        CU_FAIL_FATAL("Return value not accuracy");
}
_TEST(low) {
    double co, si;
    co = si = 0.0;
    double gnm = -1 * UINT32_MAX;

    unsigned ret = kfft_sincos_double(&co, &si, gnm);
    if (ret != KFFT_MATH_NORMAL)
        CU_FAIL_FATAL("Return not KFFT_MATH_NORMAL");
    if (((co - cos(gnm)) > TEST_ACCURACY) || ((si - sin(gnm)) > TEST_ACCURACY))
        CU_FAIL_FATAL("Return value not accuracy");
}
void
_run_suite(void) {
    CU_pSuite suite = CUnitCreateSuite("math custom sincos() func (double and float)");
    if (suite) {
        _ADD_TEST(suite, null);
        _ADD_TEST(suite, zero);
        _ADD_TEST(suite, normal);
        _ADD_TEST(suite, normal_inverse);
        _ADD_TEST(suite, hight);
        _ADD_TEST(suite, low);
    }
    return;
}
