#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "CUnit/Basic.h"
#include "util.h"

#include "kfft_custom_math.h"

#define TEST_CYCLE_COUNT 13
#define TEST_ACCURACY 0.001

_TEST(zero) { CU_ASSERT_EQUAL(0, kfft_math_sqrt(0)); }
_TEST(negative) {
    double tsum = 0;

    srand(time(NULL));
    for (unsigned i = 0; i < TEST_CYCLE_COUNT; i++)
        tsum += fabs(kfft_math_sqrt(rand() * -1));

    CU_ASSERT_EQUAL(0, tsum);
}
_TEST(normal) {
    srand(time(NULL));
    for (unsigned i = 0; i < TEST_CYCLE_COUNT; i++) {
        unsigned gnm = rand();
        if (kfft_math_sqrt(gnm) - sqrt(gnm) > TEST_ACCURACY)
            CU_FAIL_FATAL("SQRT accuracy problem");
    }
}

void
_run_suite(void) {
    CU_pSuite suite = CUnitCreateSuite("math custom sqrt() func");
    if (suite) {
        _ADD_TEST(suite, zero);
        _ADD_TEST(suite, negative);
        _ADD_TEST(suite, normal);
    }
    return;
}
