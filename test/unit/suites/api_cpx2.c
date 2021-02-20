#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "CUnit/Basic.h"

#include "defs_cpx2.h"

#define BUFSZ(X) (X) * sizeof(kfft_cpx)
#define BUFSZ2(X, Y) (X) * (Y) * sizeof(kfft_cpx)

#define TEST_PRECISION 0.001


_TEST(simple_valid_23){
}
_TEST(simple_valid_53){
}
_TEST(simple_valid_79){
}
_TEST(simple_valid_1011){
}
_TEST(simple_inverse_valid_23){
}
_TEST(simple_inverse_valid_53){
}
_TEST(simple_inverse_valid_79){
}
_TEST(simple_inverse_valid_1011){
}

void
_run_suite(void) {
    // ... all startup functionality
    CU_pSuite suite = CUnitCreateSuite("2d complex API");
    if (suite) {
        _ADD_TEST(suite, simple_valid_23);
        _ADD_TEST(suite, simple_valid_53);
        _ADD_TEST(suite, simple_valid_79);
        _ADD_TEST(suite, simple_valid_1011);
        _ADD_TEST(suite, simple_inverse_valid_23);
        _ADD_TEST(suite, simple_inverse_valid_53);
        _ADD_TEST(suite, simple_inverse_valid_79);
        _ADD_TEST(suite, simple_inverse_valid_1011);
    }
    return;
}
