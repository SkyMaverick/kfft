#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "CUnit/Basic.h"

_TEST(foo) {
//  CU_ASSERT_EQUAL(0,1);
}

_TEST(foo1) {
//  CU_ASSERT_EQUAL(1,1);
}

void 
_run_suite (void) {
    // ... all startup functionality
    CU_pSuite suite=CUnitCreateSuite("SuiteName");
    if(suite) {
        _ADD_TEST(suite, foo);
        _ADD_TEST(suite, foo1);
    }
    return;
}
