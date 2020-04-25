#pragma once

#include "kfft_math_sse.h"

#ifdef KFFT_USE_SYSMATH
    #include <math.h>
#endif
#define KFFT_CONST_PI 3.141592653589793238462643383279502884197169399375105820974944

/*
  Explanation of macros dealing with complex math:
   C_MUL(m,a,b)         : m = a*b
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */
#if defined(KFFT_HALF_SCALAR)
// TODO
#else
// TODO
#endif
