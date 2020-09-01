#pragma once

/*!
    \file

    Custom implementation trivial mathematical operations
    without libm dependency
 */

#include "kfft_types.h"

/// SinCos function return codes
enum kfft_sincos_ret {
    KFFT_MATH_NORMAL,    ///< normal evaluate
    KFFT_MATH_OVERFLOW,  ///< overflow (src > 0)
    KFFT_MATH_UNDERFLOW, ///< overflow (src < 0)
    KFFT_MATH_ISNAN,     ///< no source or dest argument
    KFFT_MATH_TLOSS      ///< loss evaluation params
};

/*!
    Calculate sinus and cosinus value for custom phase.

    \param[in] co - cosinus buffer pointer (return)
    \param[in] si - sinus buffer pointer (return)
    \param[in] x - phase argument

    \result Return code ::kfft_sincos_ret
 */
unsigned
kfft_sincos_double(double* co, double* si, double x);

/*!
    Found square root for number.
    \note This function use simpliest bisection method.
    It's not faster, but simple.

    \param[in] number - aargument
    \result Found square root
 */

static inline kfft_scalar
kfft_math_sqrt(const kfft_scalar number) {
    if (number <= 0)
        return 0;

    const kfft_scalar ACCURACY = 0.001;
    kfft_scalar lower, upper, guess;

    lower = upper = 1;

    if (number < 1) {
        lower = number;
    } else {
        upper = number;
    }

    while ((upper - lower) > ACCURACY) {
        guess = (lower + upper) / 2;
        if (guess * guess > number)
            upper = guess;
        else
            lower = guess;
    }
    return (lower + upper) / 2;
}
