#pragma once

#include "kfft_types.h"

enum {
    KFFT_MATH_NORMAL,
    KFFT_MATH_OVERFLOW,
    KFFT_MATH_UNDERFLOW,
    KFFT_MATH_ISNAN,
    KFFT_MATH_TLOSS
};

unsigned
kfft_sincos_double(double* co, double* si, double x);

static inline kfft_scalar
kfft_math_sqrt(const kfft_scalar number) {
    if (number <= 0)
        return 0;

    const kfft_scalar ACCURACY = 0.001;
    kfft_scalar lower, upper, guess;

    if (number < 1) {
        lower = number;
        upper = 1;
    } else {
        lower = 1;
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
