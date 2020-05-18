/* Sine/cosine/exponent using the special tables */
/* double precision implementation */

/* functions return */
#define M_CONST_PI 3.141592653589793238462643383279502884197169399375105820974944

#include "kfft_custom_math.h"
#include "kfft_scd.h"

/* decompose a double precision number into 5 unsigned short numbers + sign.
   Two numbers: 2: integer part + 3 numbers fractal part + sign
*/
static int
int_fract(int* sign, unsigned short* dest, double src) {
    unsigned long intp;
    unsigned short *dst, as;
    double a, *fr;
    int i, k;
    if (src < 0.) {
        *sign = -1;
        src = -src;
    } else if (src > 0.)
        *sign = 1;
    else if (src == 0.) {
        *sign = 0;
        for (i = 0; i < 5; i++)
            dest[i] = 0;
        return (KFFT_MATH_NORMAL);
    } else
        return (KFFT_MATH_ISNAN);
    if (src > (double)0xffffffff)
        return (*sign == 1 ? KFFT_MATH_OVERFLOW : KFFT_MATH_UNDERFLOW);
    intp = (unsigned long)src;
    src -= (double)intp;
    dest[0] = (unsigned short)(intp & 0xffff);
    dest[1] = (unsigned short)(intp >> 16);
    dst = dest + 2;
    fr = Frac;
    for (i = 0; i < 3; i++) {
        *dst = 0;
        as = 1;
        for (k = 0; k < 16; k++) {
            if ((a = src - (*fr)) >= 0.) {
                src = a;
                *dst |= as;
            }
            fr++;
            as <<= 1;
        }
        dst++;
    }
    return (KFFT_MATH_NORMAL);
}

/* calculate the sine and cosine */
#define M_2PI (2. * M_CONST_PI)
unsigned
kfft_sincos_double(double* co, double* si, double x) {
    int sign;
    unsigned short xx[5], *xt;
    int i, k, as;
    double st, ct, *sit, *cot;
    /* divide x by 2*pi */
    x /= M_2PI;
    /* decompose x */
    switch (int_fract(&sign, xx, x)) {
    case KFFT_MATH_OVERFLOW:
    case KFFT_MATH_UNDERFLOW: {
        *co = 1.;
        (*si) = 0.;
        return (KFFT_MATH_TLOSS);
    }
    case KFFT_MATH_ISNAN:
        return (KFFT_MATH_ISNAN);
    }
    *co = 1.;
    *si = 0.;
    cot = Cosi;
    sit = Sine;
    xt = xx + 2;
    for (i = 0; i < 3; i++) {
        as = 1;
        for (k = 0; k < 16; k++) {
            if (*xt & as) {
                ct = (*co);
                st = (*si);
                *co = ct * (*cot) - st * (*sit);
                *si = ct * (*sit) + st * (*cot);
            }
            as <<= 1;
            sit++;
            cot++;
        }
        xt++;
    }
    /* if the sign was changed, change sine result */
    if (sign < 0)
        *si = -(*si);
    return (KFFT_MATH_NORMAL);
}
