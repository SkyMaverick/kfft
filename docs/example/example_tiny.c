#include "kfft_tiny.h"
#include <stdio.h>

#define KFFT_ARRAY_SIZE 100

#define KFFT_ARRAY_X 50
#define KFFT_ARRAY_Y 50
#define KFFT_ARRAY_SIZE2 ((KFFT_ARRAY_X) * (KFFT_ARRAY_Y))

/* Complex 2D convolution example */
static kfft_return_t
kfft_complex_conv_2d(void) {
    kfft_return_t ret = KFFT_RET_ALLOC_FAIL;
    kfft_cpx *in1, *in2, *out;

    in1 = in2 = out = kfft_malloc(3 * KFFT_ARRAY_SIZE2 * sizeof(kfft_cpx));
    if (out) {
        in1 += KFFT_ARRAY_SIZE2;
        in2 += 2 * KFFT_ARRAY_SIZE2;

        /* Fill complex data to IN buffer HERE !!! */

        kfft_args_2d A = {KFFT_ARRAY_X, KFFT_ARRAY_Y};
        kfft_plan plan = kfft_configuration(KFFT_PLAN_COMPLEX_CONV2D, KFFT_FLAG_GENERIC, &A);
        if (plan) {
            ret = kfft_convolution(plan, in1, in2, out);

            /* STDOUT complex data contains OUT buffer HERE !!! */

            kfft_finalization(plan);
        }
        kfft_free(out);
    }
    return ret;
}

/* Scalar 2D convolution example */
static kfft_return_t
kfft_scalar_conv_2d(void) {
    kfft_return_t ret = KFFT_RET_ALLOC_FAIL;
    kfft_cpx *in1, *in2, *out;
    /*
       It is allowed to use one buffer,
       but kfft_evaluation() will allocate an equivalent temporary buffer.
     */
    in1 = in2 = out = kfft_malloc(3 * KFFT_ARRAY_SIZE2 * sizeof(kfft_scalar));
    if (out) {
        in1 += KFFT_ARRAY_SIZE2;
        in2 += 2 * KFFT_ARRAY_SIZE2;

        /* Fill complex data to IN buffer HERE !!! */

        kfft_args_2d A = {KFFT_ARRAY_X, KFFT_ARRAY_Y};
        kfft_plan plan = kfft_configuration(KFFT_PLAN_SCALAR_CONV2D, KFFT_FLAG_NORMAL, &A);
        if (plan) {
            ret = kfft_convolution(plan, in1, in2, out);

            /* STDOUT complex data contains OUT buffer HERE !!! */

            kfft_finalization(plan);
        }
        kfft_free(out);
    }
    return ret;
}

/* Scalar 1D buffer evaluate example */
static kfft_return_t
kfft_scalar_1d(void) {
    kfft_return_t ret = KFFT_RET_ALLOC_FAIL;
    kfft_scalar* in = NULL;
    kfft_cpx* out = NULL;

    in = kfft_malloc(KFFT_ARRAY_SIZE * sizeof(kfft_cpx) + sizeof(kfft_scalar));
    if (in && out) {
        out = (kfft_cpx*)(in + KFFT_ARRAY_SIZE);

        /* Fill complex data to IN buffer HERE !!! */

        size_t A = KFFT_ARRAY_SIZE;
        kfft_plan plan = kfft_configuration(KFFT_PLAN_SCALAR, KFFT_FLAG_NORMAL, &A);
        if (plan) {
            ret = kfft_evaluation(plan, in, out);

            /* STDOUT complex data contains OUT buffer HERE !!! */

            kfft_finalization(plan);
        }
        kfft_free(in);
        kfft_free(out);
    }
    return ret;
}

/* Complex 1D buffer evaluate example */
static kfft_return_t
kfft_complex_1d(void) {
    kfft_return_t ret = KFFT_RET_ALLOC_FAIL;
    kfft_cpx *in, *out;
    /*
       It is allowed to use one buffer,
       but kfft_evaluation() will allocate an equivalent temporary buffer.
     */
    in = out = kfft_malloc(2 * KFFT_ARRAY_SIZE * sizeof(kfft_cpx));
    if (out) {
        in += KFFT_ARRAY_SIZE;

        /* Fill complex data to IN buffer HERE !!! */

        size_t A = KFFT_ARRAY_SIZE;
        kfft_plan plan =
            kfft_configuration(KFFT_PLAN_COMPLEX, KFFT_FLAG_INVERSE | KFFT_FLAG_GENERIC_ONLY, &A);
        if (plan) {
            ret = kfft_evaluation(plan, in, out);

            /* STDOUT complex data contains OUT buffer HERE !!! */

            kfft_finalization(plan);
        }
        kfft_free(out);
    }
    return ret;
}

int
main() {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    if ((ret = kfft_complex_1d()) != KFFT_RET_SUCCESS)
        fprintf(stderr, "%s\n", kfft_strerr(ret));

    if ((ret = kfft_scalar_1d()) != KFFT_RET_SUCCESS)
        fprintf(stderr, "%s\n", kfft_strerr(ret));

    if ((ret = kfft_scalar_conv_2d()) != KFFT_RET_SUCCESS)
        fprintf(stderr, "%s\n", kfft_strerr(ret));

    if ((ret = kfft_complex_conv_2d()) != KFFT_RET_SUCCESS)
        fprintf(stderr, "%s\n", kfft_strerr(ret));
}
