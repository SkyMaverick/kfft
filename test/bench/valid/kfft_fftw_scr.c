#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <fftw3.h>
#include "kfft.h"

enum { RETURN_EQUAL = 0, RETURN_NONEQUAL = 1, RETURN_MMFAIL = 2, RETURN_ARGFAIL = 3 };

#if defined(TEST_HALF_SCALAR)
    #define FFTW(X) fftwf_##X
    #define fft_scalar float
#else
    #define FFTW(X) fftw_##X
    #define fft_scalar double
#endif

#define TEST_PRECISION 0.001

static unsigned
compare_forward(size_t size) {
    unsigned ret = RETURN_MMFAIL;

    fft_scalar* fftw_in = FFTW(malloc)(size * sizeof(fft_scalar));
    fft_scalar* kfft_in = kfft_malloc(size * sizeof(fft_scalar));

    if (fftw_in && kfft_in) {

        for (uint32_t i = 0; i < size; i++) {
            fftw_in[i] = i;
            kfft_in[i] = i;
        }

        FFTW(complex)* fftw_out = FFTW(malloc)(size * sizeof(FFTW(complex)));
        kfft_cpx* kfft_out = kfft_malloc(size * sizeof(kfft_cpx));

        if (fftw_out && kfft_out) {
            FFTW(plan)
            fftw_plan = FFTW(plan_dft_r2c_1d)(size, fftw_in, fftw_out, FFTW_ESTIMATE);

            if (fftw_plan) {
                kfft_plan_sclr* kfft_plan = kfft_config_scalar(size, KFFT_FLAG_NORMAL, NULL, NULL);

                if (kfft_plan) {
                    // processing FFTW
                    FFTW(execute)(fftw_plan);
                    // processing KFFT
                    kfft_eval_scalar(kfft_plan, kfft_in, kfft_out);

                    // COMPARE SEQUENCES
                    ret = RETURN_EQUAL;
                    for (size_t i = 0; i < size / 2; i++) {
                        //                        fprintf (stdout, "[%zu] Real %5.3f - %5.3f  ", i,
                        //                        fftw_out[i][0], kfft_out[i].r); fprintf (stdout,
                        //                        "Image %5.3f - %5.3f\n", fftw_out[i][1],
                        //                        kfft_out[i].i);
                        fprintf(stdout, "[%zu] Diff %5.3f - %5.3f\n", i,
                                fftw_out[i][0] - kfft_out[i].r, fftw_out[i][1] - kfft_out[i].i);
                        if (fabs(fftw_out[i][0] - kfft_out[i].r) > TEST_PRECISION) {
                            ret = RETURN_NONEQUAL;
                            goto bailout;
                        }
                        if (fabs(fftw_out[i][1] - kfft_out[i].i) > TEST_PRECISION) {
                            ret = RETURN_NONEQUAL;
                            goto bailout;
                        }
                    }

                bailout:
                    kfft_cleanup(kfft_plan);
                }
                FFTW(destroy_plan)(fftw_plan);
            }
            FFTW(free)(fftw_out);
            kfft_free(kfft_out);
        }
        FFTW(free)(fftw_in);
        kfft_free(kfft_in);
    }

    return ret;
}

static unsigned
compare_inverse(size_t size) {
    unsigned ret = RETURN_MMFAIL;

    FFTW(complex)* fftw_in = FFTW(malloc)(size * sizeof(FFTW(complex)));
    kfft_cpx* kfft_in = kfft_malloc(size * sizeof(kfft_cpx));

    if (fftw_in && kfft_in) {

        for (uint32_t i = 0; i < size; i++) {
            fftw_in[i][0] = size - i;
            kfft_in[i].r = size - i;
            fftw_in[i][1] = 0;
            kfft_in[i].i = 0;
        }

        fft_scalar* fftw_out = FFTW(malloc)(size * sizeof(fft_scalar));
        fft_scalar* kfft_out = kfft_malloc(size * sizeof(fft_scalar));

        if (fftw_out && kfft_out) {
            FFTW(plan)
            fftw_plan = FFTW(plan_dft_c2r_1d)(size / 2, fftw_in, fftw_out, FFTW_ESTIMATE);

            if (fftw_plan) {
                kfft_plan_sclr* kfft_plan = kfft_config_scalar(size, KFFT_FLAG_INVERSE, NULL, NULL);

                if (kfft_plan) {
                    // processing FFTW
                    FFTW(execute)(fftw_plan);
                    // processing KFFT
                    kfft_evali_scalar(kfft_plan, kfft_in, kfft_out);

                    // COMPARE SEQUENCES
                    ret = RETURN_EQUAL;
                    for (size_t i = 0; i < size / 2; i++) {
                        fprintf(stdout, "[%zu] Real %5.3f - %5.3f  ", i, fftw_out[i], kfft_out[i]);
                        //                        fprintf (stdout, "Image %5.3f - %5.3f\n",
                        //                        fftw_out[i], kfft_out[i]); fprintf (stdout, "[%zu]
                        //                        Diff %5.3f - %5.3f\n", i, fftw_out[i] -
                        //                        kfft_out[i], fftw_out[i] - kfft_out[i]);
                        if (fabs(fftw_out[i] - kfft_out[i]) > TEST_PRECISION) {
                            ret = RETURN_NONEQUAL;
                            goto bailout;
                        }
                    }

                bailout:
                    kfft_cleanup(kfft_plan);
                }
                FFTW(destroy_plan)(fftw_plan);
            }
            FFTW(free)(fftw_out);
            kfft_free(kfft_out);
        }
        FFTW(free)(fftw_in);
        kfft_free(kfft_in);
    }

    return ret;
}

// static unsigned
// compare_spectr_inv(fft_scalar* fftw_buffer, fft_scalar* kfft_buffer, size_t size) {
//     unsigned ret = RETURN_MMFAIL;
//     fft_scalar* fftw_out = FFTW(malloc)(size * sizeof(fft_scalar));
//     fft_scalar* kfft_out = kfft_malloc(size * sizeof(fft_scalar));
//
//     if (fftw_out && kfft_out) {
//
//         FFTW(plan) fftw_plan = FFTW(plan_dft_c2r_1d)(size, (FFTW(complex)*)fftw_buffer, fftw_out,
//         FFTW_ESTIMATE); if (fftw_plan) {
//             kfft_plan_sclr* kfft_plan = kfft_config_scalar(size, KFFT_FLAG_INVERSE, NULL, NULL);
//             if (kfft_plan) {
//                 // processing FFTW
//                 FFTW(execute)(fftw_plan);
//                 // processing KFFT
//                 kfft_evali_scalar(kfft_plan, (kfft_cpx*)kfft_buffer, kfft_out);
//
//                 ret = RETURN_EQUAL;
//                 // COMPARE SEQUENCES
//                 for (size_t i = 0; i < size; i++) {
//                     if (fabs(fftw_out[i] - kfft_out[i]) > TEST_PRECISION) {
//                         ret = RETURN_NONEQUAL;
//                         goto bailout;
//                     }
//                 }
//
//             bailout:
//                 kfft_cleanup(kfft_plan);
//             }
//             FFTW(destroy_plan)(fftw_plan);
//         }
//         FFTW(free)(fftw_out);
//         kfft_free(kfft_out);
//     }
//     return ret;
// }
//
// static unsigned
// compare_spectr_fwd(fft_scalar* fftw_buffer, fft_scalar* kfft_buffer, size_t size) {
//     unsigned ret = RETURN_MMFAIL;
//
//     FFTW(complex)* fftw_out = FFTW(malloc)(size * sizeof(FFTW(complex)));
//     kfft_cpx* kfft_out = kfft_malloc(size * sizeof(kfft_cpx));
//
//     if (fftw_out && kfft_out) {
//
//         FFTW(plan) fftw_plan = FFTW(plan_dft_r2c_1d)(size, fftw_buffer, fftw_out, FFTW_ESTIMATE);
//         if (fftw_plan) {
//             kfft_plan_sclr* kfft_plan = kfft_config_scalar(size, KFFT_FLAG_NORMAL, NULL, NULL);
//             if (kfft_plan) {
//                 // processing FFTW
//                 FFTW(execute)(fftw_plan);
//                 // processing KFFT
//                 kfft_eval_scalar(kfft_plan, kfft_buffer, kfft_out);
//
//                 // COMPARE SEQUENCES
//
//                 ret = RETURN_EQUAL;
//                 for (size_t i = 0; i < size; i++) {
//                     printf("BOOO\n");
//                     if (fabs(fftw_out[i][0] - kfft_out[i].r) > TEST_PRECISION) {
//                         ret = RETURN_NONEQUAL;
//                         goto bailout;
//                     }
//                     if (fabs(fftw_out[i][1] - kfft_out[i].i) > TEST_PRECISION) {
//                         ret = RETURN_NONEQUAL;
//                         goto bailout;
//                     }
//                 }
//
//             bailout:
//                 kfft_cleanup(kfft_plan);
//             }
//             FFTW(destroy_plan)(fftw_plan);
//         }
//         FFTW(free)(fftw_out);
//         kfft_free(kfft_out);
//     }
//
//     return ret;
// }

int
main(int argc, char* argv[]) {
    unsigned ret = RETURN_NONEQUAL;
    if (argc > 1) {
        size_t size = atol(argv[1]);

        ret = compare_forward(size);
        if (ret == RETURN_EQUAL)
            ret = compare_inverse(size);
        //        fft_scalar* fftw_buffer = FFTW(malloc)(size * sizeof(fft_scalar));
        //        if (fftw_buffer) {
        //            fft_scalar* kfft_buffer = kfft_malloc(size * sizeof(fft_scalar));
        //            if (kfft_buffer) {
        //                for (size_t i = 0; i < size; i++) {
        //                    kfft_buffer[i] = (fft_scalar)i;
        //                    fftw_buffer[i] = kfft_buffer[i];
        //                }
        //
        //                ret = compare_spectr_fwd(fftw_buffer, kfft_buffer, size);
        ////                if (ret == RETURN_EQUAL)
        ////                    ret = compare_spectr_inv(fftw_buffer, kfft_buffer, size);
        //                kfft_free(kfft_buffer);
        //            }
        //            FFTW(free)(fftw_buffer);
        //        }
    } else {
        fprintf(stderr, "%s\n", "Need sequense size as parameter");
        return RETURN_ARGFAIL;
    }
    return ret;
}
