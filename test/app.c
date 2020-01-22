#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "kfft.h"
#if defined FFTW_COMPARE
    #include <fftw3.h>
#endif

#ifndef SIMPLE_APP
    #ifndef TEST_COUNT
        #define TEST_COUNT 50
    #endif
#else
    #define TEST_COUNT 1
#endif

#define VALS_KFFT 0
#define VALS_FFTW 1

static double
kfft_ktest(kfft_scalar* amp_scalar, kfft_cpx* tmp_buffer, size_t size) {
#ifndef CHECK_WITHOUT_PLAN
    clock_t t_start = clock();

    kfft_t FCfg = kfft_config(size, 0, 0, NULL);
    if (kfft_isnull(FCfg)) {
        assert("Allocation FAIL");
    }
    kfft(FCfg, amp_scalar, tmp_buffer);
    kfft_free(&FCfg);

    return (clock() - t_start) * 1000 / CLOCKS_PER_SEC;
#else

    kfft_t FCfg = kfft_config(size, 0, 0, NULL);
    if (kfft_isnull(FCfg)) {
        assert("Allocation FAIL");
    }
    clock_t t_start = clock();
    kfft(FCfg, amp_scalar, tmp_buffer);
    clock_t t_ret = clock();

    kfft_free(&FCfg);

    return (t_ret - t_start) * 1000 / CLOCKS_PER_SEC;
#endif
}

#if defined FFTW_COMPARE
static double
__fftw_test(double* amp_scalar, fftw_complex* tmp_buffer, size_t size) {
    #ifndef CHECK_WITHOUT_PLAN
    clock_t t_start = clock();
    fftw_plan plan = fftw_plan_dft_r2c_1d(size, amp_scalar, tmp_buffer, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return (clock() - t_start) * 1000 / CLOCKS_PER_SEC;
    #else
    fftw_plan plan = fftw_plan_dft_r2c_1d(size, amp_scalar, tmp_buffer, FFTW_ESTIMATE);
    clock_t t_start = clock();
    fftw_execute(plan);
    clock_t t_ret = clock();
    fftw_destroy_plan(plan);
    return (t_ret - t_start) * 1000 / CLOCKS_PER_SEC;
    #endif /* CHECK_WITHOUT_PLAN */
}
#endif /* FFTW_COMPARE */

#ifndef LOG_AUTO
static void
__report_func(double vals[][2]) {
    double t_kavg = 0;
    double t_favg = 0;

    for (int32_t i = 0; i < TEST_COUNT; i++) {
        t_kavg += vals[i][VALS_KFFT];
        t_favg += vals[i][VALS_FFTW];
    }
    #ifndef CHECK_WITHOUT_PLAN
    fprintf(stdout, "%s (%d - %s):\n\n", "Algoritm work time", TEST_COUNT, "checks");
    #else
    fprintf(stdout, "%s (%d - %s):\n\n", "Algoritm work time (without initialization)", TEST_COUNT,
            "checks");
    #endif

    #if defined FFTW_COMPARE
    fprintf(stdout, "%s\n", " --STEP-----KFFT-------FFTW----");
    for (int32_t i = 0; i < TEST_COUNT; i++) {
        fprintf(stdout, "| %5d | %8.2f | %8.2f  |\n", i, vals[i][VALS_KFFT], vals[i][VALS_FFTW]);
    }
    fprintf(stdout, "%s\n", " ------------------------------");
    fprintf(stdout, "%s:\t%8.2f\n", "Average work time for KFFT", t_kavg / TEST_COUNT);
    fprintf(stdout, "%s:\t%8.2f\n", "Average work time for FFTW", t_favg / TEST_COUNT);
    #else
    fprintf(stdout, "%s\n", " --STEP------KFFT--");
    for (int32_t i = 0; i < TEST_COUNT; i++) {
        fprintf(stdout, "| %5d | %8.2f |\n", i, vals[i][VALS_KFFT]);
    }
    fprintf(stdout, "%s\n", " ------------------ ");
    fprintf(stdout, "%s:\t%8.2f\n", "Average work time for KFFT", t_kavg / TEST_COUNT);
    #endif
}

#else
static void
__report_func(double vals[][2]) {
    double t_kavg = 0;
    double t_favg = 0;

    for (int32_t i = 0; i < TEST_COUNT; i++) {
        t_kavg += vals[i][VALS_KFFT];
    #if defined FFTW_COMPARE
        t_favg += vals[i][VALS_FFTW];
    #endif
    }
    for (int32_t i = 0; i < TEST_COUNT; i++) {
        fprintf(stdout, "%8.2f", vals[i][VALS_KFFT]);
    }
    #if defined FFTW_COMPARE
    fprintf(stdout, "\n");
    for (int32_t i = 0; i < TEST_COUNT; i++) {
        fprintf(stdout, "%8.2f", vals[i][VALS_FFTW]);
    }
    #endif
    fprintf(stdout, "\n");
    fprintf(stdout, "%8.2f\n", t_kavg / TEST_COUNT);
    #if defined FFTW_COMPARE
    fprintf(stdout, "%8.2f\n", t_favg / TEST_COUNT);
    #endif
}
#endif

int
main(int argc, char* argv[]) {
    if (argc > 1) {
        double ivals[TEST_COUNT][2];
        size_t size = atoi(argv[1]);

        kfft_scalar* amp_scalar = malloc((size + 1) * sizeof(kfft_scalar));
        kfft_cpx* kfft_spectr = malloc((size + 1) * sizeof(kfft_cpx));
#if defined FFTW_COMPARE
        fftw_complex* fftw_spectr = malloc((size + 1) * sizeof(fftw_complex));
#endif

        for (int32_t i = 0; i < TEST_COUNT; i++) {
            // Zero memory
            memset(amp_scalar, 0, size * sizeof(kfft_scalar));
            memset(kfft_spectr, 0, size * sizeof(kfft_cpx));

            srand((unsigned)time(NULL));
            for (size_t j = 0; j < size; j++) {
                amp_scalar[j] = rand();
            }

            ivals[i][VALS_KFFT] = kfft_ktest(amp_scalar, kfft_spectr, size);
#if defined FFTW_COMPARE
            ivals[i][VALS_FFTW] = __fftw_test((double*)amp_scalar, fftw_spectr, size);
#endif
        }

        free(kfft_spectr);
#if defined FFTW_COMPARE
        free(fftw_spectr);
#endif
        free(amp_scalar);

        __report_func(ivals);

        return 0;
    } else {
        fprintf(stderr, "%s\n", "Need sequense size as parameter");
        return 1;
    }
}
