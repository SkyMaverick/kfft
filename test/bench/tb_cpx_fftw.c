#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>

#include <fftw3.h>

#ifndef TEST_COUNT
    #define TEST_COUNT 12
#endif

static double
fft_test(fftw_complex* tbuf, size_t size) {

    size_t memneed = size * sizeof(fftw_complex);

    fftw_complex* ftmp = fftw_malloc(size * sizeof(fftw_complex));
    if (ftmp) {
        memcpy(ftmp, tbuf, memneed);
#ifdef CHECK_WITH_PLAN
        clock_t t_start = clock();

        fftw_plan plan = fftw_plan_dft_1d(size, ftmp, tbuf, FFTW_FORWARD, FFTW_ESTIMATE);
        if (plan == NULL) {
            free(ftmp);
            return -1;
        }

        fftw_execute(plan);
        fftw_destroy_plan(plan);

        clock_t t_ret = clock();
#else
        fftw_plan plan = fftw_plan_dft_1d(size, ftmp, tbuf, FFTW_FORWARD, FFTW_ESTIMATE);
        if (plan == NULL) {
            free(ftmp);
            return -1;
        }

        clock_t t_start = clock();
        fftw_execute(plan);
        clock_t t_ret = clock();

        fftw_destroy_plan(plan);
#endif
        fftw_free(ftmp);
        return (t_ret - t_start) * 1000 / CLOCKS_PER_SEC;
    } else {
        return -1;
    }
}

void
stdout_time(double* T) {
    double eval = 0.0;
    // Unanalize first and last elements
    for (size_t i = 1; i < TEST_COUNT - 2; i++) {
        eval += T[i];
    }
    fprintf(stdout, "%7.5f\n", eval / (TEST_COUNT - 2));
}

int
main(int argc, char* argv[]) {
    if (argc > 1) {

        double ivals[TEST_COUNT];
        size_t size = atoi(argv[1]);

        fftw_complex* kfft_spectr = fftw_malloc(size * sizeof(fftw_complex));
        if (kfft_spectr) {

            for (int32_t i = 0; i < TEST_COUNT; i++) {
                memset(kfft_spectr, 0, size * sizeof(fftw_complex));

                srand(time(NULL));
                for (size_t j = 0; j < size; j++) {
                    double rand_arg = rand();
                    memcpy(&(kfft_spectr[j]), &rand_arg, sizeof(fftw_complex));
                }

                double ret = fft_test(kfft_spectr, size);
                if (ret < 0)
                    return -1;
                ivals[i] = ret;
            }
            stdout_time(ivals);

            fftw_free(kfft_spectr);
        }
        return 0;
    } else {
        fprintf(stderr, "%s\n", "Need sequense size as parameter");
        return 1;
    }
}
