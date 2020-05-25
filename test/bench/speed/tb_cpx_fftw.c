#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>

#include <fftw3.h>

#ifndef TEST_COUNT
    #define TEST_COUNT 5
#endif
#if defined(TEST_HALF_SCALAR)
    #define FFTW(X) fftwf_##X
    #define fftw_scalar float
#else
    #define FFTW(X) fftw_##X
    #define fftw_scalar double
#endif

static double
fft_test(FFTW(complex) * tbuf, FFTW(complex) * ftmp, size_t size) {

    memcpy(ftmp, tbuf, size * sizeof(FFTW(complex)));
#ifdef CHECK_WITH_PLAN
    clock_t t_start = clock();

    FFTW(plan) plan = FFTW(plan_dft_1d)(size, ftmp, tbuf, FFTW_FORWARD, FFTW_ESTIMATE);
    if (plan == NULL) {
        free(tbuf);
        return -1;
    }

    FFTW(execute)(plan);
    FFTW(destroy_plan)(plan);

    clock_t t_ret = clock();
#else
    FFTW(plan) plan = FFTW(plan_dft_1d)(size, ftmp, tbuf, FFTW_FORWARD, FFTW_ESTIMATE);
    if (plan == NULL) {
        free(tbuf);
        return -1;
    }

    clock_t t_start = clock();
    FFTW(execute)(plan);
    clock_t t_ret = clock();

    FFTW(destroy_plan)(plan);
#endif
    return (t_ret - t_start) * 1000 / CLOCKS_PER_SEC;
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

        FFTW(complex)* fft_spectr = FFTW(malloc)(2 * size * sizeof(FFTW(complex)));
        if (fft_spectr) {
            FFTW(complex)* temp = fft_spectr + size;

            for (int32_t i = 0; i < TEST_COUNT; i++) {
                memset(fft_spectr, 0, size * sizeof(FFTW(complex)));

                srand(time(NULL));
                for (size_t j = 0; j < size; j++) {
                    fft_spectr[j][0] = (fftw_scalar)rand();
                    fft_spectr[j][1] = (fftw_scalar)rand();
                }

                double ret = fft_test(fft_spectr, temp, size);
                if (ret < 0)
                    return -1;
                ivals[i] = ret;
            }
            stdout_time(ivals);

            FFTW(free)(fft_spectr);
        }
        return 0;
    } else {
        fprintf(stderr, "%s\n", "Need sequense size as parameter");
        return 1;
    }
}
