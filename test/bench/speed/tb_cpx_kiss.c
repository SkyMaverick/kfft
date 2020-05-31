#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "kiss_fft.h"

#ifndef TEST_COUNT
    #define TEST_COUNT 5
#endif

static double
fft_test(kiss_fft_cpx* tbuf, kiss_fft_cpx* ftmp, size_t size) {
    memcpy(ftmp, tbuf, size * sizeof(kiss_fft_cpx));
#ifdef CHECK_WITH_PLAN
    clock_t t_start = clock();

    kiss_fft_cfg plan = kiss_fft_alloc(size, 0, NULL, NULL);
    if (plan == NULL) {
        free(tbuf);
        return -1;
    }

    kiss_fft(plan, ftmp, tbuf);
    free(plan);

    clock_t t_ret = clock();
#else
    kiss_fft_cfg plan = kiss_fft_alloc(size, 0, NULL, NULL);
    if (plan == NULL) {
        free(tbuf);
        return -1;
    }

    clock_t t_start = clock();
    kiss_fft(plan, ftmp, tbuf);
    clock_t t_ret = clock();

    free(plan);
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
        size_t size = atol(argv[1]);

        kiss_fft_cpx* fft_spectr = calloc(size * 2, sizeof(kiss_fft_cpx));
        if (fft_spectr) {
            kiss_fft_cpx* temp = fft_spectr + size;

            for (int32_t i = 0; i < TEST_COUNT; i++) {
                memset(fft_spectr, 0, size * sizeof(kiss_fft_cpx));

                srand(time(NULL));
                for (size_t j = 0; j < size; j++) {
                    fft_spectr[j].r = rand();
                    fft_spectr[j].i = rand();
                }

                double ret = fft_test(fft_spectr, temp, size);
                if (ret < 0)
                    return -1;
                ivals[i] = ret;
            }
            stdout_time(ivals);

            free(fft_spectr);
        }
        return 0;
    } else {
        fprintf(stderr, "%s\n", "Need sequense size as parameter");
        return 1;
    }
}
