#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>

#include "kfft.h"

#ifndef TEST_COUNT
    #define TEST_COUNT 12
#endif

static double
kfft_ktest(kfft_cpx* tbuf, size_t size) {

    size_t memneed = size * sizeof(kfft_cpx);

    kfft_cpx* ftmp = kfft_malloc(size * sizeof(kfft_cpx));
    if (ftmp) {
        memcpy(ftmp, tbuf, memneed);
#ifdef CHECK_WITH_PLAN
        clock_t t_start = clock();

        kfft_comp_t* plan = kfft_config_cpx(size, KFFT_FLAG_NORMAL, NULL, NULL);
        if (plan == NULL) {
            free(ftmp);
            return -1;
        }

        kfft_eval_cpx(plan, ftmp, tbuf);
        kfft_cleanup(plan);

        clock_t t_ret = clock();
#else
        kfft_comp_t* plan = kfft_config_cpx(size, KFFT_FLAG_NORMAL, NULL, NULL);
        if (plan == NULL) {
            kfft_free(ftmp);
            return -1;
        }

        clock_t t_start = clock();
        kfft_eval_cpx(plan, ftmp, tbuf);
        clock_t t_ret = clock();

        kfft_cleanup(plan);
#endif
        kfft_free(&ftmp);
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

        kfft_cpx* kfft_spectr = kfft_malloc(size * sizeof(kfft_cpx));
        if (kfft_spectr) {

            for (int32_t i = 0; i < TEST_COUNT; i++) {
                memset(kfft_spectr, 0, size * sizeof(kfft_cpx));

                srand(time(NULL));
                for (size_t j = 0; j < size; j++) {
                    kfft_spectr[j].r = rand();
                    kfft_spectr[j].i = rand();
                }

                double ret = kfft_ktest(kfft_spectr, size);
                if (ret < 0)
                    return -1;
                ivals[i] = ret;
            }
            stdout_time(ivals);

            kfft_free(&kfft_spectr);
        }
        return 0;
    } else {
        fprintf(stderr, "%s\n", "Need sequense size as parameter");
        return 1;
    }
}
