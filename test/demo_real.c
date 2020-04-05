#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kfft.h"

static inline uint32_t
get_size(const uint32_t n) {
    size_t memneeded = 0;
    kfft_config_scalar(n, 0, 0, &memneeded);
    return memneeded;
}

static void
display_info(void) {
    kfft_info_t info;

    kfft_info(&info);

    fprintf(stdout, "LibKFFT version : %d.%d.%d\n\n", info.vmajor, info.vminor, info.vpatch);

    fprintf(stdout, "%s - %s\n", "Enable trace messages",
            (info.flags & KFFT_INFO_TRACE) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use SIMD instructions",
            (info.flags & KFFT_INFO_USE_SIMD) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use alloca() function",
            (info.flags & KFFT_INFO_USE_ALLOCA) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use system math functions",
            (info.flags & KFFT_INFO_USE_SYSMATH) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use Rader algoritm",
            (info.flags & KFFT_INFO_RADER_ALGO) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable lesser memory mode",
            (info.flags & KFFT_INFO_MEMLESS_MODE) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable half scalar mode",
            (info.flags & KFFT_INFO_HALF_SCALAR) ? "YES" : "NO");
}

int
main(int argc, char* argv[]) {
    if (argc > 1) {
        display_info();
        printf("\n");

        kfft_scalar* amp_scalar = calloc(argc, sizeof(double));
        for (int32_t i = 1; i < argc; i++) {
            amp_scalar[i - 1] = atof(argv[i]);
            printf("%5.3f ", amp_scalar[i - 1]);
        }
        printf("\n");

        kfft_cpx* FOut = calloc(argc, sizeof(kfft_cpx));

        size_t memneed = get_size(argc - 1);

        printf("Create forward config for %d len\n", argc - 1);
        kfft_sclr_t* FCfg = kfft_config_scalar(argc - 1, KFFT_FLAG_NORMAL, 0, NULL);

        printf("Forward FFT transform\n");
        kfft_eval_scalar(FCfg, amp_scalar, FOut);

        printf("\nRESULT >>> ");
        for (int32_t i = 0; i < argc - 1; i++) {
            printf("r%5.3fi%5.3f | ", FOut[i].r, FOut[i].i);
        }
        printf("\n\n");

        printf("Create inverse config for %d len\n", argc - 1);

        kfft_config_scalar(argc - 1, KFFT_FLAG_INVERSE | KFFT_FLAG_RENEW, KFFT_PLAN_ALLOCATOR(FCfg),
                           &memneed);
        memset(amp_scalar, 0, argc * sizeof(double));

        printf("Inverse FFT transform\n");
        kfft_evali_scalar(FCfg, FOut, amp_scalar);

        printf("\nRESULT >>> ");
        for (int32_t i = 0; i < argc - 1; i++) {
            printf("%5.3f | ", amp_scalar[i]);
        }

        printf("\n\n");

        free(FOut);
        free(amp_scalar);
        kfft_free(FCfg);

        return 0;
    } else {
        fprintf(stderr, "%s\n", "Need scalar array as arguments");
        return 1;
    }
}
