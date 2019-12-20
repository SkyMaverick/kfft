#include <stdlib.h>
#include <string.h>
#include "kfft.h"

int
main (int argc, char* argv[])
{
    if (argc > 1) {
        double * amp_scalar = calloc (argc, sizeof(double));
        for (int i=1; i < argc; i++) {
            amp_scalar [i-1] = atof(argv[i]);
            printf ("%5.3f ", amp_scalar[i-1]);
        }
        printf ("\n");

        kfft_cpx* FOut = calloc (argc, sizeof(kfft_cpx));

        size_t memneed = kfft_get_size(argc-1);

        printf ("Create forward config for %d len\n", argc-1);
        kfft_t FCfg = kfft_config (argc-1, 0, 0, NULL);
        
        printf ("Forward FFT transform\n");
        kfft (FCfg, amp_scalar, FOut);

        for (int i = 0; i < argc-1; i++) {
            printf("r%5.3fi%5.3f | ", FOut[i].r, FOut[i].i);
        }
        printf ("\n\n\n");
        
        printf ("Create inverse config for %d len\n", argc-1);
        kfft_config (argc-1, 1, FCfg, &memneed);
        memset (amp_scalar, 0, argc * sizeof(double));

        printf ("Inverse FFT transform\n");
        kffti (FCfg, FOut, amp_scalar);

        for (int i = 0; i < argc-1; i++) {
            printf("%5.3f | ", amp_scalar[i]);
        }

        printf ("\n");

        free (FOut);
        free (amp_scalar);
        kfft_free (&FCfg);

        return 0;
    } else {
        fprintf (stderr, "%s\n","Need scalar array as arguments");
        return 1;
    }
}
