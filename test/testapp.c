#include <stdlib.h>
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

        kiss_fft_cpx* FOut = calloc (argc, sizeof(kiss_fft_cpx));

        kiss_fft_cfg  FCfg = kiss_fft_alloc (argc-1, 0, 0, NULL, NULL);

        kiss_fft (FCfg, amp_scalar, FOut);
        
        for (int i = 0; i < argc-1; i++) {
            printf("r%5.3fi%5.3f | ", FOut[i].r, FOut[i].i);    
        }
        printf ("\n");
        

        free (FOut);
        free (amp_scalar);
        kiss_fft_free (&FCfg);
        
        return 0;
    } else {
        fprintf (stderr, "%s\n","Need scalar array as arguments");
        return 1;
    }
}
