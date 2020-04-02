#include "kfft.h"
#include "kfft_math.h"

KFFT_API void
kfft_shift(kfft_cpx* buf, const uint32_t size, const bool is_inverse) {
    ssize_t k = 0;
    uint32_t c = (uint32_t)floor((kfft_scalar)size / 2);
    kfft_cpx tmp = {0, 0};
    // For odd and for even numbers of element use different algorithm
    if (size % 2 == 0) {
        for (k = 0; k < c; k++)
            C_SWAP(tmp, buf[k], buf[k + c]);
    } else {
        if (is_inverse) {
            C_CPY(tmp, buf[size - 1]);
            for (k = c - 1; k >= 0; k--) {
                C_CPY(buf[c + k + 1], buf[k]);
                C_CPY(buf[k], buf[c + k]);
            }
            C_CPY(buf[c], tmp);
        } else {
            C_CPY(tmp, buf[0]);
            for (k = 0; k < c; k++) {
                C_CPY(buf[k], buf[c + k + 1]);
                C_CPY(buf[c + k + 1], buf[k + 1]);
            }
            C_CPY(buf[c], tmp);
        }
    }
}
