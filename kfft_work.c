static kfft_return_t
kf_work(kfft_cpx* Fout, const kfft_cpx* f, const uint32_t fstride, uint32_t in_stride,
        uint32_t* factors, const kfft_comp_t* st) {

    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* Fout_beg = Fout;
    if (st->flags & KFFT_FLAG_GENERIC_ONLY) {
        ret = kf_bfly_generic(Fout, 1, st, 1, st->nfft);
    } else {
        const uint32_t p = *factors++; /* the radix  */
        const uint32_t m = *factors++; /* stage's fft length/p */
        const kfft_cpx* Fout_end = Fout + p * m;

        kfft_trace_core(st->level, "Work: p - %u | m - %u\n", p, m);

        if (m == 1) {
            do {
                *Fout = *f;
                f += fstride * in_stride;
            } while (++Fout != Fout_end);
        } else {
            do {
                // recursive call:
                // DFT of size m*p performed by doing
                // p instances of smaller DFTs of size m,
                // each one takes a decimated version of the input
                ret = kf_work(Fout, f, fstride * p, in_stride, factors, st);

                if (ret != KFFT_RET_SUCCESS)
                    goto bailout;

                f += fstride * in_stride;
            } while ((Fout += m) != Fout_end);
        }

        Fout = Fout_beg;

        // recombine the p smaller DFTs
        switch (p) {
        case 2:
            VEXFUNC(st, kf_bfly2, Fout, fstride, st, m);
            break;
        case 3:
            VEXFUNC(st, kf_bfly3, Fout, fstride, st, m);
            break;
        case 4:
            VEXFUNC(st, kf_bfly4, Fout, fstride, st, m);
            break;
        case 5:
            VEXFUNC(st, kf_bfly5, Fout, fstride, st, m);
            break;
        default:
            ret = kf_bfly_generic(Fout, fstride, st, m, p);
            break;
        }
    }
bailout:
    return ret;
}
