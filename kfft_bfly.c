static inline void
kf_bfly2(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, uint32_t m) {
    kfft_cpx* Fout2;
    uint32_t twidx = 0;
    kfft_cpx t;
    Fout2 = Fout + m;
    do {
        C_MUL(t, *Fout2, TWIDDLE(twidx, st));
        twidx += fstride;
        C_SUB(*Fout2, *Fout, t);
        C_ADDTO(*Fout, t);
        ++Fout2;
        ++Fout;
    } while (--m);
}

static inline void
kf_bfly4(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, const uint32_t m) {
    uint32_t tw1, tw2, tw3;
    kfft_cpx scratch[6];
    uint32_t k = m;
    const uint32_t m2 = 2 * m;
    const uint32_t m3 = 3 * m;

    tw1 = tw2 = tw3 = 0;

    do {
        C_MUL(scratch[0], Fout[m], TWIDDLE(tw1, st));
        C_MUL(scratch[1], Fout[m2], TWIDDLE(tw2, st));
        C_MUL(scratch[2], Fout[m3], TWIDDLE(tw3, st));

        C_SUB(scratch[5], *Fout, scratch[1]);
        C_ADDTO(*Fout, scratch[1]);
        C_ADD(scratch[3], scratch[0], scratch[2]);
        C_SUB(scratch[4], scratch[0], scratch[2]);
        C_SUB(Fout[m2], *Fout, scratch[3]);
        tw1 += fstride;
        tw2 += fstride * 2;
        tw3 += fstride * 3;
        C_ADDTO(*Fout, scratch[3]);

        if (st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        } else {
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
        ++Fout;
    } while (--k);
}

static inline void
kf_bfly3(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, uint32_t m) {
    uint32_t k = m;
    const uint32_t m2 = 2 * m;
    uint32_t tw1, tw2;
    kfft_cpx scratch[5];
    kfft_cpx epi3;
    epi3 = TWIDDLE(fstride * m, st);

    tw1 = tw2 = 0;

    do {
        C_MUL(scratch[1], Fout[m], TWIDDLE(tw1, st));
        C_MUL(scratch[2], Fout[m2], TWIDDLE(tw2, st));

        C_ADD(scratch[3], scratch[1], scratch[2]);
        C_SUB(scratch[0], scratch[1], scratch[2]);
        tw1 += fstride;
        tw2 += fstride * 2;

        Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
        Fout[m].i = Fout->i - HALF_OF(scratch[3].i);

        C_MULBYSCALAR(scratch[0], epi3.i);

        C_ADDTO(*Fout, scratch[3]);

        Fout[m2].r = Fout[m].r + scratch[0].i;
        Fout[m2].i = Fout[m].i - scratch[0].r;

        Fout[m].r -= scratch[0].i;
        Fout[m].i += scratch[0].r;

        ++Fout;
    } while (--k);
}

static inline void
kf_bfly5(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, uint32_t m) {
    kfft_cpx *Fout0, *Fout1, *Fout2, *Fout3, *Fout4;
    uint32_t u;
    kfft_cpx scratch[13];
    kfft_cpx ya, yb;
    ya = TWIDDLE(fstride * m, st);
    yb = TWIDDLE(fstride * 2 * m, st);

    Fout0 = Fout;
    Fout1 = Fout0 + m;
    Fout2 = Fout0 + 2 * m;
    Fout3 = Fout0 + 3 * m;
    Fout4 = Fout0 + 4 * m;

    for (u = 0; u < m; ++u) {
        scratch[0] = *Fout0;

        C_MUL(scratch[1], *Fout1, TWIDDLE(u * fstride, st));
        C_MUL(scratch[2], *Fout2, TWIDDLE(2 * u * fstride, st));
        C_MUL(scratch[3], *Fout3, TWIDDLE(3 * u * fstride, st));
        C_MUL(scratch[4], *Fout4, TWIDDLE(4 * u * fstride, st));

        C_ADD(scratch[7], scratch[1], scratch[4]);
        C_SUB(scratch[10], scratch[1], scratch[4]);
        C_ADD(scratch[8], scratch[2], scratch[3]);
        C_SUB(scratch[9], scratch[2], scratch[3]);

        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r, ya.r) + S_MUL(scratch[8].r, yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i, ya.r) + S_MUL(scratch[8].i, yb.r);

        scratch[6].r = S_MUL(scratch[10].i, ya.i) + S_MUL(scratch[9].i, yb.i);
        scratch[6].i = -S_MUL(scratch[10].r, ya.i) - S_MUL(scratch[9].r, yb.i);

        C_SUB(*Fout1, scratch[5], scratch[6]);
        C_ADD(*Fout4, scratch[5], scratch[6]);

        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r, yb.r) + S_MUL(scratch[8].r, ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i, yb.r) + S_MUL(scratch[8].i, ya.r);
        scratch[12].r = -S_MUL(scratch[10].i, yb.i) + S_MUL(scratch[9].i, ya.i);
        scratch[12].i = S_MUL(scratch[10].r, yb.i) - S_MUL(scratch[9].r, ya.i);

        C_ADD(*Fout2, scratch[11], scratch[12]);
        C_SUB(*Fout3, scratch[11], scratch[12]);

        ++Fout0;
        ++Fout1;
        ++Fout2;
        ++Fout3;
        ++Fout4;
    }
}
