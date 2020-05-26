#include "kfft.h"

#if defined(KFFT_RADER_ALGO)
    #define CHECK_PLAN_NOTPRIME(S)                                                                 \
        (((S)->flags & (KFFT_FLAG_GENERIC_ONLY | KFFT_FLAG_GENERIC))) ||                           \
            (((S)->factors[0]) <= KFFT_RADER_LIMIT)
#else
    #define CHECK_PLAN_NOTPRIME(S) true
#endif

static kfft_comp_t*
kfft_config_lvlcpx(const uint32_t nfft, const uint32_t flags, const uint8_t level, kfft_pool_t* A,
                   size_t* lenmem); /* forward declaration configure function */

#ifdef KFFT_TRACE

static void
kfft_trace_plan(kfft_comp_t* P) {
    kfft_trace_core(P->level, "%s: %p", "Create KFFT complex plan", (void*)P);
    kfft_trace("\n\t %s - %u", "nfft", P->nfft);
    kfft_trace("\n\t %s - %u", "level", P->level);
    kfft_trace("\n\t %s - %u : ", "flags", P->flags);

    if (P->flags) {
        for (uint8_t i = 0; i < sizeof(uint32_t); i++) {
            if (P->flags & (1 << i)) {
                switch (i) {
                case 0:
                    kfft_trace("| %s ", "KFFT_FLAG_INVERSE");
                    break;
                case 1:
                    kfft_trace("| %s ", "KFFT_FLAG_RENEW");
                    break;
                case 2:
                    kfft_trace("| %s ", "KFFT_FLAG_GENERIC");
                    break;
                case 3:
                    kfft_trace("| %s ", "KFFT_FLAG_GENERIC_ONLY");
                    break;
                }
            }
        }
    } else {
        kfft_trace("| %s ", "KFFT_FLAG_NORMAL");
    }

    kfft_trace("\n\t %s - %d", "factors count", P->fac_count);
    if (P->fac_count)
        kfft_trace("\n\t %u %s", P->nfft, "fac -");

    for (uint32_t i = 0; i < P->fac_count; i++)
        kfft_trace(" %u", P->factors[i]);

    kfft_trace("\n\t %s - %u", "primes count", P->prm_count);

    #if !defined(KFFT_MEMLESS_MODE)
    for (uint32_t i = 0; i < P->prm_count; i++) {
        kfft_trace("\n\t %d %s", P->primes[i].prime, "qidx -");
        for (uint32_t j = 0; j < P->primes[i].prime - 1; j++)
            kfft_trace(" %u", P->primes[i].qidx[j]);
        kfft_trace("\n\t %d %s", P->primes[i].prime, "pidx -");
        for (uint32_t j = 0; j < P->primes[i].prime - 1; j++)
            kfft_trace(" %u", P->primes[i].pidx[j]);
    }

    kfft_trace("\n\t %s - %p", "twiddles", (void*)(P->twiddles));
    #endif /* KFFT_MEMLESS_MODE */
    kfft_trace("%s\n", "");
}

#endif /* KFFT_TRACE */

#if defined(KFFT_RADER_ALGO)
static inline void
kfft_rader_idxs(uint32_t* idx, const uint32_t root, const uint32_t size) {
    for (uint32_t i = 0; i < size - 1; i++) {
        idx[i] = kfft_math_modpow(root, i, size);
    }
}

    #if defined(KFFT_MEMLESS_MODE)
        #define RAD_PRIME_IDX(i, P) kfft_math_modpow((P)->q, i, (P)->prime)
        #define RAD_INVERSE_IDX(i, P) kfft_math_modpow((P)->p, i, (P)->prime)
    #else
        #define RAD_PRIME_IDX(i, P) (P)->qidx[i]
        #define RAD_INVERSE_IDX(i, P) (P)->pidx[i]
    #endif

#endif /* KFFT_RADER_ALGO */

/*
    Inline needed source modules for transformation .
    Select SIMD files HERE //TODO
 */

#include "kfft_conv.c"    /* Complex sequenses convolution */
#include "kfft_bfly.c"    /* Butterfly transformations (Cooley - Tukey)*/
#include "kfft_generic.c" /* Generic or Rader algorithm for prime-size lengt sequences*/

#include "kfft_work.c" /* Main work recursive procedure */

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static inline void
kf_factor(kfft_comp_t* st) {
    uint32_t p = 4, n = st->nfft, *facbuf = st->factors;
#if defined(KFFT_RADER_ALGO)
    kfft_splan_t* pbuf = st->primes;
#endif /* KFFT_RADER_ALGO */

    unsigned floor_sqrt;
    floor_sqrt = (unsigned)KFFT_SQRT((double)st->nfft);

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
            case 4:
                p = 2;
                break;
            case 2:
                p = 3;
                break;
            default:
                p += 2;
                break;
            }
            if (p > floor_sqrt)
                p = n; /* no more factors, skip to end */
        }
#if defined(KFFT_RADER_ALGO)
        if (st->level < KFFT_PLAN_LEVEL && p > KFFT_BFLY_LEVEL && p > KFFT_RADER_LIMIT) {
            pbuf->prime = p;
            st->prm_count++;
            pbuf++;
        }
#endif /* KFFT_RADER_ALGO */
        st->fac_count++;

        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

static inline kfft_return_t
kfft_kinit(kfft_comp_t* st) {

    kfft_return_t ret = KFFT_RET_SUCCESS;
    /* Generate twiddles  */
#if defined(KFFT_MEMLESS_MODE)
    (void)st; // unused warning disable
#else
    #if defined(KFFT_RADER_ALGO)
    if (CHECK_PLAN_NOTPRIME(st))
    #endif /* KFFT_RADER_ALGO */
        for (uint32_t i = 0; i < st->nfft; ++i) {
            st->twiddles[i] = kfft_kernel_twiddle(i, st->nfft, st->flags & KFFT_FLAG_INVERSE);
        }
#endif     /* KFFT_MEMLESS_MODE */

#if defined(KFFT_RADER_ALGO)
    if (st->prm_count > 0) {
        if (!(st->flags & KFFT_FLAG_GENERIC_ONLY)) {

            for (uint32_t i = 0; i < st->prm_count; i++) {
                kfft_splan_t* sP = &(st->primes[i]);
                uint32_t len = sP->prime - 1;

                sP->q = kfft_math_prmn(sP->prime);
                sP->p = kfft_math_prmni(sP->q, sP->prime);

    #if !defined(KFFT_MEMLESS_MODE)
                sP->qidx = kfft_pool_alloc(st->object.mmgr, sizeof(uint32_t) * len);
                sP->pidx = kfft_pool_alloc(st->object.mmgr, sizeof(uint32_t) * len);

                if (sP->qidx && sP->pidx) {

                    kfft_rader_idxs(sP->qidx, sP->q, sP->prime);
                    kfft_rader_idxs(sP->pidx, sP->p, sP->prime);
    #endif /* not KFFT_MEMLESS_MODE */
                    sP->splan =
                        kfft_config_lvlcpx(len, KFFT_CHECK_FLAGS(st->flags & (~KFFT_FLAG_INVERSE)),
                                           st->level + 1, st->object.mmgr, NULL);
                    sP->splani =
                        kfft_config_lvlcpx(len, (KFFT_CHECK_FLAGS(st->flags | KFFT_FLAG_INVERSE)),
                                           st->level + 1, st->object.mmgr, NULL);

                    sP->shuffle_twiddles = kfft_pool_alloc(st->object.mmgr, sizeof(kfft_cpx) * len);
                    if (sP->shuffle_twiddles) {

                        for (uint32_t j = 0; j < len; j++) {
                            uint32_t ip = RAD_INVERSE_IDX(j, sP);

                            sP->shuffle_twiddles[j] =
                                kfft_kernel_twiddle(ip, sP->prime, st->flags & KFFT_FLAG_INVERSE);
                        }

                        ret = kfft_eval_cpx(sP->splan, sP->shuffle_twiddles, sP->shuffle_twiddles);
                    } else {
                        ret = KFFT_RET_ALLOC_FAIL;
                    }
    #if !defined(KFFT_MEMLESS_MODE)
                } else {
                    ret = KFFT_RET_ALLOC_FAIL;
                }
    #endif /* not KFFT_MEMLESS_MODE */
            }
        }
    }
#endif /* KFFT_RADER_ALGO */
    return ret;
}

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t flags, const uint8_t level, kfft_comp_t* st) {

    kf_factor(st);
    size_t ret = sizeof(kfft_comp_t);

#if !defined(KFFT_RADER_ALGO)
    KFFT_UNUSED_VAR(level);
    KFFT_UNUSED_VAR(st);
#endif /* KFFT_RADER_ALGO */

#if !defined(KFFT_MEMLESS_MODE)
    #if defined(KFFT_RADER_ALGO)
    if (CHECK_PLAN_NOTPRIME(st))
    #endif /* KFFT_RADER_ALGO */
        ret += sizeof(kfft_cpx) * nfft;
#else
    KFFT_UNUSED_VAR(nfft);
#endif

    if (!(flags & KFFT_FLAG_GENERIC_ONLY)) {
#if defined(KFFT_RADER_ALGO)
        /* Recursive calculate memory for plan and all subplans */
        for (size_t i = 0; i < st->prm_count; i++) {
            size_t snfft = st->primes[i].prime - 1;

            size_t delta_mem = 0;
            kfft_config_lvlcpx(snfft, flags, level + 1, NULL, &delta_mem);
            delta_mem *= 2; // splan and splani

            delta_mem += sizeof(uint32_t) * snfft * 2; // index table
            delta_mem += sizeof(kfft_cpx) * snfft;     // shuffle twiddles

            ret += delta_mem;
        }
#endif /* KFFT_RADER_ALGO */
    }
    return ret;
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kfft-specific function.
 * */

static kfft_comp_t*
kfft_config_lvlcpx(const uint32_t nfft, const uint32_t flags, const uint8_t level, kfft_pool_t* A,
                   size_t* lenmem) {
    kfft_comp_t* st = NULL;

    kfft_comp_t tmp;
    KFFT_ZEROMEM(&tmp, sizeof(kfft_comp_t));

    tmp.nfft = nfft;
    tmp.flags = flags;
    tmp.level = level;

    size_t memneeded = kfft_calculate(nfft, flags, level, &tmp);

    KFFT_ALGO_PLAN_PREPARE(st, flags, kfft_comp_t, memneeded, A, lenmem);
    if (st) {
        memcpy(&(tmp.object), &(st->object), sizeof(kfft_object_t));
        memcpy(st, &tmp, sizeof(kfft_comp_t));

#if !defined(KFFT_MEMLESS_MODE)

    #if defined(KFFT_RADER_ALGO)
        if (CHECK_PLAN_NOTPRIME(st)) {
    #endif /* KFFT_RADER_ALGO */

            st->twiddles = kfft_pool_alloc(st->object.mmgr, sizeof(kfft_cpx) * st->nfft);
            if (st->twiddles == NULL) {
                KFFT_ALGO_PLAN_TERMINATE(st, A);
                return NULL;
            }

    #if defined(KFFT_RADER_ALGO)
        }
    #endif /* KFFT_RADER_ALGO */

#endif /* not KFFT_MEMLESS_MODE */
        if (kfft_kinit(st) != KFFT_RET_SUCCESS) {
            KFFT_ALGO_PLAN_TERMINATE(st, A);
            return NULL;
        }

#ifdef KFFT_TRACE
        kfft_trace_plan(st);
#endif
    }
    return st;
}

kfft_comp_t*
kfft_config_cpx(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    return kfft_config_lvlcpx(nfft, flags, 0, A, lenmem);
}

kfft_return_t
kfft_eval_cpx(kfft_comp_t* cfg, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    if (cfg->flags & KFFT_FLAG_GENERIC_ONLY) {
        if (fin != fout)
            memcpy(fout, fin, sizeof(kfft_cpx) * cfg->nfft);
        ret = kf_work(fout, NULL, 1, 1, 0, cfg);
    } else {
        if (fin == fout) {
            // NOTE: this is not really an in-place FFT algorithm.
            // It just performs an out-of-place FFT into a temp buffer
            kfft_cpx* tmpbuf =
                (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * cfg->nfft, KFFT_PLAN_ALIGN(cfg));
            if (tmpbuf) {
                KFFT_ALLOCA_CLEAR(tmpbuf, sizeof(kfft_cpx) * cfg->nfft);

                ret = kf_work(tmpbuf, fin, 1, 1, cfg->factors, cfg);

                if (ret == KFFT_RET_SUCCESS)
                    memcpy(fout, tmpbuf, sizeof(kfft_cpx) * cfg->nfft);

                KFFT_TMP_FREE(tmpbuf, KFFT_PLAN_ALIGN(cfg));
            } else {
                ret = KFFT_RET_BUFFER_FAIL;
            }
        } else {
            ret = kf_work(fout, fin, 1, 1, cfg->factors, cfg);
        }
    }

    if (ret == KFFT_RET_SUCCESS) {
        for (uint32_t i = 0; i < cfg->nfft; i++)
            if (cfg->flags & KFFT_FLAG_INVERSE)
                C_DIVBYSCALAR(fout[i], cfg->nfft);
    }

    return ret;
}
