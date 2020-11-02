#include "kfft.h"

#if defined(KFFT_RADER_ALGO)
    #define CHECK_PLAN_NOTPRIME(S)                                                                 \
        (((S)->flags & (KFFT_FLAG_GENERIC_ONLY | KFFT_FLAG_GENERIC))) ||                           \
            (((S)->factors[0]) <= KFFT_RADER_LIMIT)
#else
    #define CHECK_PLAN_NOTPRIME(S) true
#endif

static kfft_plan_cpx*
kfft_config_lvlcpx(const uint32_t nfft, const uint32_t flags, const uint8_t level, kfft_pool_t* A,
                   size_t* lenmem); /* forward declaration configure function */

#ifdef KFFT_TRACE

static void
kfft_trace_plan(kfft_plan_cpx* P) {
    KFFT_OMP(omp critical(trace_log)) {
        kfft_trace_core(P->level, "%s: %p", "Create KFFT complex plan", (void*)P);
        kfft_trace_raw("\n\t %s - %u", "nfft", P->nfft);
        kfft_trace_raw("\n\t %s - %u", "level", P->level);
        kfft_trace_raw("\n\t %s - %u : ", "flags", P->flags);

        if (P->flags) {
            for (uint8_t i = 0; i < sizeof(uint32_t); i++) {
                if (P->flags & (1 << i)) {
                    switch (i) {
                    case 0:
                        kfft_trace_raw("| %s ", "KFFT_FLAG_INVERSE");
                        break;
                    case 1:
                        kfft_trace_raw("| %s ", "KFFT_FLAG_RENEW");
                        break;
                    case 2:
                        kfft_trace_raw("| %s ", "KFFT_FLAG_GENERIC");
                        break;
                    case 3:
                        kfft_trace_raw("| %s ", "KFFT_FLAG_GENERIC_ONLY");
                        break;
                    }
                }
            }
        } else {
            kfft_trace_raw("| %s ", "KFFT_FLAG_NORMAL");
        }

        kfft_trace_raw("\n\t %s - %d", "factors count", P->fac_count);
        if (P->fac_count)
            kfft_trace_raw("\n\t %u %s", P->nfft, "fac -");

        for (uint32_t i = 0; i < P->fac_count; i++)
            kfft_trace_raw(" %u", P->factors[i]);

        kfft_trace_raw("\n\t %s - %u", "primes count", P->prm_count);

    #if !defined(KFFT_MEMLESS_MODE)
        for (uint32_t i = 0; i < P->prm_count; i++) {
            kfft_trace_raw("\n\t %d %s", P->primes[i].prime, "qidx -");
            for (uint32_t j = 0; j < P->primes[i].prime - 1; j++)
                kfft_trace_raw(" %u", P->primes[i].qidx[j]);
            kfft_trace_raw("\n\t %d %s", P->primes[i].prime, "pidx -");
            for (uint32_t j = 0; j < P->primes[i].prime - 1; j++)
                kfft_trace_raw(" %u", P->primes[i].pidx[j]);
        }

        kfft_trace_raw("\n\t %s - %p", "twiddles", (void*)(P->twiddles));
    #endif /* KFFT_MEMLESS_MODE */
        kfft_trace_raw("%s\n", "");
    }
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

#include "kfft_bfly.c"    /* Butterfly transformations (Cooley - Tukey)*/
#include "kfft_generic.c" /* Generic or Rader algorithm for prime-size lengt sequences*/

#include "kfft_work.c" /* Main work recursive procedure */

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static inline void
kf_factor(kfft_plan_cpx* P) {
    uint32_t p = 4, n = P->nfft, *facbuf = P->factors;
#if defined(KFFT_RADER_ALGO)
    kfft_plan_rader* pbuf = P->primes;
#endif /* KFFT_RADER_ALGO */

    unsigned floor_sqrt;
    floor_sqrt = (unsigned)KFFT_SQRT((double)P->nfft);

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
        if (P->level < KFFT_PLAN_LEVEL && p > KFFT_BFLY_LEVEL && p >= KFFT_RADER_LIMIT) {
            pbuf->prime = p;
            P->prm_count++;
            pbuf++;
        }
#endif /* KFFT_RADER_ALGO */
        P->fac_count++;

        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

static inline kfft_return_t
kfft_kinit(kfft_plan_cpx* P) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
/* Generate twiddles  */
#if defined(KFFT_MEMLESS_MODE)
    (void)P; // unused warning disable
#else
    //     #if defined(KFFT_RADER_ALGO)
    //     if (CHECK_PLAN_NOTPRIME(P))
    //     #endif /* KFFT_RADER_ALGO */
    for (uint32_t i = 0; i < P->nfft; ++i) {
        P->twiddles[i] = kfft_kernel_twiddle(i, P->nfft, P->flags & KFFT_FLAG_INVERSE);
    }
#endif /* KFFT_MEMLESS_MODE */

#if defined(KFFT_RADER_ALGO)
    if (P->prm_count > 0) {
        if (__likely__(!(P->flags & KFFT_FLAG_GENERIC_ONLY))) {
            for (uint32_t i = 0; i < P->prm_count; i++) {
                kfft_plan_rader* sP = &(P->primes[i]);
                uint32_t len = sP->prime - 1;

                sP->q = kfft_math_prmn(sP->prime);
                sP->p = kfft_math_prmni(sP->prime, sP->q);

    #if !defined(KFFT_MEMLESS_MODE)
                sP->qidx = kfft_pool_alloc(P->object.mmgr, sizeof(uint32_t) * len);
                sP->pidx = kfft_pool_alloc(P->object.mmgr, sizeof(uint32_t) * len);

                if (__likely__(sP->qidx && sP->pidx)) {
                    kfft_rader_idxs(sP->qidx, sP->q, sP->prime);
                    kfft_rader_idxs(sP->pidx, sP->p, sP->prime);
    #endif /* not KFFT_MEMLESS_MODE */
                    sP->plan =
                        kfft_config_lvlcpx(len, KFFT_CHECK_FLAGS(P->flags & (~KFFT_FLAG_INVERSE)),
                                           P->level + 1, P->object.mmgr, NULL);
                    sP->plan_inv =
                        kfft_config_lvlcpx(len, (KFFT_CHECK_FLAGS(P->flags | KFFT_FLAG_INVERSE)),
                                           P->level + 1, P->object.mmgr, NULL);

                    sP->shuffle_twiddles = kfft_pool_alloc(P->object.mmgr, sizeof(kfft_cpx) * len);
                    if (__likely__(sP->shuffle_twiddles)) {
                        for (uint32_t j = 0; j < len; j++) {
                            uint32_t ip = RAD_INVERSE_IDX(j, sP);

                            sP->shuffle_twiddles[j] =
                                kfft_kernel_twiddle(ip, sP->prime, P->flags & KFFT_FLAG_INVERSE);
                        }

                        ret = kfft_eval_cpx(sP->plan, sP->shuffle_twiddles, sP->shuffle_twiddles);
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
kfft_calculate(const uint32_t nfft, const uint32_t flags, const uint8_t level, kfft_plan_cpx* P) {
    size_t ret = sizeof(kfft_plan_cpx);

#if !defined(KFFT_RADER_ALGO)
    KFFT_UNUSED_VAR(level);
    KFFT_UNUSED_VAR(P);
#endif /* KFFT_RADER_ALGO */

#if !defined(KFFT_MEMLESS_MODE)
    //     #if defined(KFFT_RADER_ALGO)
    //     if (CHECK_PLAN_NOTPRIME(P))
    //     #endif /* KFFT_RADER_ALGO */
    ret += sizeof(kfft_cpx) * nfft;
#else
    KFFT_UNUSED_VAR(nfft);
#endif /* not KFFT_MEMLESS_MODE */

    if (__likely__(!(flags & KFFT_FLAG_GENERIC_ONLY))) {
#if defined(KFFT_RADER_ALGO)
        /* Recursive calculate memory for plan and all subplans */
        for (size_t i = 0; i < P->prm_count; i++) {
            size_t snfft = P->primes[i].prime - 1;

            size_t delta_mem = 0;
            kfft_config_lvlcpx(snfft, flags, level + 1, NULL, &delta_mem);
            delta_mem *= 2; // plan and plan_inv

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

static kfft_plan_cpx*
kfft_config_lvlcpx(const uint32_t nfft, const uint32_t flags, const uint8_t level, kfft_pool_t* A,
                   size_t* lenmem) {
    kfft_plan_cpx* P = NULL;

    kfft_plan_cpx tmp;
    KFFT_ZEROMEM(&tmp, sizeof(kfft_plan_cpx));

    tmp.nfft = nfft;
    tmp.flags = flags;
    tmp.level = level;

    kf_factor(&tmp);
    size_t memneeded = kfft_calculate(nfft, flags, level, &tmp);

    KFFT_ALGO_PLAN_PREPARE(P, flags, kfft_plan_cpx, memneeded, A, lenmem);
    if (__likely__(P)) {
        memcpy(&(tmp.object), &(P->object), sizeof(kfft_object_t));
        memcpy(P, &tmp, sizeof(kfft_plan_cpx));

#if !defined(KFFT_MEMLESS_MODE)
        //
        //     #if defined(KFFT_RADER_ALGO)
        //         if (CHECK_PLAN_NOTPRIME(P)) {
        //     #endif /* KFFT_RADER_ALGO */
        //
        P->twiddles = kfft_pool_alloc(P->object.mmgr, sizeof(kfft_cpx) * P->nfft);
        if (__unlikely__(P->twiddles == NULL)) {
            KFFT_ALGO_PLAN_TERMINATE(P, A);
            return NULL;
        }

//    #if defined(KFFT_RADER_ALGO)
//        }
//    #endif /* KFFT_RADER_ALGO */
//
#endif /* not KFFT_MEMLESS_MODE */
        if (__unlikely__(kfft_kinit(P) != KFFT_RET_SUCCESS)) {
            KFFT_ALGO_PLAN_TERMINATE(P, A);
            return NULL;
        }

#ifdef KFFT_TRACE
        kfft_trace_plan(P);
#endif
    }
    return P;
}

kfft_plan_cpx*
kfft_config_cpx(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    return kfft_config_lvlcpx(nfft, flags, 0, A, lenmem);
}

kfft_return_t
kfft_eval_cpx(kfft_plan_cpx* plan, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    if (__unlikely__(plan->flags & KFFT_FLAG_GENERIC_ONLY)) {
        if (__likely__(fin != fout))
            memcpy(fout, fin, sizeof(kfft_cpx) * plan->nfft);
        ret = kf_work(fout, NULL, 1, 1, 0, plan);
    } else {
        if (__unlikely__(fin == fout)) {
            // NOTE: this is not really an in-place FFT algorithm.
            // It just performs an out-of-place FFT into a temp buffer
            kfft_cpx* tmpbuf =
                (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan));
            if (__likely__(tmpbuf)) {
                KFFT_ALLOCA_CLEAR(tmpbuf, sizeof(kfft_cpx) * plan->nfft);

                ret = kf_work(tmpbuf, fin, 1, 1, plan->factors, plan);

                if (__likely__(ret == KFFT_RET_SUCCESS))
                    memcpy(fout, tmpbuf, sizeof(kfft_cpx) * plan->nfft);

                KFFT_TMP_FREE(tmpbuf, KFFT_PLAN_ALIGN(plan));
            } else {
                ret = KFFT_RET_BUFFER_FAIL;
            }
        } else {
            ret = kf_work(fout, fin, 1, 1, plan->factors, plan);
        }
    }

    if (__likely__(ret == KFFT_RET_SUCCESS)) {
        if (plan->flags & KFFT_FLAG_INVERSE) {
            if (__likely__(plan->level == 0)) {
                if (!(plan->flags & KFFT_FLAG_DISABLE_NORM)) {
                    for (uint32_t i = 0; i < plan->nfft; i++)
                        C_DIVBYSCALAR(fout[i], plan->nfft);
                }
            } else {
                for (uint32_t i = 0; i < plan->nfft; i++)
                    C_DIVBYSCALAR(fout[i], plan->nfft);
            }
        }
    }

    return ret;
}
