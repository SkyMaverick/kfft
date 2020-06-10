#include "kfa_defs.h"

static unsigned
work_complex_2d(kfft_cpx* in, state_t* st) {
    KCB_TYPE(eval2_cpx) eval = KFFT_CALLBACK(st, eval2_cpx);
    KCB_TYPE(config2_cpx) config = KFFT_CALLBACK(st, config2_cpx);
    KCB_TYPE(shift2_cpx) shift = KFFT_CALLBACK(st, shift2_cpx);

    if ((eval == NULL) || (config == NULL) || (shift == NULL))
        return KFA_RET_FAIL_UNREAL;

    kfft_comp2_t* plan = config(st->dims.x, st->dims.y, sm2kfl(st), NULL, NULL);
    if (plan) {
        unsigned ret;
        if (KFA_CHECK(st, SHIFT)) {
            if (KFA_CHECK(st, INVERSE)) {
                shift(in, NULL, st->dims.x, st->dims.y, true, KFFT_PLAN_MMGR(plan));
                ret = eval(plan, in, in);
            } else {
                ret = eval(plan, in, in);
                if (ret == KFFT_RET_SUCCESS)
                    shift(in, NULL, st->dims.x, st->dims.y, false, KFFT_PLAN_MMGR(plan));
            }
        } else {
            ret = eval(plan, in, in);
        }
        KRNL_FUNCS(st).cb_cleanup(plan);

        if (ret != KFFT_RET_SUCCESS) {
            KRNL_FUNCS(st).cb_strerr(ret);
        } else {
            return KFA_RET_SUCCESS;
        }
    }
    return KFA_RET_FAIL_INTRNL;
}

static unsigned
work_complex_sparse(kfft_cpx* in, state_t* st) {
    KCB_TYPE(eval_sparse_cpx) eval = KFFT_CALLBACK(st, eval_sparse_cpx);
    KCB_TYPE(config_sparse_cpx) config = KFFT_CALLBACK(st, config_sparse_cpx);
    KCB_TYPE(shift_sparse_cpx) shift = KFFT_CALLBACK(st, shift_sparse_cpx);

    if ((eval == NULL) || (config == NULL) || (shift == NULL))
        return KFA_RET_FAIL_UNREAL;

    size_t cpx_len = st->lenght / 2;

    kfft_csparse_t* plan = config(cpx_len, sm2kfl(st), st->sparse.dx, st->sparse.sx, NULL, NULL);
    if (plan) {
        unsigned ret;
        if (KFA_CHECK(st, SHIFT)) {
            if (KFA_CHECK(st, INVERSE)) {
                shift(in, NULL, cpx_len, st->sparse.dx, st->sparse.sx, true, KFFT_PLAN_MMGR(plan));
                ret = eval(plan, in, in);
            } else {
                ret = eval(plan, in, in);
                if (ret == KFFT_RET_SUCCESS)
                    shift(in, NULL, cpx_len, st->sparse.dx, st->sparse.sx, false,
                          KFFT_PLAN_MMGR(plan));
            }
        } else {
            ret = eval(plan, in, in);
        }
        KRNL_FUNCS(st).cb_cleanup(plan);

        if (ret != KFFT_RET_SUCCESS) {
            KRNL_FUNCS(st).cb_strerr(ret);
        } else {
            return KFA_RET_SUCCESS;
        }
    }
    return KFA_RET_FAIL_INTRNL;
}

static unsigned
work_complex_normal(kfft_cpx* in, state_t* st) {
    KCB_TYPE(eval_cpx) eval = KFFT_CALLBACK(st, eval_cpx);
    KCB_TYPE(config_cpx) config = KFFT_CALLBACK(st, config_cpx);
    KCB_TYPE(shift_cpx) shift = KFFT_CALLBACK(st, shift_cpx);

    if ((eval == NULL) || (config == NULL) || (shift == NULL))
        return KFA_RET_FAIL_UNREAL;

    size_t cpx_len = st->lenght / 2;

    kfft_comp_t* plan = config(cpx_len, sm2kfl(st), NULL, NULL);
    if (plan) {
        unsigned ret;
        if (KFA_CHECK(st, SHIFT)) {
            if (KFA_CHECK(st, INVERSE)) {
                shift(in, cpx_len, true, KFFT_PLAN_MMGR(plan));
                ret = eval(plan, in, in);
            } else {
                ret = eval(plan, in, in);
                if (ret == KFFT_RET_SUCCESS)
                    shift(in, cpx_len, false, KFFT_PLAN_MMGR(plan));
            }
        } else {
            ret = eval(plan, in, in);
        }
        KRNL_FUNCS(st).cb_cleanup(plan);

        if (ret != KFFT_RET_SUCCESS) {
            KRNL_FUNCS(st).cb_strerr(ret);
        } else {
            return KFA_RET_SUCCESS;
        }
    }
    return KFA_RET_FAIL_INTRNL;
}

unsigned
work_complex_plan(kfft_scalar* in, state_t* st) {
    unsigned ret = KFA_RET_FAIL_INTRNL;

    if (KFA_CHECK(st, SPARSE))
        ret = work_complex_sparse((kfft_cpx*)in, st);
    if (KFA_CHECK(st, 2D))
        ret = work_complex_2d((kfft_cpx*)in, st);

    ret = work_complex_normal((kfft_cpx*)in, st);

    if (ret == KFA_RET_SUCCESS)
        write_stdout(in, st);

    return ret;
}
