static inline unsigned
scalar_sparse_forward(kfft_scalar* in, kfft_ssparse_t* plan, state_t* st) {
    unsigned ret = KFA_RET_FAIL_INTRNL;

    KCB_TYPE(eval_sparse_scalar) eval = KFFT_CALLBACK(st, eval_sparse_scalar);
    KCB_TYPE(shift_sparse_cpx) shift = KFFT_CALLBACK(st, shift_sparse_cpx);

    kfft_cpx* out = KRNL_FUNCS(st).cb_malloc(st->out_count * sizeof(kfft_cpx));
    if (out) {

        for (uint32_t i = 0; i < st->in_count; i++)
            out[i].r = in[i];

        unsigned ret = eval(plan, in, out);

        if (ret == KFFT_RET_SUCCESS) {
            if (KFA_CHECK(st, SHIFT))
                shift(out, NULL, st->out_count, st->sparse.dx, st->sparse.sx, false,
                      KFFT_PLAN_MMGR(plan));

            write_stdout((kfft_scalar*)out, st);
            ret = KFA_RET_SUCCESS;
        } else {
            KRNL_FUNCS(st).cb_strerr(ret);
            ret = KFA_RET_FAIL_INTRNL;
        }
        KRNL_FUNCS(st).cb_free_null((void**)&out);
    }
    return ret;
}
static inline unsigned
scalar_sparse_inverse(kfft_scalar* in, kfft_ssparse_t* plan, state_t* st) {
    unsigned ret = KFA_RET_FAIL_INTRNL;

    KCB_TYPE(evali_sparse_scalar) evali = KFFT_CALLBACK(st, evali_sparse_scalar);
    KCB_TYPE(shift_sparse_cpx) shift = KFFT_CALLBACK(st, shift_sparse_cpx);

    kfft_scalar* out = KRNL_FUNCS(st).cb_malloc(st->out_count * sizeof(kfft_scalar));
    if (out) {
        for (uint32_t i = 0; i < st->in_count; i++)
            out[i] = ((kfft_cpx*)in)[i].r;

        if (KFA_CHECK(st, SHIFT))
            shift((kfft_cpx*)in, NULL, st->out_count, st->sparse.dx, st->sparse.sx, true,
                  KFFT_PLAN_MMGR(plan));

        unsigned ret = evali(plan, (kfft_cpx*)in, out);

        if (ret == KFFT_RET_SUCCESS) {
            write_stdout(out, st);
            ret = KFA_RET_SUCCESS;
        } else {
            KRNL_FUNCS(st).cb_strerr(ret);
            ret = KFA_RET_FAIL_INTRNL;
        }
        KRNL_FUNCS(st).cb_free_null((void**)&out);
    }
    return ret;
}
static unsigned
work_scalar_sparse(kfft_scalar* in, state_t* st) {
    KCB_TYPE(config_sparse_scalar) config = KFFT_CALLBACK(st, config_sparse_scalar);
    if (config == NULL)
        return KFA_RET_FAIL_UNREAL;

    kfft_ssparse_t* plan =
        config(st->in_count, sm2kfl(st), st->sparse.dx, st->sparse.sx, NULL, NULL);
    if (plan) {
        unsigned ret = (KFA_CHECK(st, INVERSE)) ? scalar_sparse_inverse(in, plan, st)
                                                : scalar_sparse_forward(in, plan, st);
        KRNL_FUNCS(st).cb_cleanup(plan);
        return ret;
    }
    return KFA_RET_FAIL_INTRNL;
}

static inline unsigned
scalar_2d_forward(kfft_scalar* in, kfft_sclr2_t* plan, state_t* st) {
    unsigned ret = KFA_RET_FAIL_INTRNL;

    KCB_TYPE(eval2_scalar) eval = KFFT_CALLBACK(st, eval2_scalar);
    KCB_TYPE(shift2_cpx) shift = KFFT_CALLBACK(st, shift2_cpx);

    kfft_cpx* out = KRNL_FUNCS(st).cb_malloc(st->out_count * sizeof(kfft_cpx));
    if (out) {
        unsigned ret = eval(plan, in, out);

        if (ret == KFFT_RET_SUCCESS) {
            if (KFA_CHECK(st, SHIFT))
                shift(out, NULL, st->dims.x, st->dims.y, false, KFFT_PLAN_MMGR(plan));
            write_stdout((kfft_scalar*)out, st);
            ret = KFA_RET_SUCCESS;
        } else {
            KRNL_FUNCS(st).cb_strerr(ret);
            ret = KFA_RET_FAIL_INTRNL;
        }
        KRNL_FUNCS(st).cb_free_null((void**)&out);
    }
    return ret;
}
static inline unsigned
scalar_2d_inverse(kfft_scalar* in, kfft_sclr2_t* plan, state_t* st) {
    unsigned ret = KFA_RET_FAIL_INTRNL;

    KCB_TYPE(evali2_scalar) evali = KFFT_CALLBACK(st, evali2_scalar);
    KCB_TYPE(shift2_cpx) shift = KFFT_CALLBACK(st, shift2_cpx);

    kfft_scalar* out = KRNL_FUNCS(st).cb_malloc(st->out_count * sizeof(kfft_scalar));
    if (out) {
        if (KFA_CHECK(st, SHIFT))
            shift((kfft_cpx*)in, NULL, st->dims.x, st->dims.y, true, KFFT_PLAN_MMGR(plan));
        unsigned ret = evali(plan, (kfft_cpx*)in, out);

        if (ret == KFFT_RET_SUCCESS) {
            write_stdout(out, st);
            ret = KFA_RET_SUCCESS;
        } else {
            KRNL_FUNCS(st).cb_strerr(ret);
            ret = KFA_RET_FAIL_INTRNL;
        }
        KRNL_FUNCS(st).cb_free_null((void**)&out);
    }
    return ret;
}
static unsigned
work_scalar_2d(kfft_scalar* in, state_t* st) {
    KCB_TYPE(config2_scalar) config = KFFT_CALLBACK(st, config2_scalar);
    if (config == NULL)
        return KFA_RET_FAIL_UNREAL;

    kfft_sclr2_t* plan = config(st->dims.x, st->dims.y, sm2kfl(st), NULL, NULL);
    if (plan) {
        unsigned ret = (KFA_CHECK(st, INVERSE)) ? scalar_2d_inverse(in, plan, st)
                                                : scalar_2d_forward(in, plan, st);
        KRNL_FUNCS(st).cb_cleanup(plan);
        return ret;
    }
    return KFA_RET_FAIL_INTRNL;
}
static inline unsigned
scalar_normal_forward(kfft_scalar* in, kfft_sclr_t* plan, state_t* st) {
    unsigned ret = KFA_RET_FAIL_INTRNL;

    KCB_TYPE(eval_scalar) eval = KFFT_CALLBACK(st, eval_scalar);
    KCB_TYPE(shift_cpx) shift = KFFT_CALLBACK(st, shift_cpx);

    kfft_cpx* out = KRNL_FUNCS(st).cb_malloc(st->out_count * sizeof(kfft_cpx));
    if (out) {
        unsigned ret = eval(plan, in, out);

        if (ret == KFFT_RET_SUCCESS) {
            if (KFA_CHECK(st, SHIFT))
                shift(out, st->out_count, false, KFFT_PLAN_MMGR(plan));
            write_stdout((kfft_scalar*)out, st);
            ret = KFA_RET_SUCCESS;
        } else {
            KRNL_FUNCS(st).cb_strerr(ret);
            ret = KFA_RET_FAIL_INTRNL;
        }
        KRNL_FUNCS(st).cb_free_null((void**)&out);
    }
    return ret;
}
static inline unsigned
scalar_normal_inverse(kfft_scalar* in, kfft_sclr_t* plan, state_t* st) {
    unsigned ret = KFA_RET_FAIL_INTRNL;

    KCB_TYPE(evali_scalar) evali = KFFT_CALLBACK(st, evali_scalar);
    KCB_TYPE(shift_cpx) shift = KFFT_CALLBACK(st, shift_cpx);

    kfft_scalar* out = KRNL_FUNCS(st).cb_malloc(st->out_count * sizeof(kfft_scalar));
    if (out) {
        if (KFA_CHECK(st, SHIFT))
            shift((kfft_cpx*)in, st->in_count, true, KFFT_PLAN_MMGR(plan));

        unsigned ret = evali(plan, (kfft_cpx*)in, out);

        if (ret == KFFT_RET_SUCCESS) {
            write_stdout(out, st);
            ret = KFA_RET_SUCCESS;
        } else {
            KRNL_FUNCS(st).cb_strerr(ret);
            ret = KFA_RET_FAIL_INTRNL;
        }
        KRNL_FUNCS(st).cb_free_null((void**)&out);
    }
    return ret;
}

static unsigned
work_scalar_normal(kfft_scalar* in, state_t* st) {
    KCB_TYPE(config_scalar) config = KFFT_CALLBACK(st, config_scalar);
    if (config == NULL)
        return KFA_RET_FAIL_UNREAL;

    kfft_sclr_t* plan = config(st->in_count, sm2kfl(st), NULL, NULL);
    if (plan) {
        unsigned ret = (KFA_CHECK(st, INVERSE)) ? scalar_normal_inverse(in, plan, st)
                                                : scalar_normal_forward(in, plan, st);
        KRNL_FUNCS(st).cb_cleanup(plan);
        return ret;
    }
    return KFA_RET_FAIL_INTRNL;
}

unsigned
work_scalar_plan(kfft_scalar* in, state_t* st) {
    unsigned ret = KFA_RET_FAIL_INTRNL;

    if (KFA_CHECK(st, SPARSE)) {
        ret = work_scalar_sparse(in, st);
    } else {
        if (KFA_CHECK(st, 2D)) {
            ret = work_scalar_2d(in, st);
        } else {
            ret = work_scalar_normal(in, st);
        }
    }

    return ret;
}
