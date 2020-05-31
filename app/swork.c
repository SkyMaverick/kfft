#if defined(KFFT_2D_ENABLE)

static kfft_return_t
work_scalar2_forward(kfft_scalar* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* ftmp = kfft_malloc(M->len * sizeof(kfft_cpx));
    if (buf && ftmp) {
        kfft_sclr2_t* plan = kfft_config2_scalar(M->x, M->y, M->flags, 0, NULL);
        if (plan) {
            ret = kfft_eval2_scalar(plan, buf, ftmp);
            if (ret == KFFT_RET_SUCCESS) {
                if (M->is_shift)
                    kfft_shift2_cpx(ftmp, NULL, M->x, M->y, false, KFFT_PLAN_MMGR(plan));
            }
            kfft_cleanup(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        }
        if (ret == KFFT_RET_SUCCESS) {
            write_stdout((kfft_scalar*)ftmp, M->len * 2);
            fprintf(stdout, "%s\n", "");
        }
        kfft_free(&ftmp);
    }
    return ret;
}
static kfft_return_t
work_scalar2_inverse(kfft_cpx* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_scalar* ftmp = kfft_malloc(M->len * sizeof(kfft_scalar));
    if (buf && ftmp) {
        kfft_sclr2_t* plan = kfft_config2_scalar(M->x, M->y, M->flags, 0, NULL);
        if (plan) {
            if (M->is_shift)
                kfft_shift2_cpx(buf, NULL, M->x, M->y, true, KFFT_PLAN_MMGR(plan));

            ret = kfft_evali2_scalar(plan, buf, ftmp);
            kfft_cleanup(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        }
        if (ret == KFFT_RET_SUCCESS) {
            write_stdout((kfft_scalar*)ftmp, M->len);
            fprintf(stdout, "%s\n", "");
        }
        kfft_free(&ftmp);
    }
    return ret;
}
#endif /* KFFT_2D_ENABLE */
#if defined(KFFT_SPARSE_ENABLE)
static kfft_return_t
work_scalar_sparse_forward(kfft_scalar* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* ftmp = kfft_malloc(M->len * sizeof(kfft_cpx));
    if (buf && ftmp) {
        kfft_ssparse_t* plan =
            kfft_config_sparse_scalar(M->len, M->flags, M->dim, M->step, 0, NULL);
        if (plan) {
            ret = kfft_eval_sparse_scalar(plan, buf, ftmp);
            if (ret == KFFT_RET_SUCCESS) {
                if (M->is_shift)
                    kfft_shift_sparse_cpx(ftmp, NULL, M->len, M->dim, M->step, true,
                                          KFFT_PLAN_MMGR(plan));
            }
            kfft_cleanup(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        } /* plan != NULL */

        if (ret == KFFT_RET_SUCCESS) {
            write_stdout((kfft_scalar*)ftmp, M->len * 2);
            fprintf(stdout, "%s\n", "");
        }
        kfft_free(&ftmp);
    } /* in && ftmp */
    return ret;
}

static kfft_return_t
work_scalar_sparse_inverse(kfft_cpx* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_scalar* ftmp = kfft_malloc(M->len * sizeof(kfft_scalar));
    if (buf && ftmp) {
        kfft_ssparse_t* plan =
            kfft_config_sparse_scalar(M->len, M->flags, M->dim, M->step, 0, NULL);
        if (plan) {
            if (M->is_shift)
                kfft_shift_sparse_cpx(buf, NULL, M->len, M->dim, M->step, true,
                                      KFFT_PLAN_MMGR(plan));

            ret = kfft_evali_sparse_scalar(plan, buf, ftmp);
            kfft_cleanup(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        } /* plan != NULL */

        if (ret == KFFT_RET_SUCCESS) {
            write_stdout(ftmp, M->len);
            fprintf(stdout, "%s\n", "");
        }
        kfft_free(&ftmp);
    } /* in && ftmp */
    return ret;
}
#endif /* KFFT_SPARSE_ENABLE */

static kfft_return_t
work_scalar_forward(kfft_scalar* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* ftmp = kfft_malloc(M->len * sizeof(kfft_cpx));
    if (buf && ftmp) {
        kfft_sclr_t* plan = kfft_config_scalar(M->len, M->flags, 0, NULL);
        if (plan) {
            ret = kfft_eval_scalar(plan, buf, ftmp);
            if (ret == KFFT_RET_SUCCESS) {
                if (M->is_shift)
                    kfft_shift_cpx(ftmp, M->len, false, KFFT_PLAN_MMGR(plan));
            }
            kfft_cleanup(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        } /* plan != NULL */

        if (ret == KFFT_RET_SUCCESS) {
            write_stdout((kfft_scalar*)ftmp, M->len * 2);
            fprintf(stdout, "%s\n", "");
        }
        kfft_free(&ftmp);
    } /* in && ftmp */
    return ret;
}

static kfft_return_t
work_scalar_inverse(kfft_cpx* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_scalar* ftmp = kfft_malloc(M->len * sizeof(kfft_scalar));
    if (buf && ftmp) {
        kfft_sclr_t* plan = kfft_config_scalar(M->len, M->flags, 0, NULL);
        if (plan) {
            if (M->is_shift)
                kfft_shift_cpx(buf, M->len, true, KFFT_PLAN_MMGR(plan));

            ret = kfft_evali_scalar(plan, buf, ftmp);
            kfft_cleanup(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        } /* plan != NULL */

        if (ret == KFFT_RET_SUCCESS) {
            write_stdout(ftmp, M->len);
            fprintf(stdout, "%s\n", "");
        }
        kfft_free(&ftmp);
    } /* in && ftmp */
    return ret;
}

static inline kfft_return_t
work_scalar(char* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    void* fin = NULL;
    M->len = parse_buffer(&fin, buf, (M->flags & KFFT_FLAG_INVERSE) ? true : false);
    if (fin != NULL) {
        if (M->len > 0) {
#if defined(KFFT_2D_ENABLE)
            if (M->is_2d) {
                if (prepare_2d(M) == 0) {
                    ret = (M->flags & KFFT_FLAG_INVERSE)
                              ? work_scalar2_inverse((kfft_cpx*)fin, M)
                              : work_scalar2_forward((kfft_scalar*)fin, M);
                } else
                    ret = KFFT_RET_BADARGUMENTS;
            } else
#endif /* KFFT_2D_ENABLE */
                if (M->is_sparse) {
                ret = (M->flags & KFFT_FLAG_INVERSE)
                          ? work_scalar_sparse_inverse((kfft_cpx*)fin, M)
                          : work_scalar_sparse_forward((kfft_scalar*)fin, M);
            } else {
                ret = (M->flags & KFFT_FLAG_INVERSE) ? work_scalar_inverse((kfft_cpx*)fin, M)
                                                     : work_scalar_forward((kfft_scalar*)fin, M);
            }
        }
        kfft_free(&fin);
    }
    return ret;
}
