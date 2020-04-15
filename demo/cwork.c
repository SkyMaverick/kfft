static kfft_return_t
work_cpx_internal2(kfft_cpx* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    // TODO
    return ret;
}

static kfft_return_t
work_cpx_internal(kfft_cpx* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* ftmp = calloc(M->len, sizeof(kfft_cpx));
    if (buf && ftmp) {
        if ((!(M->flags & KFFT_FLAG_INVERSE)) && (M->is_shift))
            kfft_shift_cpx(ftmp, M->len, true);

        kfft_comp_t* plan = kfft_config_cpx(M->len, M->flags, 0, 0, NULL);
        if (plan) {
            ret = kfft_eval_cpx(plan, buf, ftmp);
            kfft_free(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        } /* plan != NULL */

        if (ret == KFFT_RET_SUCCESS) {
            if ((M->flags & KFFT_FLAG_INVERSE) && (M->is_shift))
                kfft_shift_cpx(ftmp, M->len, true);

            write_stdout((kfft_scalar*)ftmp, M->len * 2);
            fprintf(stdout, "%s\n", "");
        }
        free(ftmp);
    } /* in && ftmp */
    return ret;
}

static inline kfft_return_t
work_cpx(char* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* fin = NULL;
    M->len = parse_buffer((void**)(&fin), buf, true);
    if (fin != NULL) {
        if (M->len > 0) {
            if (M->is_2d) {
                if (prepare_2d(M) == 0) {
                    ret = work_cpx_internal2(fin, M);
                } else {
                    ret = KFFT_RET_BADARGUMENTS;
                }
            } else {
                ret = work_cpx_internal(fin, M);
            }
        }
        free(fin);
    }
    return ret;
}
