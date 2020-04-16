#if defined(KFFT_2D_ENABLE)

static kfft_return_t
work_cpx2_internal(kfft_cpx* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_ALLOC_FAIL;

    kfft_cpx* ftmp = calloc(M->len, sizeof(kfft_cpx));
    if (buf && ftmp) {
        // TODO
        //        if ((!(M->flags & KFFT_FLAG_INVERSE)) && (M->is_shift))
        //            kfft_shift_cpx(ftmp, M->len, true);

        kfft_comp2_t* plan = kfft_config2_cpx(M->x, M->y, M->flags, 0, NULL);
        if (plan) {
            ret = kfft_eval2_cpx(plan, buf, ftmp);
            kfft_free(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        } /* plan != NULL */

        if (ret == KFFT_RET_SUCCESS) {
            // TODO
            //            if ((M->flags & KFFT_FLAG_INVERSE) && (M->is_shift))
            //                kfft_shift_cpx(ftmp, M->len, true);

            write_stdout((kfft_scalar*)ftmp, M->len * 2);
            fprintf(stdout, "%s\n", "");
        }
        free(ftmp);
    } /* in && ftmp */
    return ret;
}
#endif /* KFFT_2D_ENABLE */

static kfft_return_t
work_cpx_internal(kfft_cpx* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_ALLOC_FAIL;

    kfft_cpx* ftmp = calloc(M->len, sizeof(kfft_cpx));
    if (buf && ftmp) {
        if ((M->flags & KFFT_FLAG_INVERSE) && (M->is_shift))
            kfft_shift_cpx(buf, M->len, true);

        kfft_comp_t* plan = kfft_config_cpx(M->len, M->flags, 0, 0, NULL);
        if (plan) {
            ret = kfft_eval_cpx(plan, buf, ftmp);
            kfft_free(plan);
        } else {
            ret = KFFT_RET_ALLOC_FAIL;
        } /* plan != NULL */

        if (ret == KFFT_RET_SUCCESS) {
            if ((!(M->flags & KFFT_FLAG_INVERSE)) && (M->is_shift))
                kfft_shift_cpx(ftmp, M->len, false);

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
#if defined(KFFT_2D_ENABLE)
            if (M->is_2d) {
                if (prepare_2d(M) == 0) {
                    ret = work_cpx2_internal(fin, M);
                } else {
                    ret = KFFT_RET_BADARGUMENTS;
                }
            } else
#endif /* KFFT_2D_ENABLE */
                ret = work_cpx_internal(fin, M);
        }
        free(fin);
    }
    return ret;
}
