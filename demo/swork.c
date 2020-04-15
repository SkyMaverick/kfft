static kfft_return_t
work_scalar_forward(kfft_scalar* buf, app_mode_t* M) {

    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* ftmp = calloc(M->len, sizeof(kfft_cpx));
    if (buf && ftmp) {
        if (M->is_2d) {
            // TODO
        } else {
            if (M->is_shift)
                kfft_shift_scalar(buf, M->len, false);

            kfft_sclr_t* plan = kfft_config_scalar(M->len, M->flags, 0, NULL);
            if (plan) {
                ret = kfft_eval_scalar(plan, buf, ftmp);
                kfft_free(plan);
            } else {
                ret = KFFT_RET_ALLOC_FAIL;
            } /* plan != NULL */
        }     /* is_2d */
        if (ret == KFFT_RET_SUCCESS) {
            write_stdout((kfft_scalar*)ftmp, M->len * 2);
            fprintf(stdout, "%s\n", "");
        }
        free(ftmp);
    } /* in && ftmp */
    return ret;
}

static kfft_return_t
work_scalar_inverse(kfft_cpx* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_scalar* ftmp = calloc(M->len, sizeof(kfft_scalar));
    if (buf && ftmp) {
        if (M->is_2d) {
            // TODO
        } else {
            kfft_sclr_t* plan = kfft_config_scalar(M->len, M->flags, 0, NULL);
            if (plan) {
                ret = kfft_evali_scalar(plan, buf, ftmp);
                kfft_free(plan);
            } else {
                ret = KFFT_RET_ALLOC_FAIL;
            } /* plan != NULL */
        }     /* is_2d */
        if (ret == KFFT_RET_SUCCESS) {
            if (M->is_shift)
                kfft_shift_scalar(ftmp, M->len, true);

            write_stdout(ftmp, M->len);
            fprintf(stdout, "%s\n", "");
        }
        free(ftmp);
    } /* in && ftmp */
    return ret;
}

static int
prepare_2d(app_mode_t* M) {
    if ((M->len > 0) && (M->x > 0) && (M->len % M->x > 0))
        return 1;

    M->y = M->len / M->x;
    return 0;
}

static inline kfft_return_t
work_scalar(char* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    void* fin = NULL;
    M->len = parse_buffer(&fin, buf, (M->flags & KFFT_FLAG_INVERSE) ? true : false);

    if (M->len > 0) {
        if ((M->is_2d) && (prepare_2d(M))) {
            ret = KFFT_RET_BADARGUMENTS;
            goto bailout;
        }

        ret = (M->flags & KFFT_FLAG_INVERSE) ? work_scalar_inverse((kfft_cpx*)fin, M)
                                             : work_scalar_forward((kfft_scalar*)fin, M);
    }

bailout:
    return ret;
}
