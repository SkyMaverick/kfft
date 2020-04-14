static inline kfft_return_t
work_scalar(char* buf, app_mode_t* M) {

#define _VC(X) ((kfft_cpx*)(X))
#define _VS(X) ((kfft_scalar*)(X))

    kfft_return_t ret = KFFT_RET_SUCCESS;

    void* fin = NULL;
    ssize_t nfft = parse_buffer(&fin, buf, (M->flags & KFFT_FLAG_INVERSE) ? true : false);
    if (nfft > 0) {
        void* ftmp = (M->flags & KFFT_FLAG_INVERSE) ? calloc(nfft, sizeof(kfft_scalar))
                                                    : calloc(nfft, sizeof(kfft_cpx));
        if (fin && ftmp) {
            if (M->is_2d) {

            } else {
                kfft_sclr_t* plan = kfft_config_scalar(nfft, M->flags, 0, NULL);
                if (plan) {
                    ret = (M->flags & KFFT_FLAG_INVERSE)
                              ? kfft_evali_scalar(plan, _VC(fin), _VS(ftmp))
                              : kfft_eval_scalar(plan, _VS(fin), _VC(ftmp));
                    kfft_free(plan);
                } else {
                    ret = KFFT_RET_ALLOC_FAIL;
                } /* plan != NULL */
            }     /* is_2d */
            if (ret == KFFT_RET_SUCCESS) {
                if (M->flags & KFFT_FLAG_INVERSE) {
                    write_stdout(_VS(ftmp), nfft);
                } else {
                    write_stdout(_VS(ftmp), nfft * 2);
                }
                fprintf(stdout, "%s\n", "");
            }
            free(ftmp);
        } /* in && ftmp */
    }
    return ret;
#undef _VC
#undef _VS
}
