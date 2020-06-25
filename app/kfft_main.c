#include "kfa_defs.h"

#include "loader.c"
#include "cmdline.c"

static inline void
write_stdout(kfft_scalar* in, state_t* st) {
    for (size_t i = 0; i < st->out_len; i++) {
        (i > 0) ? fprintf(stdout, " %.3f", in[i]) : fprintf(stdout, "%.3f", in[i]);
    }
    fprintf(stdout, "%s\n", "");
    fflush(stdout);
}

#include "complex.c"
#include "scalar.c"

int
main(int argc, char* argv[]) {
    unsigned ret = KFA_RET_SUCCESS;

    state_t* k_state = calloc(1, sizeof(state_t));
    if (k_state == NULL)
        return -1;

    ret = load_kfft_core(k_state);
    if (ret == KFA_RET_SUCCESS) {

        kfft_scalar* buffer = cmd_line_parse(argc, argv, k_state);
        if (buffer) {
            ret = (k_state->mode & KFA_MODE_SCALAR) ? work_scalar_plan(buffer, k_state)
                                                    : work_complex_plan(buffer, k_state);
            KRNL_FUNCS(k_state).cb_free_null((void*)(&buffer));
        } else {
            display_help();
        }
        unload_kfft_core(k_state);
    }
    free(k_state);
    return ret;
}
