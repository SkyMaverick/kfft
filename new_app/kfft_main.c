#include "kfa_defs.h"

#include "loader.c"
#include "cmdline.c"

static inline void
write_stdout_binary(kfft_scalar* in, state_t* st) {
    // TODO Copy buffer
    fprintf(stdout, "%s\n", (char*)in);
    fflush(stdout);
}

static inline void
write_stdout_manual(kfft_scalar* in, state_t* st) {
    char buf[80];
    size_t out_size = 1;

    char* out = calloc(1, sizeof(STDOUT_BUF_SIZE));
    if (out) {
        out[0] = '\0';
        for (size_t i = 0; i < st->out_lenght; i++) {
            sprintf(buf, "%.3f ", in[i]);
            out_size += strlen(buf);
            char* old = out;
            out = realloc(out, out_size);
            if (out == NULL) {
                free(old);
            }
            strcat(out, buf);
        }
        out[out_size - 2] = '\0'; // erase last space

        fprintf(stdout, "%s\n", out);
        fflush(stdout);

        free(out);
    } /* out allocated */
}

void
write_stdout(kfft_scalar* in, state_t* st) {
    (st->mode & KFA_MODE_BINARY) ? write_stdout_binary(in, st) : write_stdout_manual(in, st);
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
        //        unload_kfft_core(k_state);
    }
    free(k_state);
    return ret;
}
