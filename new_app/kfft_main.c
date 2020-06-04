#include "kfa_defs.h"

#include "loader.c"
#include "cmdline.c"

int
main(int argc, char* argv[]) {
    unsigned ret = KFA_RET_SUCCESS;

    state_t* k_state = calloc(1, sizeof(state_t));
    if (k_state == NULL)
        return -1;

    ret = load_kfft_core(k_state);
    if (ret == KFA_RET_SUCCESS) {

        char* buffer = cmd_line_parse(argc, argv, k_state);

        unload_kfft_core(k_state);
    }
    free(k_state);
    return ret;
}
