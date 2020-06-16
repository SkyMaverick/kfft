#include "kfa_defs.h"

#if defined(KFFT_OS_WINDOWS)
static inline unsigned
load_kfft_library(state_t* st) {
    st->lib.handle = LoadLibraryA(KFFT_LIBRARY_NAME);
    return (st->lib.handle) ? KFA_RET_SUCCESS : KFA_RET_FAIL_LOAD;
}

static inline void
unload_kfft_library(state_t* st) {
    if (st->lib.handle)
        FreeLibrary(st->lib.handle);
}
#else
static inline unsigned
load_kfft_library(state_t* st) {
    st->lib.handle = dlopen(KFFT_LIBRARY_NAME, RTLD_LAZY);
    return (st->lib.handle) ? KFA_RET_SUCCESS : KFA_RET_FAIL_LOAD;
}
static inline void
unload_kfft_library(state_t* st) {
    if (st->lib.handle)
        dlclose(st->lib.handle);
}
#endif

unsigned
read_core_info(state_t* st) {
    kfft_info_t info = {0};
    KRNL_FUNCS(st).cb_info(&info);

    // TODO Add new info if it's need
    st->scsz = (info.flags & KFFT_INFO_HALF_SCALAR) ? sizeof(float) : sizeof(double);

    return KFA_RET_SUCCESS;
}

unsigned
load_kfft_core(state_t* st) {
    unsigned ret = load_kfft_library(st);
    if (ret != KFFT_RET_SUCCESS)
        return ret;

    KRNL_FUNCS(st).cb_info = KFFT_CALLBACK(st, info);
    KRNL_FUNCS(st).cb_next_fast_size = KFFT_CALLBACK(st, next_fast_size);
    KRNL_FUNCS(st).cb_malloc = KFFT_CALLBACK(st, malloc);
    KRNL_FUNCS(st).cb_free_null = KFFT_CALLBACK(st, free_null);
    KRNL_FUNCS(st).cb_cleanup = KFFT_CALLBACK(st, cleanup);
    KRNL_FUNCS(st).cb_strerr = KFFT_CALLBACK(st, strerr);

    return read_core_info(st);
}

void
unload_kfft_core(state_t* st) {
    unload_kfft_library(st);
}
