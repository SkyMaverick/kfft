#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kfft.h"
#if defined(KFFT_OS_WINDOWS)
    #include <Windows.h>
    #include "getopt_win.h"
#else
    #include <dlfcn.h>
    #include <unistd.h>
    #include <getopt.h>
#endif

#if defined(KFFT_OS_WINDOWS)
    #define KFFT_LIBRARY_NAME "kfft.dll"
#else
    #define KFFT_LIBRARY_NAME "libkfft.so"
#endif

#include "const.h"
#include "config.h"

#define KCB_TYPE(X) kfft_callback_##X

typedef struct {
    KCB_TYPE(info) cb_info;
    KCB_TYPE(next_fast_size) cb_next_fast_size;
    KCB_TYPE(malloc) cb_malloc;
    KCB_TYPE(free_null) cb_free_null;
    KCB_TYPE(cleanup) cb_cleanup;
    KCB_TYPE(strerr) cb_strerr;
} vtm_t;

enum {
    KFA_RET_SUCCESS = 0,
    KFA_RET_FAIL_INTRNL = 1,
    KFA_RET_FAIL_PARSE = 2,
    KFA_RET_FAIL_LOAD = 3,
    KFA_RET_FAIL_ARGS = 4,
    KFA_RET_FAIL_UNREAL = 5,
};

enum {
    KFA_MODE_SCALAR = 1 << 0,
    KFA_MODE_INVERSE = 1 << 2,
    KFA_MODE_GENERIC = 1 << 3,
    KFA_MODE_GENONLY = 1 << 4,
    KFA_MODE_SHIFT = 1 << 5,
    KFA_MODE_2D = 1 << 6,
    KFA_MODE_SPARSE = 1 << 7,
    KFA_MODE_STDIN = 1 << 8,
    KFA_MODE_BINARY = 1 << 9,
};

typedef struct {
    struct {
#if defined(KFFT_OS_WINDOWS)
        HMODULE handle;
#else
        void* handle;
#endif
        vtm_t sfuns;
    } lib;
    struct {
        size_t lenght;
    } buf;
    struct {
        uint32_t x;
        uint32_t y;
    } dims;
    struct {
        uint32_t dx;
        uint32_t sx;
    } sparse;

    size_t in_len;
    size_t in_count;
    size_t out_len;
    size_t out_count;

    uint8_t scsz;
    uint32_t mode;
} state_t;

#define KFA_CHECK(S, X) (S)->mode& KFA_MODE_##X

#if defined(KFFT_OS_WINDOWS)
    #define KFFT_CALLBACK(S, X) (KCB_TYPE(X)) GetProcAddress((S)->lib.handle, "kfft_" #X)
#else
    #define KFFT_CALLBACK(S, X) (KCB_TYPE(X)) dlsym((S)->lib.handle, "kfft_" #X)
#endif /* KFFT_OS_WINDOWS */

#define KRNL_FUNCS(S) (S)->lib.sfuns
#define KRNL_LIB(S) (S)->lib.handle

static inline uint32_t
sm2kfl(state_t* st) {
    uint32_t ret = KFFT_FLAG_NORMAL;

    if (st->mode & KFA_MODE_INVERSE)
        ret |= KFFT_FLAG_INVERSE;
    if (st->mode & KFA_MODE_GENERIC)
        ret |= KFFT_FLAG_GENERIC;
    if (st->mode & KFA_MODE_GENONLY)
        ret |= KFFT_FLAG_GENERIC_ONLY;

    return ret;
}
