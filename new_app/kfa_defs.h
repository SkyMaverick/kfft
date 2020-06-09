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

typedef struct {
    kfft_callback_info cb_info;
    kfft_callback_next_fast_size cb_next_fast_size;
    kfft_callback_malloc cb_malloc;
    kfft_callback_free_null cb_free_null;
    kfft_callback_cleanup cb_cleanup;
    kfft_callback_strerr cb_strerr;
} vtm_t;

enum {
    KFA_RET_SUCCESS = 0,
    KFA_RET_FAIL_INTERNAL = 1,
    KFA_RET_FAIL_PARSE = 2,
    KFA_RET_FAIL_LOAD = 3,
    KFA_RET_FAIL_ARGS = 4
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

    size_t lenght;
    size_t out_lenght;

    uint32_t mode;
} state_t;

#if defined(KFFT_OS_WINDOWS)
    #define KFFT_CALLBACK(S, X) GetProcAddress((S)->lib.handle, "kfft_" X)
#else
    #define KFFT_CALLBACK(S, X) dlsym((S)->lib.handle, "kfft_" X)
#endif /* KFFT_OS_WINDOWS */

#define KRNL_FUNCS(S) (S)->lib.sfuns
#define KRNL_LIB(S) (S)->lib.handle
