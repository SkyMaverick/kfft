#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/ioctl.h> /* FIONREAD ioctl() */

#if defined(KFFT_OS_WINDOWS)
    #include "getopt_win.h"
#else
    #include <getopt.h>
#endif

#include "kfft.h"
#include "const.h"

typedef struct {
    void* buf;

    bool is_shift;
    bool is_cpx;
    bool is_2d;
    bool is_stdin;

    size_t len;
    size_t x;
    size_t y;

    uint32_t flags;
} app_mode_t;

static void
display_info(void) {
    kfft_info_t info;
    kfft_info(&info);

    fprintf(stdout, "LibKFFT version : %d.%d.%d\n\n", info.vmajor, info.vminor, info.vpatch);
    fprintf(stdout, "%s - %s\n", "Enable trace messages",
            (info.flags & KFFT_INFO_TRACE) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use SIMD instructions",
            (info.flags & KFFT_INFO_USE_SIMD) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use alloca() function",
            (info.flags & KFFT_INFO_USE_ALLOCA) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use system math functions",
            (info.flags & KFFT_INFO_USE_SYSMATH) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use OpenMP functions",
            (info.flags & KFFT_INFO_USE_OPENMP) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable use Rader algoritm",
            (info.flags & KFFT_INFO_RADER_ALGO) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable lesser memory mode",
            (info.flags & KFFT_INFO_MEMLESS_MODE) ? "YES" : "NO");
    fprintf(stdout, "%s - %s\n", "Enable half scalar mode",
            (info.flags & KFFT_INFO_HALF_SCALAR) ? "YES" : "NO");
}

static void
display_help(void) {
    kfft_info_t info;
    kfft_info(&info);

    fprintf(stdout, "Use libkfft version : %d.%d.%d\n\n", info.vmajor, info.vminor, info.vpatch);

    fprintf(stdout, "%s", help_msg);
}

static size_t
parse_buffer(void** out, char* buf, bool as_cpx) {
    char* args = buf;

    size_t len = 1;
    // Analize
    while ((args = strchr(args, ' ')) != NULL)
        len++, *args = '\0', args++;

    args = buf;

    kfft_scalar* tmp = calloc(len, sizeof(kfft_cpx));
    if (tmp) {
        for (size_t i = 0; i < len; i++) {
            tmp[i] = (kfft_scalar)atof(args);
            args += strlen(args) + 1;
        }
    }
    *out = tmp;
    return (as_cpx) ? len / 2 : len;
}

int
stdin_check(void) {
    fd_set rd;
    struct timeval tv = {1, 0};
    int ret;

    FD_ZERO(&rd);
    FD_SET(STDIN_FILENO, &rd);
    ret = select(1, &rd, NULL, NULL, &tv);

    return (ret > 0);
}

static char*
read_stdin(void) {
    char buf[STDIN_BUF_SIZE];
    size_t ret_size = 1;
    char* ret = NULL;

    if (stdin_check() > 0) {
        ret = malloc(STDIN_BUF_SIZE * sizeof(char));
        if (ret) {
            ret[0] = '\0';
            size_t n = 0;
            while ((n = fread(buf, 1, STDIN_BUF_SIZE, stdin)) > 0) {
                char* old = ret;
                ret_size += n;
                ret = realloc(ret, ret_size + 1);
                if (ret == NULL) {
                    free(old);
                    return NULL;
                };
                strcat(ret, buf);
            }
            ret[ret_size] = '\0';
        } /* ret allocated */
    }
    return ret;
}

static void
write_stdout(kfft_scalar* in, size_t sz) {
    char buf[80];
    size_t out_size = 1;

    char* out = calloc(1, sizeof(STDOUT_BUF_SIZE));
    if (out) {
        out[0] = '\0';
        for (size_t i = 0; i < sz; i++) {
            sprintf(buf, "%.3f ", in[i]);
            out_size += strlen(buf);
            char* old = out;
            out = realloc(out, out_size);
            if (out == NULL) {
                free(old);
            }
            strcat(out, buf);
        }
        out[out_size] = '\0';

        fprintf(stdout, "%s", out);
        fflush(stdout);

        free(out);
    } /* out allocated */
}

static char*
cmd_line_parse(int argc, char* argv[], app_mode_t* mode) {
    kfft_info_t info;
    char* ret = NULL;

    int opt = 0;
    while ((opt = getopt(argc, argv, "gGidsSx:f:vV?")) != -1) {
        switch (opt) {
        case 'g':
            mode->flags |= KFFT_FLAG_GENERIC;
            break;
        case 'G':
            mode->flags |= KFFT_FLAG_GENERIC_ONLY;
            break;
        case 'i':
            mode->flags |= KFFT_FLAG_INVERSE;
            break;
        case 's':
            mode->is_shift = true;
            break;
        case 'S':
            mode->is_cpx = false;
            break;
        case '?':
            display_help();
            exit(0);
        case 'v':
            kfft_info(&info);
            fprintf(stdout, "%d.%d.%d\n", info.vmajor, info.vminor, info.vpatch);
            exit(0);
        case 'V':
            display_info();
            exit(0);
        }
    }
    if (mode->is_stdin) {
        ret = read_stdin();
    }
    if (ret == NULL) {
        if (optind < argc) {
            size_t tmp_size = 1;
            char* buf = calloc(1, tmp_size);
            if (buf) {
                do {
                    char* old = buf;
                    tmp_size += strlen(argv[optind]) + 1;

                    buf = realloc(buf, tmp_size);
                    if (buf == NULL) {
                        free(old);
                        goto bailout;
                    }

                    strcat(buf, argv[optind]);
                    strcat(buf, " ");
                } while (argc > ++optind);

                buf[strlen(buf) - 1] = '\0';
            } /* buf allocated */
            ret = buf;
        }
    } /* ret == NULL */
bailout:
    return ret;
}

static kfft_return_t
work_cpx(char* buf, app_mode_t* M) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    if (M->is_2d) {
    } else {
    }
    return ret;
}

static kfft_return_t
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

int
main(int argc, char* argv[]) {
    // clang-format off
    app_mode_t mode = {
            .buf = NULL,
            .is_shift = true,
            .is_cpx = true,
            .is_2d = false,
            .is_stdin = true,
            
            .len = 0,
            .x = 0,
            .y = 0,
            
            .flags = 0};
    // clang-format on

    kfft_return_t ret = KFFT_RET_SUCCESS;
    char* buffer = cmd_line_parse(argc, argv, &mode);

    if (buffer) {
        ret = (mode.is_cpx) ? work_cpx(buffer, &mode) : work_scalar(buffer, &mode);
        free(buffer);
    } else {
        display_help();
    }
    fflush(stdout);
    fclose(stdout);

    return ret;
}
