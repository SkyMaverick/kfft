#include "util.h"

#include <CUnit/Basic.h>
#include <dirent.h>
#include <dlfcn.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define PATH_MAX 4096

typedef void (*cb_suite)(void);

typedef struct modules_list {
    void* handle;
    cb_suite start_func;

    char* module_path;
    struct modules_list* next;
} module_t;

size_t
lookup_suites(char* lookup_dir, module_t** found_suites) {
    size_t suites_count = 0;
    DIR* suites_path = opendir(lookup_dir);
    struct dirent* ls = NULL;
    cb_suite func_suite = NULL;

    if (!suites_path) {
        fprintf(stderr, "%s error: %s\n", lookup_dir, strerror(errno));
        return -1;
    }

    module_t* tmp_module = NULL;
    while ((ls = readdir(suites_path))) {
        if ((strncmp(ls->d_name, ".", 1) == 0) || (strstr(ls->d_name, ".so") == 0))
            continue;

        char buf_path[PATH_MAX];
        snprintf(buf_path, PATH_MAX, "%s/%s", lookup_dir, ls->d_name);

        void* h_module = dlopen(buf_path, RTLD_LAZY);
        char* error = NULL;
        error = dlerror();
        if (!h_module) {
            fprintf(stderr, "Couldn't load module: %s\n", error);
            continue;
        }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic" /*ignore dlsym (void*) casting warning */
        func_suite = dlsym(h_module, "_run_suite");
#pragma GCC diagnostic pop

        error = dlerror();
        if (error) {
            fprintf(stderr, "Couldn't load module: %s\n", error);
            dlclose(h_module);
            continue;
        }

        module_t* new_module = (module_t*)malloc(sizeof(module_t));
        memset(new_module, 0, sizeof(module_t));
        if (tmp_module) {
            tmp_module->next = new_module;
            tmp_module = new_module;
        } else {
            tmp_module = new_module;
            *found_suites = new_module;
        }

        fprintf(stdout, "Found test suite: %s\n", ls->d_name);

        new_module->handle = h_module;
        new_module->start_func = func_suite;
        new_module->module_path = strdup(buf_path);

        suites_count++;
    }

    closedir(suites_path);

    return suites_count;
}

int
main(int argc, char* argv[]) {
    (void)argc;

    char app_path[PATH_MAX];
    char suites_path[PATH_MAX];

    if (!realpath(argv[0], app_path)) {
        snprintf(app_path, PATH_MAX, "%s", argv[0]);
    }
    char* e = strrchr(app_path, '/');
    if (e)
        *e = 0;

    snprintf(suites_path, PATH_MAX, "%s/%s", app_path, "suites");

    struct stat d_info;
    module_t* suites_list = NULL;

    if (stat(suites_path, &d_info) < 0) {
        if (!S_ISDIR(d_info.st_mode)) {
            fprintf(stderr, "%s is not directory\n", suites_path);
            return 1;
        }
    }

    if (access(suites_path, R_OK | X_OK)) {
        fprintf(stderr, "%s is not accessible\n", suites_path);
        return 1;
    }

    size_t suites_count = lookup_suites(suites_path, &suites_list);
    if (suites_count) {
        fprintf(stdout, "Found %zu test suites\n", suites_count);
    } else {
        fprintf(stdout, "Suites modules not found\n");
        return 0;
    }

    CUnitInitialize();

    for (module_t* suite = suites_list; suite; suite = suite->next) {
        suite->start_func();
    }

    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CUnitUInitialize();

    // Cleanup
    while (suites_list != NULL) {
        module_t* tmp_suites = suites_list;
        if (suites_list->next) {
            suites_list = suites_list->next;
        } else
            suites_list = NULL;

        fprintf(stdout, "Free module: %s\n", tmp_suites->module_path);

        dlclose(tmp_suites->handle);
        free(tmp_suites->module_path);
        free(tmp_suites);
    }
    return CU_get_error();
}
