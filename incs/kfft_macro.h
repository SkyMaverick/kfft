#pragma once

// *** OPERATIONS SYSTEM MACRO ****************************************

#if defined(__WIN32__) || defined(_WIN32) || defined(WIN32)
    #define KFFT_OS_WIN32
#elif defined(_WIN64) || defined(_M_X64) || defined(_M_AMD64)
    #define KFFT_OS_WIN64
#elif defined(__linux) || defined(__linux__)
    #define KFFT_OS_LINUX
#elif defined(__ANDROID__)
    #define KFFT_OS_ANDROID
#elif defined(__MSYS__)
    #define KFFT_OS_MSYS
#elif defined(__CYGWIN__)
    #define KFFT_OS_CYGWIN
#elif defined(__APPLE__) && (defined(__GNUC__) || defined(__xlC__) || defined(__xlc__))
    #define KFFT_OS_DARWIN
#elif defined(__FreeBSD__)
    #define KFFT_OS_FREEBSD
    #define KFFT_OS_BSD4
#elif defined(__DragonFly__)
    #define KFFT_OS_DRAGONFLY
    #define KFFT_OS_BSD4
#elif defined(__NetBSD__)
    #define KFFT_OS_NETBSD
    #define KFFT_OS_BSD4
#elif defined(__OpenBSD__)
    #define KFFT_OS_OPENBSD
    #define KFFT_OS_BSD4
#endif

#if defined(KFFT_OS_MSYS) || defined(KFFT_OS_CYGWIN)
    #define KFFT_OS_WINEMULATOR
#endif

#if defined(KFFT_OS_WIN32) || defined(KFFT_OS_WIN64)
    #define KFFT_OS_WINDOWS
#endif

#if defined(KFFT_OS_BSD4) || defined(KFFT_OS_DARWIN) || defined(KFFT_OS_LINUX) ||                  \
    defined(KFFT_OS_ANDROID) || defined(KFFT_OS_CYGWIN) || defined(KFFT_OS_MSYS)

    #define KFFT_OS_POSIX
#endif

// *** COMPILATOR MACRO ****************************************

#if defined(_MSC_VER)
    #define KFFT_CC_MSVC
    #if defined(__INTEL_COMPILER)
        #define KFFT_CC_INTEL
    #endif
    #if defined(__clang__)
        #define KFFT_CC_CLANG
    #endif
#elif defined(__GNUC__)
    #define KFFT_CC_GNU
    #if defined(__MINGW32__)
        #define KFFT_CC_MINGW
    #endif
    #if defined(__INTEL_COMPILER)
        #define KFFT_CC_INTEL
    #endif
    #if defined(__clang__)
        #define KFFT_CC_CLANG
    #endif
    #if defined(_CRAYC)
        #define KFFT_CC_CRAY
    #endif
#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)
    #define KFFT_CC_SUN
#elif defined(__xlc__) || defined(__xlC__)
    #define KFFT_CC_XLC
#elif defined(__HP_cc) || defined(__HP_aCC)
    #define KFFT_CC_HP
#elif defined(__DECC) || defined(__DECCXX)
    #define KFFT_CC_DEC
#elif (defined(__sgi) || defined(sgi)) &&                                                          \
    (defined(_COMPILER_VERSION) || defined(_SGI_COMPILER_VERSION))
    #define KFFT_CC_MIPS
#elif defined(__USLC__) && defined(__SCO_VERSION__)
    #define KFFT_CC_USLC
#elif defined(__WATCOMC__)
    #define KFFT_CC_WATCOM
#elif defined(__BORLANDC__)
    #define KFFT_CC_BORLAND
#elif defined(__INTEL_COMPILER)
    #define KFFT_CC_INTEL
#elif defined(__PGI)
    #define KFFT_CC_PGI
#elif defined(_CRAYC)
    #define KFFT_CC_CRAY
#endif

// === IMPORT / EXPORT ====================

#ifndef __has_attribute
    #define __has_attribute(x) (0)
#endif

#ifndef __kfft_export
    #if defined(KFFT_OS_WINDOWS) || defined(KFFT_OS_WINEMULATOR)
        #if defined(KFFT_CC_GNU) || __has_attribute(dllexport)
            #define __kfft_export __attribute__((dllexport))
        #elif defined(KFFT_CC_MSVC)
            #define __kfft_export __declspec(dllexport)
        #else
            #define __kfft_export
        #endif
    #elif defined(KFFT_CC_GNU) || __has_attribute(visibility)
        #define __kfft_export __attribute__((visibility("default")))
    #else
        #define __kfft_export
    #endif
#endif /* __kfft_export */

#ifndef __kfft_import
    #if defined(KFFT_OS_WINDOWS) || defined(KFFT_OS_WINEMULATOR)
        #if defined(KFFT_CC_GNU) || __has_attribute(dllimport)
            #define __kfft_import __attribute__((dllimport))
        #elif defined(KFFT_CC_MSVC)
            #define __kfft_import __declspec(dllimport)
        #else
            #define __kfft_import
        #endif
    #else
        #define __kfft_import
    #endif
#endif /* __kfft_import */

#if defined(KFFT_EXPORTS)
    #define KFFT_API __kfft_export
#elif defined(KFFT_IMPORTS)
    #define KFFT_API __kfft_import
#else
    #define KFFT_API
#endif

#define KFFT_UNUSED_VAR(X) (void)(X)

// Use: #pragma message WARN("My message")
#if defined(KFFT_CC_MSVC)
    #define FILE_LINE_LINK __FILE__ "(" STRINGISE(__LINE__) ") : "
    #define WARN(exp) (FILE_LINE_LINK "WARNING: " exp)
#else //__GNUC__ - may need other defines for different compilers
    #define WARN(exp) ("WARNING: " exp)
#endif
