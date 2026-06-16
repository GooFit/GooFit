#pragma once

#include <cfenv>

#if defined(__APPLE__) && defined(__MACH__)

// macOS does not ship feenableexcept/fedisableexcept, so provide them.

#if defined(__arm64__) || defined(__aarch64__)

// On Apple Silicon the trap-enable bits live in FPCR and are the FE_* flags
// shifted left by 8 (FE_INVALID 0x01 -> __fpcr_trap_invalid 0x100, etc.).
// Note: current Apple M-series cores ignore these bits and never raise SIGFPE,
// so this is effectively a no-op at runtime, but it keeps the API and compiles.

inline auto feenableexcept(unsigned int excepts) -> int {
    fenv_t fenv;
    if(fegetenv(&fenv)) {
        return -1;
    }
    unsigned int old_excepts = (fenv.__fpcr >> 8) & FE_ALL_EXCEPT;
    fenv.__fpcr |= static_cast<unsigned long long>(excepts & FE_ALL_EXCEPT) << 8;
    return fesetenv(&fenv) ? -1 : static_cast<int>(old_excepts);
}

inline auto fedisableexcept(unsigned int excepts) -> int {
    fenv_t fenv;
    if(fegetenv(&fenv)) {
        return -1;
    }
    unsigned int old_excepts = (fenv.__fpcr >> 8) & FE_ALL_EXCEPT;
    fenv.__fpcr &= ~(static_cast<unsigned long long>(excepts & FE_ALL_EXCEPT) << 8);
    return fesetenv(&fenv) ? -1 : static_cast<int>(old_excepts);
}

#else

// Public domain polyfill for feenableexcept on x86 OS X
// http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c
// Fenv wrapper from:
// https://github.com/ArduPilot/ardupilot/blob/master/libraries/AP_Common/missing/fenv.h

inline auto feenableexcept(unsigned int excepts) -> int {
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;

    // previous masks
    unsigned int old_excepts;

    if(fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr &= ~(new_excepts << 7);

    return fesetenv(&fenv) ? -1 : old_excepts;
}

inline auto fedisableexcept(unsigned int excepts) -> int {
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    // all previous masks
    unsigned int old_excepts;

    if(fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // mask
    fenv.__control |= new_excepts;
    fenv.__mxcsr |= new_excepts << 7;

    return fesetenv(&fenv) ? -1 : old_excepts;
}

#endif

#endif
