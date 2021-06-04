#pragma once

#include <cfenv>

#if defined(__APPLE__) && defined(__MACH__)

// Public domain polyfill for feenableexcept on OS X
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
