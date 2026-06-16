#pragma once

#include <string>

// Self-contained replacement for the (x86-only, unmaintained) FeatureDetector
// submodule. Warns when the host CPU advertises a SIMD feature that the
// compiler was not told to use (i.e. the build is missing -march=native or an
// equivalent flag).
//
// Host detection uses the GCC/Clang __builtin_cpu_supports() builtin, so there
// is no external dependency. On non-x86 targets (e.g. aarch64, where NEON is
// part of the baseline) and on compilers without the builtin (MSVC) the check
// is a no-op.

namespace GooFit {

namespace detail {

inline auto feature_warning(const std::string &name) -> std::string {
    return "Warning: " + name + " supported by the host but not enabled during compilation\n";
}

// Each warning fires only when the host advertises the feature (runtime) and
// the matching compile-time macro is absent. The macro half is folded in at
// preprocess time so the host check is skipped entirely on a fully optimized build.
#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#ifdef __SSE__
#define GOOFIT_COMP_SSE true
#else
#define GOOFIT_COMP_SSE false
#endif
#ifdef __SSE2__
#define GOOFIT_COMP_SSE2 true
#else
#define GOOFIT_COMP_SSE2 false
#endif
#ifdef __SSE3__
#define GOOFIT_COMP_SSE3 true
#else
#define GOOFIT_COMP_SSE3 false
#endif
#ifdef __SSE4_1__
#define GOOFIT_COMP_SSE41 true
#else
#define GOOFIT_COMP_SSE41 false
#endif
#ifdef __SSE4_2__
#define GOOFIT_COMP_SSE42 true
#else
#define GOOFIT_COMP_SSE42 false
#endif
#ifdef __AVX__
#define GOOFIT_COMP_AVX true
#else
#define GOOFIT_COMP_AVX false
#endif
#ifdef __AVX2__
#define GOOFIT_COMP_AVX2 true
#else
#define GOOFIT_COMP_AVX2 false
#endif
#ifdef __AVX512F__
#define GOOFIT_COMP_AVX512F true
#else
#define GOOFIT_COMP_AVX512F false
#endif
#endif

} // namespace detail

/// Build a string of warnings for SIMD features the host supports but the
/// compiler did not enable. Empty when fully optimized or on an unsupported target.
inline auto cpu_feature_warnings() -> std::string {
    std::string output;

#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
    if(!GOOFIT_COMP_SSE && __builtin_cpu_supports("sse"))
        output += detail::feature_warning("SSE");
    if(!GOOFIT_COMP_SSE2 && __builtin_cpu_supports("sse2"))
        output += detail::feature_warning("SSE2");
    if(!GOOFIT_COMP_SSE3 && __builtin_cpu_supports("sse3"))
        output += detail::feature_warning("SSE3");
    if(!GOOFIT_COMP_SSE41 && __builtin_cpu_supports("sse4.1"))
        output += detail::feature_warning("SSE4.1");
    if(!GOOFIT_COMP_SSE42 && __builtin_cpu_supports("sse4.2"))
        output += detail::feature_warning("SSE4.2");
    if(!GOOFIT_COMP_AVX && __builtin_cpu_supports("avx"))
        output += detail::feature_warning("AVX");
    if(!GOOFIT_COMP_AVX2 && __builtin_cpu_supports("avx2"))
        output += detail::feature_warning("AVX2");
    if(!GOOFIT_COMP_AVX512F && __builtin_cpu_supports("avx512f"))
        output += detail::feature_warning("AVX512F");
#endif

    return output;
}

} // namespace GooFit

#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#undef GOOFIT_COMP_SSE
#undef GOOFIT_COMP_SSE2
#undef GOOFIT_COMP_SSE3
#undef GOOFIT_COMP_SSE41
#undef GOOFIT_COMP_SSE42
#undef GOOFIT_COMP_AVX
#undef GOOFIT_COMP_AVX2
#undef GOOFIT_COMP_AVX512F
#endif
