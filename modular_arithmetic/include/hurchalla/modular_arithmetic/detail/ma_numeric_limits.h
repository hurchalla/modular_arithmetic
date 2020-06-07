
#ifndef HURCHALLA_MODULAR_ARITHMETIC_NUMERIC_LIMITS_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_NUMERIC_LIMITS_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include <limits>
#include <climits>

namespace hurchalla { namespace modular_arithmetic {


// Rationale for ma_numeric_limits:
// Specializations for __int128_t and __uint128_t  are the reason we use
// ma_numeric_limits throughout the modular arithmetic library instead of
// std::numeric_limits.  Std::numeric_limits is too problematic on its own when
// it comes to the 128 bit int compiler extensions.  More specifically,
// gcc/clang/icc provide __int128_t and __uint128_t, but these compilers do not
// always specialize std::numeric_limits for those types in all situations where
// the types are valid, available, and usable.  For example, clang does not
// provide specializations if it is passed the flag -std=c++11 (or any other
// version xx in -std=c++xx), despite the fact that __int128_t and __uint128_t
// are still available.  Without being passed the flag -std=c++xx, clang
// usually does provide the std::numeric_limits specializations.  For more
// information see
// https://quuxplusone.github.io/blog/2019/02/28/is-int128-integral/
//
// By providing our own traits struct ma_numeric_limits, we can ensure that
// ma_numeric_limits is always specialized for __int128_t and __uint128_t when
// those types are valid and usable.



// Primary template:
template <typename T>
struct ma_numeric_limits : std::numeric_limits<T> {};



#if (HURCHALLA_COMPILER_HAS_UINT128_T())
// Explicit specializations for __int128_t and __uint128_t.
// Note: we only define members that the modular arithmetic library uses (rather
// than defining all the members that would be in std::numeric_limits).
template<>
struct ma_numeric_limits<__int128_t> {
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = true;
    static constexpr bool is_modulo = false;
    static constexpr int digits = 127;
    static constexpr __int128_t max() noexcept
    {
        // Return (2^127)-1 without overflow
        return ((static_cast<__int128_t>(1) << 126) - 1) +
                                            (static_cast<__int128_t>(1) << 126);
        // The return value relies on the strong evidence that  (2^127)-1  must
        // be the max value of __int128_t.
        // For two's complement, ones' complement, signed magnitude, or offset
        // binary/excess-k, this would be the correct value.  It would be wrong
        // for base -2 (negative two) representation, but it seems exceedingly
        // unlikely that a compiler would use that representation or anything
        // different from two's complement/ones' complement/signed magnitude.
        // When I look at the generated asm from the latest versions (as of
        // 6/7/20) of gcc clang and icc, adding or subtracting two __int128_t
        // values, or casting a single __int128_t to int64_t, all appear to rely
        // on two's complement.
        //
        // Nevertheless we can't be absolutely certain this value is correct for
        // all compilers and compiler versions.  Unfortunately I don't see any
        // good reliable way to get the max __int128_t value from the compiler,
        // considering that we already know clang/gcc won't specialize
        // std::numeric_limits for __int128_t  under a number of circumstances.
        //
        // Note: Perhaps one way to get the correct value when the compiler does
        // not provide it, would be to compile and run a small program with gcc
        // extensions enabled (which should get the compiler to provide a
        // specialization std::numeric_limits<__int128_t>) and print the max
        // value in a dedicated header file.  Then during any build of modular
        // arithmetic, that utility program would build and run first, and then
        // in a subsequent build, this ma_numeric_limits header file would
        // #include the header that the utility program wrote.  This adds an
        // awkward extra step to building, and extra complication overall.  I
        // believe the very strong guess used above is preferable to the added
        // complexity.
    }
};
template<>
struct ma_numeric_limits<__uint128_t> {
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = false;
    static constexpr bool is_integer = true;
    static constexpr bool is_modulo = true;
    static constexpr int digits = 128;
    static constexpr __uint128_t max() noexcept
    {
        // rely upon modulo (2^128) wrap-around on subtraction below
        return static_cast<__uint128_t>(0) - static_cast<__uint128_t>(1);
    }
};
#endif


}}  // end namespace

#endif
