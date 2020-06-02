
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
    static constexpr __uint128_t max() noexcept
    {
        constexpr __int128_t one = static_cast<__int128_t>(1);
        constexpr __int128_t twoPow126 = one << 126;
        // Return (2^127)-1 without overflow
        return (twoPow126 - 1) + twoPow126;
        // Note that this relies on a guess that  (2^127)-1  is the max value
        // of __int128_t (and that it's representable in __int128_t).  For
        // two's complement this would be true.  But this isn't a given.
        // Unfortunately I don't see any reliable way to get the max value from
        // the compiler - considering that we already know clang/gcc won't
        // specialize std::numeric_limits for __int128_t  under a number of
        // circumstances.  *Maybe* stdint.h would have this information, but at
        // the moment I'm assuming even if it did, it would be as unreliable as
        // std::numeric_limits.
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
