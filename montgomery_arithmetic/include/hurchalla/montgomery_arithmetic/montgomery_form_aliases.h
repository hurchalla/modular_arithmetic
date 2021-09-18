// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_ALIASES_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_ALIASES_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>

namespace hurchalla {


// This file has an alias for a higher performance MontgomeryForm type called
// MontgomeryQuarter, which is available for use with smaller modulus sizes.
// But unless you wish to squeeze out every possible performance advantage, you
// will likely find it is more convenient to simply use MontgomeryForm<T>.
// This file also has a special purpose alias MontgomeryStandardMathWrapper
// which is described further below.
//
// USAGE:
// The suffix "Quarter" in the alias name MontgomeryQuarter indicates the size
// limit for the modulus that you are allowed to use to construct a
// MontgomeryQuarter object: you may use the smallest quarter of the range of
// all possible values of type T for the modulus.  More specifically, if we
// let R = (1 << ut_numeric_limits<T>::digits), MontgomeryQuarter<T> allows any
// modulus < R/4.  For example, MontgomeryQuarter<uint64_t> allows any modulus
// less than (1 << 64)/4.  It is undefined behavior to use a modulus that is not
// within the allowed range.  The modulus you use must also be odd, which is
// always required for montgomery arithmetic.
// In contrast, the default class MontgomeryForm<T> has no restriction on the
// its modulus size, though it too requires the modulus must be odd.
// You can expect that MontgomeryQuarter<T> will perform better (very often) or
// at worst the same as MontgomeryForm<T>, if both are given the same modulus.
//
// Note that this file also has an alias called MontgomeryStandardMathWrapper.
// This alias maps to a class that uses the MontgomeryForm interface but that
// internally performs all calculations with standard modular arithmetic rather
// than any montgomery arithmetic.  This can be useful as an aid for comparing
// performance between montgomery and non-montgomery modular arithmetic.
// Since MontgomeryStandardMathWrapper does not use montgomery arithmetic, its
// modulus is allowed to be either even or odd.
//
// PERFORMANCE DETAILS:
// The MontgomeryQuarter<T> alias can offer a significant performance
// improvement over MontgomeryForm<T>.  If you know either at compile-time or
// via run-time checks that your modulus will be small enough to allow you to
// use this alias, then you might expect it to provide perhaps a 5-20%
// performance gain over MontgomeryForm<T>.
// MontgomeryStandardMathWrapper<T> usually will perform worse than
// MontgomeryQuarter<T> and MontgomeryForm<T>, and often it performs much worse.
// However, on some modern systems with extremely fast dividers it is possible
// that it could outperform MontgomeryForm<T> and/or
// MontgomeryStandardMathWrapper<T>.
// With all performance details, you need to measure on your system to know what
// to truly expect.


template <typename, template <typename> class> class MontyAliasHelper;


// The MontgomeryQuarter alias
template <typename T>
using MontgomeryQuarter = MontgomeryForm<T,
                typename MontyAliasHelper<T, detail::MontyQuarterRange>::type>;

// The MontgomeryStandardMathWrapper alias provides the MontgomeryForm interface
// but uses no montgomery arithmetic.  All arithmetic is done with standard
// modular arithmetic instead.  This can be useful to compare performance of
// standard modular arithmetic with montgomery arithmetic - some systems with
// extremely fast divide operations could in theory perform better in some
// situations with standard modular arithmetic than with montgomery arithmetic.
// This wrapper lets you simply switch your type instead of rewriting code, when
// you want to compare performance.
// Ordinarily, montgomery arithmetic tends to be considerably faster than
// standard modular arithmetic whenever a large amount of modular multiplication
// is needed, and so this is probably unlikely to be an alias you would expect
// to use.  However, CPU architectures vary and evolve, and what is true today
// may not be true tomorrow - the only way to truly know what is fastest on a
// given system is to measure and compare performance.
template <typename T>
using MontgomeryStandardMathWrapper =
                MontgomeryForm<T, detail::MontyWrappedStandardMath<T>>;




// You should not use this class (it's intended for the alias implementations)
template <typename T, template <typename> class M>
class MontyAliasHelper final {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    static constexpr int bitsT = ut_numeric_limits<T>::digits;
    static constexpr int target_bits = HURCHALLA_TARGET_BIT_WIDTH;
public:
    using type =
        typename std::conditional<
            (bitsT <= target_bits - 2),
            detail::MontyQuarterRange<typename sized_uint<target_bits>::type>,
            M<T>
        >::type;
};


} // end namespace

#endif
