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


// These are aliases for high performance MontgomeryForm types that are tailored
// to the size of the modulus you plan to use.  But unless you wish to squeeze
// out every possible performance advantage, you will likely find it is more
// convenient to simply use MontgomeryForm<T>, instead of these aliases.
//
// USAGE:
// The suffix (Full, Quarter) in each alias indicates the size limit allowed for
// the modulus that you can use to construct an object of each alias type.  Full
// allows you to use the full range of type T for the modulus.  Quarter allows
// you to use the smallest quarter of the range of all possible T values, or
// more specifically, if we let R = 2^(ut_numeric_limits<T>::digits), Quarter
// allows any modulus < R/4.  For example, MontgomeryQuarter<uint64_t> allows
// any modulus less than (2^64)/4.  It is undefined behavior to use a modulus
// that is not within the allowed range for your alias type.  As with all
// montgomery arithmetic, the modulus you use must also be odd.
// You can expect that an alias type with a smaller modulus limit will perform
// better (very often) or at worst the same as an alias with a larger limit, if
// both are given the same modulus.
// Note that this file also has an alias called MontgomeryStandardMathWrapper.
// This alias maps to a class that uses the MontgomeryForm interface but that
// internally performs all calculations with standard modular arithmetic rather
// than any montgomery arithmetic.  This can be useful as an aid for comparing
// performance between montgomery and non-montgomery modular arithmetic.
//
// PERFORMANCE DETAILS:
// Note that these aliases do not all improve upon the performance of
// MontgomeryForm<T> - the performance gain (or lack of gain) that you can
// expect from these aliases depends upon the particular alias:
// (1) The MontgomeryFull<T> alias maps to the same class as MontgomeryForm<T>.
// It therefore offers no improvement on performance over MontgomeryForm<T>.
// (2) The MontgomeryQuarter<T> alias can offer a significant performance
// improvement over MontgomeryForm<T>.  If you know either at compile-time or
// via run-time checks that your modulus will be small enough to allow you to
// use this alias, then you might expect it to provide perhaps a 5-20%
// performance gain over MontgomeryForm<T>.


template <typename, template <typename> class> class MontyAliasHelper;


// The Montgomery Form Aliases
// ---------------------------

template <typename T>
using MontgomeryFull = MontgomeryForm<T,
                typename MontyAliasHelper<T, detail::MontyFullRange>::type>;

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
class MontyAliasHelper {
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
