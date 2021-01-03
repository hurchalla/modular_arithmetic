// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_ALIASES_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_ALIASES_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySixthRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>

namespace hurchalla { namespace montgomery_arithmetic {


// These are aliases for high performance MontgomeryForm types that are tailored
// to the size of the modulus you plan to use.  But unless you wish to squeeze
// out every possible performance advantage, you will likely find it is more
// convenient to simply use MontgomeryForm<T>, instead of these aliases.
//
// USAGE:
// The suffix (Full, Half, Quarter, Sixth) in each alias indicates the size 
// limit allowed for the modulus you can use to construct an object of each
// alias type.  Full allows you to use the full range of type T for the modulus.
// Half allows you to use the smaller half of the range of possible T values, or
// more specifically, if we let R = 2^(ut_numeric_limits<T>::digits), Half
// allows any modulus < R/2.  Quarter allows any modulus < R/4, and Sixth allows
// any modulus < R/6.  To give an example, MontgomeryQuarter<uint64_t> allows
// any modulus less than (2^64)/4.  It is undefined behavior to use a modulus
// that is not within the allowed range for your alias type.  As with all
// montgomery arithmetic, the modulus you use must also be odd.
// Relatively speaking, an alias type with a smaller modulus limit will often
// perform better than an alias with a larger limit; you can expect that the
// smaller alias will never perform worse than the larger alias, if given the
// same modulus.
// (This file also provides an alias called MontgomeryStandardMathWrapper.  This
// alias maps to a class that uses the MontgomeryForm interface but that
// internally performs all calculations with standard modular arithmetic rather
// than any montgomery arithmetic.  This can be useful as an aid for comparing
// performance between montgomery and non-montgomery modular arithmetic.)
//
// PERFORMANCE DETAILS:
// Note that these aliases do not always improve upon the performance of
// MontgomeryForm<T> - the performance gain (or lack of gain) you can expect
// from these aliases depends on the particular alias:
// (1) The MontgomeryFull<T> alias maps to the same class as MontgomeryForm<T>.
// Thus it offers no improvement on performance over MontgomeryForm<T>.
// (2) The MontgomeryQuarter<T> and MontgomerySixth<T> aliases offer significant
// performance improvements over MontgomeryForm<T>.  If you know either at
// compile-time or via run-time checking that your modulus will be small enough
// to allow you to use one of these aliases, then you might expect these aliases
// to provide perhaps a 5-20% performance gain over MontgomeryForm<T>.
// (3) The alias MontgomerySixth<T> improves upon MontgomeryQuarter<T> only for
// the famul() member function; otherwise it's the same as MontgomeryQuarter<T>.
// (4) The alias MontgomeryHalf<T> improves upon MontgomeryFull<T> only for the
// famul() member function; it's otherwise equivalent to MontgomeryFull<T>.


template <typename, template <typename> class> class MontyAliasHelper;


// The Montgomery Form Aliases
// ---------------------------

template <typename T>
using MontgomeryFull = MontgomeryForm<T,
                typename MontyAliasHelper<T, detail::MontyFullRange>::type>;

template <typename T>
using MontgomeryHalf = MontgomeryForm<T,
                typename MontyAliasHelper<T, detail::MontyHalfRange>::type>;

template <typename T>
using MontgomeryQuarter = MontgomeryForm<T,
                typename MontyAliasHelper<T, detail::MontyQuarterRange>::type>;

template <typename T>
using MontgomerySixth = MontgomeryForm<T,
                typename MontyAliasHelper<T, detail::MontySixthRange>::type>;

// The MontgomeryStandardMathWrapper alias provides the MontgomeryForm interface
// but uses no montgomery arithmetic.  All arithmetic is done with standard
// modular arithmetic instead.  This can be useful to compare performance of
// standard modular arithmetic with montgomery arithmetic - some systems with
// extremely fast divide operations may perform better in some situations with
// standard modular arithmetic than with montgomery arithmetic.  This wrapper
// lets you simply switch your type instead of rewriting code, when you want to
// compare performance.
// Ordinarily, this is probably not an alias you should expect to use.
template <typename T>
using MontgomeryStandardMathWrapper =
                MontgomeryForm<T, detail::MontyWrappedStandardMath<T>>;




// You should not use this class (it's intended for the alias implementations)
template <typename T, template <typename> class M>
class MontyAliasHelper {
    static_assert(util::ut_numeric_limits<T>::is_integer, "");
    static_assert(!util::ut_numeric_limits<T>::is_signed, "");
    static_assert(util::ut_numeric_limits<T>::is_modulo, "");
    static constexpr int bitsT = util::ut_numeric_limits<T>::digits;
    static constexpr int target_bits = HURCHALLA_TARGET_BIT_WIDTH;
public:
    using type =
        typename std::conditional<
            std::is_same<M<T>, detail::MontySixthRange<T>>::value,
            M<T>,
        typename std::conditional<
            util::sized_uint<bitsT*2>::is_valid && (bitsT*2 < target_bits),
            detail::MontySixthRange<typename util::sized_uint<bitsT*2>::type>,
        typename std::conditional<
            (bitsT <= target_bits - 3) &&
                        (std::is_same<M<T>, detail::MontyFullRange<T>>::value ||
                         std::is_same<M<T>, detail::MontyHalfRange<T>>::value),
            detail::MontySixthRange<typename hurchalla::util::sized_uint<target_bits>::type>,
            M<T>
        >::type>::type>::type;
};


}} // end namespace

#endif
