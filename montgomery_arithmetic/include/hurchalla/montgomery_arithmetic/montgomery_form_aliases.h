// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_ALIASES_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_ALIASES_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/MontyFullRangeMasked.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>

namespace hurchalla {


// This file has aliases for higher performance MontgomeryForm types called
// MontgomeryQuarter, MontgomeryHalf, and MontgomeryFull, each of which are
// usable within a specific range of modulus sizes.  But unless you wish to
// squeeze out every possible performance advantage, you will likely find it is
// more convenient to simply use MontgomeryForm<T>.  This file also has a
// special purpose alias MontgomeryStandardMathWrapper which is described
// further below.
//
// USAGE:
// The suffix "Quarter" in the alias name MontgomeryQuarter indicates the size
// limit for the modulus that you are allowed to use to construct a
// MontgomeryQuarter object: you may use the smallest quarter of the range of
// all possible (odd) values of the unsigned version of type T, for the modulus.
// More specifically, let U = extensible_make_unsigned<T>::type  and let U2 be
// a theoretical unlimited precision integer type.  We let the theoretical value
// R = (U2)1 << ut_numeric_limits<U>::digits.  MontgomeryQuarter<T> allows any
// odd positive modulus < R/4.  For example, MontgomeryQuarter<uint64_t> allows
// any odd modulus satisfying 0 < modulus < ((U2)1 << 64)/4.  And likewise,
// MontgomeryQuarter<int64_t> also allows any odd modulus satisfying
// 0 < modulus < ((U2)1 << 64)/4.  It is undefined behavior to use a modulus
// that is not within the allowed range.  The modulus you use must be odd, which
// is always required for montgomery arithmetic.
// In contrast, the default class MontgomeryForm<T> has no restriction on its
// modulus size, though it still requires that its modulus must be a positive
// odd number.
// You can expect that MontgomeryQuarter<T> will perform better (very often) or
// at worst the same as MontgomeryForm<T>, if both are given the same modulus.
//
// The suffix "Half" in the alias name MontgomeryHalf indicates the size limit
// for the modulus that you are allowed to use to construct a MontgomeryHalf
// object: you may use the smallest half of the range of all possible (odd)
// values of the unsigned version of type T, for the modulus.  More
// specifically, let U = extensible_make_unsigned<T>::type  and let U2 be a
// theoretical unlimited precision integer type.  We let the theoretical value
// R = (U2)1 << ut_numeric_limits<U>::digits.  MontgomeryHalf<T> allows any odd
// positive modulus < R/2.  For example, MontgomeryHalf<uint64_t> allows any odd
// modulus satisfying 0 < modulus < ((U2)1 << 64)/2.  MontgomeryHalf<int64_t>
// also allows any odd modulus satisfying 0 < modulus < ((U2)1 << 64)/2.  We can
// note that for any signed integer type T, this therefore permits all positive
// (and odd) values of T to be used for the modulus.  It is undefined behavior
// to use a modulus that is not within the allowed range.  The modulus you use
// must be odd, which is always required for montgomery arithmetic.
// In contrast, the default class MontgomeryForm<T> has no restriction on its
// modulus size, though it still requires that its modulus must be a positive
// odd number.
// You can expect MontgomeryHalf<T> to perform worse than MontgomeryQuarter<T>,
// but usually better than MontgomeryForm<T>, if you use the same modulus for
// all these Montgomery types.
// For a type T that is the same size as the CPU integer registers (e.g.
// uint64_t on a 64 bit computer) or a type T that is smaller than the register
// size, you can expect that MontgomeryHalf<T> will perform better (very often)
// or at worst the same as MontgomeryForm<T>, if both are given the same
// modulus.  It is possible that plain add() and subtract() may perform slightly
// worse, but if so, this would ordinarily be overcome by the improved
// performance of the multiply, square, fused-multiply/square-add/sub functions.
// However, for a type T that is larger than the CPU integer register size, it
// is not unusual for MontgomeryHalf<T> to perform worse overall than
// MontgomeryForm<T>.
// In general, if your modulus is small enough to allow MontgomeryQuarter<T>,
// you should prefer to use it instead of MontgomeryHalf<T>.
//
// The suffix "Full" in the alias name MontgomeryFull indicates that any
// positive odd-valued modulus is permissable to use to construct a
// MontgomeryFull object: you may use the full range of all positive odd values
// of type T for the modulus.  MontgomeryFull utilizes the standard, normal
// Montgomery algorithms without any interesting or unusal optimizations to the
// algorithms.  Usually MontgomeryForm<T> maps to the same underlying class as
// MontgomeryFull<T>, and so they often perform the same.  However
// MontgomeryForm<T> can map to more efficient classes in some cases.  For this
// reason, you should usually prefer to use MontgomeryForm<T>.
//
// Note that this file also has an alias called MontgomeryStandardMathWrapper.
// This alias maps to a class that uses the MontgomeryForm interface but that
// internally performs all calculations with standard modular arithmetic rather
// than any montgomery arithmetic.  This can be useful as a convenient aid to
// compare performance between montgomery and non-montgomery modular arithmetic,
// since the fact that these aliases both provide the exact same API allows you
// to use them interchangably via a function/struct/class template parameter.
// Since MontgomeryStandardMathWrapper does not use montgomery arithmetic, its
// modulus is allowed to be either even or odd.  Its modulus must be positive.
//
// PERFORMANCE DETAILS:
// The MontgomeryQuarter<T> and MontgomeryHalf<T> can offer a notable
// performance improvement over MontgomeryForm<T>.  If you know either at
// compile-time or via run-time checks that your modulus will be small enough to
// allow you to use one of these aliases, then you might roughly expect
// performance gains perhaps in the range of 5-20% over MontgomeryForm<T>.
// MontgomeryStandardMathWrapper<T> usually will perform worse than all the
// other classes and aliases mentioned here, and often it performs much worse.
// However, on some modern systems with extremely fast dividers it is possible
// that it could outperform both MontgomeryForm<T> and the normal aliases.
//
// With all performance details, you need to measure on your system to know what
// to expect in reality.


template <typename, template <typename> class> class MontyAliasHelper;


// For the template parameter InlineAllFunctions, you should usually accept the
// default rather than specify an argument.  However if you wish to reduce
// compilation times or potentially reduce your executable's size, you can set
// it to false, which may help.

template <typename T, bool InlineAllFunctions = true>
using MontgomeryQuarter = MontgomeryForm<T, InlineAllFunctions,
         detail::MontyQuarterRange<typename extensible_make_unsigned<T>::type>>;

template <typename T, bool InlineAllFunctions = true>
using MontgomeryHalf = MontgomeryForm<T, InlineAllFunctions,
              typename MontyAliasHelper<T, detail::MontyHalfRange>::type>;

template <typename T, bool InlineAllFunctions = true>
using MontgomeryFull = MontgomeryForm<T, InlineAllFunctions,
              typename MontyAliasHelper<T, detail::MontyFullRange>::type>;


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
// may not be true tomorrow - you will need to measure and compare performance
// on your system to know for certain.
template <typename T, bool InlineAllFunctions = true>
using MontgomeryStandardMathWrapper = MontgomeryForm<T, InlineAllFunctions,
  detail::MontyWrappedStandardMath<typename extensible_make_unsigned<T>::type>>;




// You should not use this class (it's intended for the alias implementations)
template <typename T, template <typename> class M>
class MontyAliasHelper final {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    using U = typename extensible_make_unsigned<T>::type;
    static constexpr int bitsT = ut_numeric_limits<T>::digits;
    static constexpr int target_bits = HURCHALLA_TARGET_BIT_WIDTH;
public:
    using type =
        typename std::conditional<
            (bitsT <= target_bits - 2),
            detail::MontyQuarterRange<typename sized_uint<target_bits>::type>,
            M<U>
        >::type;
};
// In perf testing on a 64 bit system, I found that MontyQuarterRange<uint64_t>
// performed the same or slightly better than MontyHalfRange<uint32_t> or
// MontyFullRange<uint32_t>, when given the same (uint32_t) modulus.  This is
// reflected above by MontyAliasHelper mapping to MontyQuarterRange if it gets a
// small enough type T.


// experimental alias - you should not use this
template <typename T, bool InlineAllFunctions = true>
using MontgomeryMasked = MontgomeryForm<T, InlineAllFunctions,
              typename MontyAliasHelper<T, detail::MontyFullRangeMasked>::type>;



} // end namespace

#endif
