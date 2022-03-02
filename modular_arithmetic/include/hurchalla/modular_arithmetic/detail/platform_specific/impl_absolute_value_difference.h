// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED


#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// note: uses a static member function to disallow ADL.
struct default_impl_absdiff_unsigned {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
      // result = (a < b) ? static_cast<T>(b - a) : static_cast<T>(a - b)
    T result;
    HURCHALLA_CSELECT(result, a < b, static_cast<T>(b-a), static_cast<T>(a-b));

    // POSTCONDITION:
    // This function returns absolute_value(a-b).
    HPBC_POSTCONDITION(result<=a || result<=b);
    return result;
  }
};


// primary template
template <typename T>
struct impl_absolute_value_difference_unsigned {
  HURCHALLA_FORCE_INLINE static T call(T a, T b)
  {
    return default_impl_absdiff_unsigned::call(a, b);
  }
};


// MSVC doesn't support inline asm so we skip it.
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_ABSDIFF)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

template <>
struct impl_absolute_value_difference_unsigned<std::uint32_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint32_t call(std::uint32_t a, std::uint32_t b)
  {
    using std::uint32_t;
    uint32_t diff = static_cast<uint32_t>(b - a);
    uint32_t tmp = a;  // we prefer not to overwrite an input (a)
    __asm__ ("subl %[b], %[tmp] \n\t"       /* tmp = a - b */
             "cmovbl %[diff], %[tmp] \n\t"  /* tmp = (a < b) ? diff : tmp */
             : [tmp]"+&r"(tmp)
             : [b]"r"(b), [diff]"r"(diff)
             : "cc");
    uint32_t result = tmp;

    HPBC_POSTCONDITION2(result<=a || result<=b);
    HPBC_POSTCONDITION2(result == default_impl_absdiff_unsigned::call(a, b));
    return result;
  }
};

template <>
struct impl_absolute_value_difference_unsigned<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t a, std::uint64_t b)
  {
    using std::uint64_t;
    uint64_t diff = static_cast<uint64_t>(b - a);
    uint64_t tmp = a;  // we prefer not to overwrite an input (a)
    __asm__ ("subq %[b], %[tmp] \n\t"       /* tmp = a - b */
             "cmovbq %[diff], %[tmp] \n\t"  /* tmp = (a < b) ? diff : tmp */
             : [tmp]"+&r"(tmp)
             : [b]"r"(b), [diff]"r"(diff)
             : "cc");
    uint64_t result = tmp;

    HPBC_POSTCONDITION2(result<=a || result<=b);
    HPBC_POSTCONDITION2(result == default_impl_absdiff_unsigned::call(a, b));
    return result;
  }
};

#endif


// version for unsigned T
template <typename T, bool = ut_numeric_limits<T>::is_signed>
struct impl_absolute_value_difference {
  HURCHALLA_FORCE_INLINE static T call(T a, T b)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    return impl_absolute_value_difference_unsigned<T>::call(a, b);
  }
};

// version for signed T
template <typename T>
struct impl_absolute_value_difference<T, true> {
  HURCHALLA_FORCE_INLINE static T call(T a, T b)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::is_signed, "");
    static_assert(static_cast<T>(-1) == ~(static_cast<T>(0)),
                                  "T must use two's complement representation");
    using U = typename extensible_make_unsigned<T>::type;
    static_assert(static_cast<T>(static_cast<U>(static_cast<T>(-1))) ==
                  static_cast<T>(-1), "Casting a signed T value to unsigned and"
                               " back again must result in the original value");
    HPBC_PRECONDITION2(a >= 0);
    HPBC_PRECONDITION2(b >= 0);
#if defined(HURCHALLA_AVOID_CSELECT)
    static_assert((static_cast<T>(-1) >> 1) == static_cast<T>(-1),
                          "Arithmetic right shift is required but unavailable");
    T diff = static_cast<T>(a - b);
    // if diff is negative, create a bit mask of all 1s.  Otherwise all 0s.
    U mask = static_cast<U>(diff >> ut_numeric_limits<T>::digits);
    // We now calculate the absolute value of diff.  This method comes from
    // https://graphics.stanford.edu/~seander/bithacks.html#IntegerAbs
    // The formula they give of  abs(v) = (v ^ mask) - mask
    // works because if v >= 0, then mask is all 0s, and the result is v, as
    // desired.  If v < 0, then mask is all 1s, and so the xor inverts all the
    // bits of v, followed by a subtraction of -1, which is the same as an
    // addition of 1;  these operations exactly produce the two's complement
    // (of the negative integer v), via the very well-known method of
    // "The two's complement is calculated by inverting the bits and adding one"
    // [quoted from https://en.wikipedia.org/wiki/Two%27s_complement]
    U tmp = static_cast<U>(static_cast<U>(diff) ^ mask);
    U result = static_cast<U>(tmp - mask);
    HPBC_ASSERT2(result == impl_absolute_value_difference_unsigned<U>::call(
                                         static_cast<U>(a), static_cast<U>(b)));
#else
    U result = impl_absolute_value_difference_unsigned<U>::call(
                                          static_cast<U>(a), static_cast<U>(b));
#endif
    return static_cast<T>(result);
  }
};


}}  // end namespace

#endif
