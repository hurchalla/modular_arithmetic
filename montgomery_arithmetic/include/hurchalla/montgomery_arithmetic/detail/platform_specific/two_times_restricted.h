// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_TWO_TIMES_RESTRICTED_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_TWO_TIMES_RESTRICTED_H_INCLUDED


#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace detail {


// This function:  two_times_restricted::call(a, modulus)  returns 
// (a + a) % modulus, performed as if 'a' and 'modulus' are infinite precision.
// The function has two critical restrictions that must be satisfied -
//   1) it requires modulus < R/2, where conceptually, the unlimited precision
// integer R is  R = 1 << ut_numeric_limits<make_unsigned<T>::type>digits.
//   2) the input 'a' must be prereduced, i.e. we must have  a < modulus.


// Fyi: the purpose of having structs with static member functions is to
// disallow ADL and to make specializations simple and easy.


#if !defined(__clang__)
// Note: gcc(x64 and arm32/64 and risc-V) and MSVC(x64 and arm64) tend to
// compile better machine code from the #if version.  clang(on x64 and arm32/64)
// seems to produce better code from the #else version, despite the extra
// complexity.

struct default_two_times_restricted_unsigned {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(a < modulus);  // i.e. the input must be prereduced

    // note: because modulus < R/2, and a < modulus, we know a+a won't overflow.
    T sum = static_cast<T>(a + a);
    T result = static_cast<T>(sum - modulus);
      // result = (sum < modulus) ? sum : result
    result = ::hurchalla::conditional_select(sum < modulus, sum, result);

    HPBC_CLOCKWORK_POSTCONDITION2(static_cast<T>(0) <= result && result < modulus);
    return result;
  }
};
#else

struct default_two_times_restricted_unsigned {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(a < modulus);  // the input must be prereduced

    T sum = static_cast<T>(a + a);
    T tmp = static_cast<T>(a - modulus);
    T result = static_cast<T>(a + tmp);
      // result = (result >= a) ? sum : result
    result = ::hurchalla::conditional_select(result >= a, sum, result);

    HPBC_CLOCKWORK_POSTCONDITION2(static_cast<T>(0) <= result && result < modulus);
    return result;
  }
};
#endif




// primary template
template <typename T>
struct two_times_restricted_unsigned {
  HURCHALLA_FORCE_INLINE static T call(T a, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    return default_two_times_restricted_unsigned::call(a, modulus);
  }
};



// These inline asm functions implement optimizations of the default function

// MSVC doesn't support inline asm, so we skip it.
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_TWOTIMES)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

# if (HURCHALLA_COMPILER_HAS_UINT128_T())
template <>
struct two_times_restricted_unsigned<__uint128_t> {
  HURCHALLA_FORCE_INLINE static
  __uint128_t call(__uint128_t a, __uint128_t modulus)
  {
    using std::uint64_t;
    using T = __uint128_t;
    constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(a < modulus);  // __uint128_t guarantees a>=0.

    __uint128_t sum = static_cast<__uint128_t>(a + a);
    uint64_t sumlo = static_cast<uint64_t>(sum);
    uint64_t sumhi = static_cast<uint64_t>(sum >> 64);
    uint64_t tmplo = sumlo;
    uint64_t tmphi = sumhi;
    uint64_t mlo = static_cast<uint64_t>(modulus);
    uint64_t mhi = static_cast<uint64_t>(modulus >> 64);
    __asm__ ("subq %[mlo], %[sumlo] \n\t"       /* tmp2 = sum - m */
             "sbbq %[mhi], %[sumhi] \n\t"
             "cmovbq %[tmplo], %[sumlo] \n\t"   /* sum = (sum<m) ? tmp : tmp2 */
             "cmovbq %[tmphi], %[sumhi] \n\t"
             : [sumlo]"+&r"(sumlo), [sumhi]"+&r"(sumhi)
#  if defined(__clang__)       /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [mlo]"r"(mlo), [mhi]"r"(mhi), [tmplo]"r"(tmplo), [tmphi]"r"(tmphi)
#  else
             : [mlo]"rm"(mlo), [mhi]"rm"(mhi), [tmplo]"rm"(tmplo), [tmphi]"rm"(tmphi)
#  endif
             : "cc");
    __uint128_t result = (static_cast<__uint128_t>(sumhi) << 64) | sumlo;

    HPBC_CLOCKWORK_POSTCONDITION2(result < modulus);  // __uint128_t guarantees result>=0.
    HPBC_CLOCKWORK_POSTCONDITION2(result ==
                       default_two_times_restricted_unsigned::call(a, modulus));
    return result;
  }
};
# endif

template <>
struct two_times_restricted_unsigned<std::uint64_t> {
  HURCHALLA_FORCE_INLINE static
  std::uint64_t call(std::uint64_t a, std::uint64_t modulus)
  {
    using std::uint64_t;
    using T = uint64_t;
    constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(a < modulus);  // uint64_t guarantees a>=0.

    uint64_t sum = static_cast<uint64_t>(a + a);
    uint64_t tmp = sum;
    __asm__ ("subq %[m], %[sum] \n\t"         /* tmp2 = sum - m */
             "cmovbq %[tmp], %[sum] \n\t"     /* sum = (sum<m) ? tmp : tmp2 */
             : [sum]"+&r"(sum)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [m]"r"(modulus), [tmp]"r"(tmp)
# else
             : [m]"rm"(modulus), [tmp]"rm"(tmp)
# endif
             : "cc");
    uint64_t result = sum;

    HPBC_CLOCKWORK_POSTCONDITION2(result < modulus);  // uint64_t guarantees result>=0.
    HPBC_CLOCKWORK_POSTCONDITION2(result ==
                       default_two_times_restricted_unsigned::call(a, modulus));
    return result;
  }
};

template <>
struct two_times_restricted_unsigned<std::uint32_t> {
  HURCHALLA_FORCE_INLINE static
  std::uint32_t call(std::uint32_t a, std::uint32_t modulus)
  {
    using std::uint32_t;
    using T = uint32_t;
    constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(a < modulus);  // uint32_t guarantees a>=0.

    uint32_t sum = static_cast<uint32_t>(a + a);
    uint32_t tmp = sum;
    __asm__ ("subl %[m], %[sum] \n\t"         /* tmp2 = sum - m */
             "cmovbl %[tmp], %[sum] \n\t"     /* sum = (sum<m) ? tmp : tmp2 */
             : [sum]"+&r"(sum)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [m]"r"(modulus), [tmp]"r"(tmp)
# else
             : [m]"rm"(modulus), [tmp]"rm"(tmp)
# endif
             : "cc");
    uint32_t result = sum;

    HPBC_CLOCKWORK_POSTCONDITION2(result < modulus);  // uint32_t guarantees result>=0.
    HPBC_CLOCKWORK_POSTCONDITION2(result ==
                       default_two_times_restricted_unsigned::call(a, modulus));
    return result;
  }
};

// end of inline asm functions for x86_64
#endif





// ARM64
// MSVC doesn't support inline asm so we skip it.
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_TWOTIMES)) && \
    defined(HURCHALLA_TARGET_ISA_ARM_64) && !defined(_MSC_VER)

# if (HURCHALLA_COMPILER_HAS_UINT128_T())
template <>
struct two_times_restricted_unsigned<__uint128_t> {
  HURCHALLA_FORCE_INLINE
  static __uint128_t call(__uint128_t a, __uint128_t modulus)
  {
    using std::uint64_t;
    using T = __uint128_t;
    constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(a < modulus);  // __uint128_t guarantees a>=0.

    __uint128_t sum = static_cast<__uint128_t>(a + a);
    uint64_t sumlo = static_cast<uint64_t>(sum);
    uint64_t sumhi = static_cast<uint64_t>(sum >> 64);
    uint64_t mlo = static_cast<uint64_t>(modulus);
    uint64_t mhi = static_cast<uint64_t>(modulus >> 64);
    uint64_t reslo;
    uint64_t reshi;
    __asm__ ("subs %[reslo], %[sumlo], %[mlo] \n\t"        /* res = sum - m */
             "sbcs %[reshi], %[sumhi], %[mhi] \n\t"
             "csel %[reslo], %[sumlo], %[reslo], lo \n\t"  /* res = (sum<m) ? sum : res */
             "csel %[reshi], %[sumhi], %[reshi], lo \n\t"
             : [reslo]"=&r"(reslo), [reshi]"=&r"(reshi)
             : [mlo]"r"(mlo), [mhi]"r"(mhi), [sumlo]"r"(sumlo), [sumhi]"r"(sumhi)
             : "cc");
    __uint128_t result = (static_cast<__uint128_t>(reshi) << 64) | reslo;

    HPBC_CLOCKWORK_POSTCONDITION2(result < modulus);  // __uint128_t guarantees result>=0.
    HPBC_CLOCKWORK_POSTCONDITION2(result ==
                       default_two_times_restricted_unsigned::call(a, modulus));
    return result;
  }
};
# endif

template <>
struct two_times_restricted_unsigned<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t a, std::uint64_t modulus)
  {
    using std::uint64_t;
    using T = std::uint64_t;
    constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(a < modulus);  // uint64_t guarantees a>=0.

    uint64_t sum = static_cast<uint64_t>(a + a);
    uint64_t res;
    __asm__ ("subs %[res], %[sum], %[m] \n\t"        /* res = sum - m */
             "csel %[res], %[sum], %[res], lo \n\t"  /* res = (sum<m) ? sum : res */
             : [res]"=&r"(res)
             : [m]"r"(modulus), [sum]"r"(sum)
             : "cc");
    uint64_t result = res;

    HPBC_CLOCKWORK_POSTCONDITION2(result < modulus);  // uint64_t guarantees result>=0.
    HPBC_CLOCKWORK_POSTCONDITION2(result ==
                       default_two_times_restricted_unsigned::call(a, modulus));
    return result;
  }
};

template <>
struct two_times_restricted_unsigned<std::uint32_t> {
  using U = std::uint32_t;
  HURCHALLA_FORCE_INLINE static U call(U a, U modulus)
  {
    std::uint64_t result = two_times_restricted_unsigned
                                              <std::uint64_t>::call(a, modulus);
    return static_cast<U>(result);
  }
};

// end of inline asm functions for ARM_64
#endif




template <>
struct two_times_restricted_unsigned<std::uint16_t> {
  using U = std::uint16_t;
  HURCHALLA_FORCE_INLINE static U call(U a, U modulus)
  {
    std::uint32_t result = two_times_restricted_unsigned
                                              <std::uint32_t>::call(a, modulus);
    return static_cast<U>(result);
  }
};
template <>
struct two_times_restricted_unsigned<std::uint8_t> {
  using U = std::uint8_t;
  HURCHALLA_FORCE_INLINE static U call(U a, U modulus)
  {
    std::uint32_t result = two_times_restricted_unsigned
                                              <std::uint32_t>::call(a, modulus);
    return static_cast<U>(result);
  }
};




// version for unsigned T
template <typename T, bool = ut_numeric_limits<T>::is_signed>
struct two_times_restricted {
  HURCHALLA_FORCE_INLINE static T call(T a, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    constexpr T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(a < modulus);  // i.e. the input must be prereduced

    return two_times_restricted_unsigned<T>::call(a, modulus);
  }
};

// version for signed T
template <typename T>
struct two_times_restricted<T, true> {
  HURCHALLA_FORCE_INLINE static T call(T a, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::is_signed, "");
    static_assert(static_cast<T>(-1) == ~(static_cast<T>(0)),
                                  "T must use two's complement representation");
    using U = typename extensible_make_unsigned<T>::type;
    static_assert(static_cast<T>(static_cast<U>(static_cast<T>(-1))) ==
                  static_cast<T>(-1), "Casting a signed T value to unsigned and"
                               " back again must result in the original value");
    constexpr U Rdiv2 = static_cast<U>(static_cast<U>(1) <<
                                            (ut_numeric_limits<U>::digits - 1));
    HPBC_CLOCKWORK_PRECONDITION2(0 < modulus && modulus < Rdiv2);
    HPBC_CLOCKWORK_PRECONDITION2(0 <= a && a < modulus);

#if defined(HURCHALLA_AVOID_CSELECT)
    static_assert((static_cast<T>(-1) >> 1) == static_cast<T>(-1),
                          "Arithmetic right shift is required but unavailable");
    T tmp = static_cast<T>(a - modulus);
    HPBC_CLOCKWORK_ASSERT2(tmp < 0);
    tmp = static_cast<T>(tmp + a);
    // if tmp is negative, use a bit mask of all 1s.  Otherwise use all 0s.
    U mask = static_cast<U>(tmp >> ut_numeric_limits<T>::digits);
    U masked_modulus = static_cast<U>(mask & static_cast<U>(modulus));
    U result = static_cast<U>(static_cast<U>(tmp) + masked_modulus);
    HPBC_CLOCKWORK_ASSERT2(result == two_times_restricted_unsigned<U>::call(
                                   static_cast<U>(a), static_cast<U>(modulus)));
#else
    U result = two_times_restricted_unsigned<U>::call(static_cast<U>(a),
                                                       static_cast<U>(modulus));
#endif
    return static_cast<T>(result);
  }
};


}}  // end namespace

#endif
