// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace detail {


// With regard to the LowlatencyTag vs. the LowuopsTag functions:
// If neither 'b' nor 'modulus' was set/modified recently before the call of
// this modular addition function, then the LowlatencyTag versions will likely
// provide lower latency than the LowuopsTag versions.  Note that LowlatencyTag
// will typically use more uops and create more pressure on the ALU than
// LowuopsTag, unless the compiler can loop hoist the extra instruction(s)
// involving 'b' and 'modulus'.

// Fyi: the purpose of having structs with static member functions is to
// disallow ADL and to make specializations simple and easy.

// primary template for default implementation
template <class PTAG>
struct default_impl_modadd_unsigned {
};

// --- Version #0 (for low uops and low ALU use) ---
template <>
struct default_impl_modadd_unsigned<LowuopsTag> {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // i.e. the input must be prereduced
    HPBC_PRECONDITION2(b<modulus);  // i.e. the input must be prereduced

    T sum = static_cast<T>(a + b);
    T result = static_cast<T>(sum - modulus);
      // result = (sum < modulus) ? sum : result;
    result = ::hurchalla::conditional_select(sum < modulus, sum, result);

    HPBC_POSTCONDITION2(static_cast<T>(0) <= result && result < modulus);
    return result;
  }
};

#if !defined(__clang__)
// Note: only clang (on x64 and arm32/64) seems to produce optimal code from
// the theoretically preferable Version 2 further below.  And so for gcc(x64 and
// arm32/64 and risc-V), MSVC (x64 and arm64), and icc, Version #1 here tends to
// compile to better machine code in practice.  Gcc in particular tends to
// generate conditional branches for Version 2, which we don't want.

// --- LowlatencyTag Version #1 ---
// This is a relatively straightforward and easy to understand version.
//
// However, on x86, Version #1's use of (modulus - b) will often require two
// uops because the compiler will see 'modulus' stays constant over a long
// period of time, and hence it will likely choose not to overwrite its
// register - necessitating an extra uop for a copy of the register.  This may
// be inconsequential, for example if the compiler is able to hoist both the
// copy and the subtraction out of a loop, but if 'b' does not stay constant
// within a loop, then hoisting is not possible and we would expect Version #1
// to require 2 uops for the subtraction in a perhaps critical loop.
// In contrast, Version #2 uses (b - modulus) which typically requires only one
// uop when 'b' does not stay constant.  This fact is why we (in theory)
// slightly prefer Version #2, although the results will vary from compiler to
// compiler.  Note that if 'b' stays constant for a period of time, or if we are
// compiling for ARM, we would generally expect total uops to be the same
// between the two function versions.  We would generally expect the latency of
// the two versions to be the same, but as always this depends on whether the
// compiler generates good (or not good) machine code.
template <>
struct default_impl_modadd_unsigned<LowlatencyTag> {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // i.e. the input must be prereduced
    HPBC_PRECONDITION2(b<modulus);  // i.e. the input must be prereduced

    // We want essentially-  result = (a+b < modulus) ? a+b : a+b-modulus
    //   But due to the potential for overflow on a+b, we need to instead test
    //   the alternative predicate (a < modulus-b), which gives us our desired
    //   result without any problem of overflow.  So we can and should use:
    //   result = (a < modulus-b) ? a+b : a+b-modulus
    T tmp = static_cast<T>(modulus - b);
    T sum = static_cast<T>(a + b);
    T result = static_cast<T>(a - tmp);
      // result = (a < tmp) ? sum : result
    result = ::hurchalla::conditional_select(a < tmp, sum, result);

    HPBC_POSTCONDITION2(static_cast<T>(0) <= result && result < modulus);
    return result;
  }
};
#else

// --- LowlatencyTag Version #2 ---
// This is a more difficult to understand default implementation version.  The
// proof of this function's correctness is given by the theorem in the comments
// at the end of this file.  See the notes at the end of those comments to
// understand the implementation details.
template <>
struct default_impl_modadd_unsigned<LowlatencyTag> {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // the input must be prereduced
    HPBC_PRECONDITION2(b<modulus);  // the input must be prereduced

    T sum = static_cast<T>(a + b);
    T tmp = static_cast<T>(b - modulus);
    T result = static_cast<T>(a + tmp);
      // result = (result >= a) ? sum : result
    result = ::hurchalla::conditional_select(result >= a, sum, result);

    HPBC_POSTCONDITION2(static_cast<T>(0) <= result && result < modulus);
    return result;
  }
};
#endif




// primary template
template <typename T, class PTAG>
struct impl_modular_addition_unsigned {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    return default_impl_modadd_unsigned<PTAG>::call(a, b, modulus);
  }
};


// These inline asm functions implement optimizations of the default function
// versions #0 and #2 (above), for LowuopsTag and LowlatencyTag respectively.

// MSVC doesn't support inline asm, so we skip it.

#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MODADD)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)


// These LowlatencyTag functions contain the calculation "b - modulus".  If
// neither 'b' nor 'modulus' was recently set/modified, then "b - modulus"
// will usually be calculated at the same time as earlier work by the CPU, or
// in a loop it could potentially be loop hoisted by the compiler.  Either way,
// this potentially allows lower latency than the LowuopsTag version.
template <>
struct impl_modular_addition_unsigned<std::uint32_t, LowlatencyTag> {
  HURCHALLA_FORCE_INLINE static
  std::uint32_t call(std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
  {
    using std::uint32_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    // By calculating tmp outside of the __asm__, we allow the compiler to
    // potentially loop hoist tmp, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    // Even without loop hoisting, tmp can potentially be calculated at the same
    // time as earlier work by the CPU via instruction level parallelism, if
    // we assume that neither b nor modulus was recently modified.
    uint32_t sum = static_cast<uint32_t>(a + b);
    uint32_t tmp = static_cast<uint32_t>(b - modulus);
    uint32_t tmp2 = a;  // we prefer not to overwrite an input (a)
    __asm__ ("addl %[tmp], %[tmp2] \n\t"     /* tmp2 = a + tmp */
             "cmovael %[sum], %[tmp2] \n\t"  /* tmp2 = (tmp2>=a) ? sum : tmp2 */
             : [tmp2]"+&r"(tmp2)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [tmp]"r"(tmp), [sum]"r"(sum)
# else
             : [tmp]"rm"(tmp), [sum]"rm"(sum)
# endif
             : "cc");
    uint32_t result = tmp2;

    HPBC_POSTCONDITION2(result < modulus);  // uint32_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
              default_impl_modadd_unsigned<LowlatencyTag>::call(a, b, modulus));
    return result;
  }
};

template <>
struct impl_modular_addition_unsigned<std::uint64_t, LowlatencyTag> {
  HURCHALLA_FORCE_INLINE static
  std::uint64_t call(std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    uint64_t sum = static_cast<uint64_t>(a + b);
    uint64_t tmp = static_cast<uint64_t>(b - modulus);
    uint64_t tmp2 = a;  // we prefer not to overwrite an input (a)
    __asm__ ("addq %[tmp], %[tmp2] \n\t"     /* tmp2 = a + tmp */
             "cmovaeq %[sum], %[tmp2] \n\t"  /* tmp2 = (tmp2>=a) ? sum : tmp2 */
             : [tmp2]"+&r"(tmp2)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [tmp]"r"(tmp), [sum]"r"(sum)
# else
             : [tmp]"rm"(tmp), [sum]"rm"(sum)
# endif
             : "cc");
    uint64_t result = tmp2;

    HPBC_POSTCONDITION2(result < modulus);  // uint64_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
              default_impl_modadd_unsigned<LowlatencyTag>::call(a, b, modulus));
    return result;
  }
};

#ifdef HURCHALLA_ENABLE_INLINE_ASM_128_BIT
template <>
struct impl_modular_addition_unsigned<__uint128_t, LowlatencyTag> {
  HURCHALLA_FORCE_INLINE static
  __uint128_t call(__uint128_t a, __uint128_t b, __uint128_t modulus)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // __uint128_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // __uint128_t guarantees b>=0.

    __uint128_t tmp = static_cast<__uint128_t>(b - modulus);
    __uint128_t sum = static_cast<__uint128_t>(a + b);
    uint64_t alo = static_cast<uint64_t>(a);
    uint64_t ahi = static_cast<uint64_t>(a >> 64);
    uint64_t tmplo = static_cast<uint64_t>(tmp);
    uint64_t tmphi = static_cast<uint64_t>(tmp >> 64);
    uint64_t sumlo = static_cast<uint64_t>(sum);
    uint64_t sumhi = static_cast<uint64_t>(sum >> 64);
    __asm__ ("addq %[tmplo], %[alo] \n\t"    /* tmp2 = a + tmp */
             "adcq %[tmphi], %[ahi] \n\t"
             "cmovaeq %[sumlo], %[alo] \n\t" /* tmp2 = (tmp2>=a) ? sum : tmp2 */
             "cmovaeq %[sumhi], %[ahi] \n\t"
             : [alo]"+&r"(alo), [ahi]"+&r"(ahi)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [tmplo]"r"(tmplo), [tmphi]"r"(tmphi), [sumlo]"r"(sumlo), [sumhi]"r"(sumhi)
# else
             : [tmplo]"rm"(tmplo), [tmphi]"rm"(tmphi), [sumlo]"rm"(sumlo), [sumhi]"rm"(sumhi)
# endif
             : "cc");
    __uint128_t result = (static_cast<__uint128_t>(ahi) << 64) | alo;

    HPBC_POSTCONDITION2(result < modulus);  // __uint128_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
              default_impl_modadd_unsigned<LowlatencyTag>::call(a, b, modulus));
    return result;
  }
};
#endif


// These LowuopsTag versions should have the lowest ALU use (one add, one sub),
// and often the lowest uop count.
template <>
struct impl_modular_addition_unsigned<std::uint32_t, LowuopsTag> {
  HURCHALLA_FORCE_INLINE static
  std::uint32_t call(std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
  {
    using std::uint32_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    uint32_t sum = static_cast<uint32_t>(a + b);
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

    HPBC_POSTCONDITION2(result < modulus);  // uint32_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
                 default_impl_modadd_unsigned<LowuopsTag>::call(a, b, modulus));
    return result;
  }
};

template <>
struct impl_modular_addition_unsigned<std::uint64_t, LowuopsTag> {
  HURCHALLA_FORCE_INLINE static
  std::uint64_t call(std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    uint64_t sum = static_cast<uint64_t>(a + b);
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

    HPBC_POSTCONDITION2(result < modulus);  // uint64_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
                 default_impl_modadd_unsigned<LowuopsTag>::call(a, b, modulus));
    return result;
  }
};

#ifdef HURCHALLA_ENABLE_INLINE_ASM_128_BIT
template <>
struct impl_modular_addition_unsigned<__uint128_t, LowuopsTag> {
  HURCHALLA_FORCE_INLINE static
  __uint128_t call(__uint128_t a, __uint128_t b, __uint128_t modulus)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // __uint128_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // __uint128_t guarantees b>=0.

    __uint128_t sum = static_cast<__uint128_t>(a + b);
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
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [mlo]"r"(mlo), [mhi]"r"(mhi), [tmplo]"r"(tmplo), [tmphi]"r"(tmphi)
# else
             : [mlo]"rm"(mlo), [mhi]"rm"(mhi), [tmplo]"rm"(tmplo), [tmphi]"rm"(tmphi)
# endif
             : "cc");
    __uint128_t result = (static_cast<__uint128_t>(sumhi) << 64) | sumlo;

    HPBC_POSTCONDITION2(result < modulus);  // __uint128_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
              default_impl_modadd_unsigned<LowuopsTag>::call(a, b, modulus));
    return result;
  }
};
#endif

// end of inline asm functions for x86_64
#endif


// You must use either LowlatencyTag or LowuopsTag for PTAG.  See comment at
// top of this file for details.

// version for unsigned T
template <typename T, class PTAG, bool = ut_numeric_limits<T>::is_signed>
struct impl_modular_addition {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    return impl_modular_addition_unsigned<T,PTAG>::call(a, b, modulus);
  }
};

// version for signed T
template <typename T, class PTAG>
struct impl_modular_addition<T, PTAG, true> {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::is_signed, "");
    static_assert(static_cast<T>(-1) == ~(static_cast<T>(0)),
                                  "T must use two's complement representation");
    using U = typename extensible_make_unsigned<T>::type;
    static_assert(static_cast<T>(static_cast<U>(static_cast<T>(-1))) ==
                  static_cast<T>(-1), "Casting a signed T value to unsigned and"
                               " back again must result in the original value");
    HPBC_PRECONDITION2(modulus > 0);
    HPBC_PRECONDITION2(0 <= a && a < modulus);
    HPBC_PRECONDITION2(0 <= b && b < modulus);

#if defined(HURCHALLA_AVOID_CSELECT)
    static_assert((static_cast<T>(-1) >> 1) == static_cast<T>(-1),
                          "Arithmetic right shift is required but unavailable");
    T tmp = static_cast<T>(b - modulus);
    HPBC_ASSERT2(tmp < 0);
    tmp = static_cast<T>(tmp + a);
    // if tmp is negative, use a bit mask of all 1s.  Otherwise use all 0s.
    U mask = static_cast<U>(tmp >> ut_numeric_limits<T>::digits);
    U masked_modulus = static_cast<U>(mask & static_cast<U>(modulus));
    U result = static_cast<U>(static_cast<U>(tmp) + masked_modulus);
    HPBC_ASSERT2(result == impl_modular_addition_unsigned<U,PTAG>::call(
                static_cast<U>(a), static_cast<U>(b), static_cast<U>(modulus)));
#else
    U result = impl_modular_addition_unsigned<U,PTAG>::call(static_cast<U>(a),
                                    static_cast<U>(b), static_cast<U>(modulus));
#endif
    return static_cast<T>(result);
  }
};




// ---------Theorem and proof for Version #2--------
// The constant "R" used below represents the value
// R = 1<<(ut_numeric_limits<T>::digits).  For example, if T is uint64_t,
// then R = 1<<64.
// We'll use a psuedo-cast notation of (Z)x to indicate when we are treating x
// as an infinite precision signed integer - i.e. a member of the set Z of
// mathematical integers.
// The notation "%%" used below should be interpreted as a conceptual modulo
// operator that will always produce a non-negative remainder for the result.
// This is slightly different from the actual C/C++ modulo operator "%", which
// produces a negative remainder if the dividend is negative.

// Theorem:  Require 0 <= a < m, 0 <= b < m, and 0 < m < R; and let
// tmp = ((Z)b - (Z)m) %% R.
//    If ((Z)a + (Z)tmp >= R), then
//        ((Z)a + (Z)b) %% m == a+b-m.
//    else
//        ((Z)a + (Z)b) %% m == a+b.
//
// Proof:
// As a precondition, we know  0 <= b < m.
// [1] Therefore  -(Z)m <= (Z)b - (Z)m < 0
// As a precondition, we know  m < R,  thus
// [2]  -R < -(Z)m
// Combining [1] and [2],  -R < (Z)b - (Z)m < 0.
// Adding R to all parts,  0 < (Z)b - (Z)m + R < R.
// We can see this expression is bound between 0 and R, so
// (Z)b - (Z)m + R == ((Z)b - (Z)m + R) %% R.
// [3] Thus  ((Z)b - (Z)m) %% R == (Z)b - (Z)m + R.
//
// C/C++ (and assembly) performs unsigned addition and subtraction modulo R, and
// of course produces an unsigned non-negative result.
// Therefore, tmp = b - m == ((Z)b - (Z)m) %% R.  Thus 0 <= tmp < R, and by [3],
// [4]  tmp = b - m == (Z)b - (Z)m + R.
//
// We would like to test whether  (Z)a + (Z)b >= (Z)m,  but this test can't be
// directly evaluated in C/C++/assembly due to potential overflow on  a+b.
// However, we can re-express this test as  (Z)b - (Z)m + R >= R - (Z)a,  and
// combining this with [4], this test becomes  (Z)tmp >= R - (Z)a.  We can then
// rearrange this last test into the test  (Z)a + (Z)tmp >= R.
// These tests are all equivalent, and so
// [5]  (Z)a + (Z)tmp >= R  implies  (Z)a + (Z)b >= (Z)m.  And
// [6]  (Z)a + (Z)tmp < R   implies  (Z)a + (Z)b < (Z)m.
// Note that we can easily evaluate the test  (Z)a + (Z)tmp >= R  using
// C/C++/assembly, by performing the addition  a + tmp,  and detecting whether
// or not the add overflows.
//
// As preconditions, we know  0 <= a < m,  0 <= b < m,  and  m < R.
// [7]  Therefore  0 <= (Z)a + (Z)b < (Z)m + (Z)m.
// [8]  Assume  (Z)a + (Z)b >= (Z)m:
//    Then by [7]          (Z)m <= (Z)a + (Z)b < (Z)m + (Z)m,  and
//                            0 <= (Z)a + (Z)b - (Z)m < (Z)m.  Thus,
//    ((Z)a + (Z)b - (Z)m) %% m == (Z)a + (Z)b - (Z)m,  and thus
//           ((Z)a + (Z)b) %% m == (Z)a + (Z)b - (Z)m.
//    Since (Z)m < R,         0 <= (Z)a + (Z)b - (Z)m < R.  Thus,
//    ((Z)a + (Z)b - (Z)m) %% R == (Z)a + (Z)b - (Z)m,  and therefore in C/C++,
//                    a + b - m == (Z)a + (Z)b - (Z)m.  This gives us
//           ((Z)a + (Z)b) %% m == a + b - m.
// [9]  Assume  (Z)a + (Z)b < (Z)m:
//    Then by [7]             0 <= (Z)a + (Z)b < (Z)m.  Thus,
//           ((Z)a + (Z)b) %% m == (Z)a + (Z)b
//    Since (Z)m < R,         0 <= (Z)a + (Z)b < R.  Thus,
//           ((Z)a + (Z)b) %% R == (Z)a + (Z)b,  and therefore in C/C++,
//                        a + b == (Z)a + (Z)b.  This gives us
//           ((Z)a + (Z)b) %% m == a + b.
//
// [10]  Combining [5] with [8],  if ((Z)a + (Z)tmp >= R)  then
//                               ((Z)a + (Z)b) %% m == a+b-m.
// [11]  Combining [6] with [9],  if ((Z)a + (Z)tmp < R)  then
//                               ((Z)a + (Z)b) %% m == a+b.
//
// Implementation notes:
// As stated above, in C/C++/assembly, we can test if  (Z)a + (Z)tmp >= R  by
//    detecting if the addition a + tmp  overflows.
// To detect overflow on the addition (a+tmp):
//    In assembly, we can perform the add and then look at the carry flag to see
//    if it overflowed.  In C, we can test if (a+tmp >= a).  If false, then the
//    add overflowed.  Note: with luck the compiler will recognize this C idiom
//    and produce assembly that is the same as what we would write.  Clang
//    seems to do fairly well at this, gcc less well, and icc least well.
//
// Putting together [4], [10], and [11], in C++ we could write
// template <typename T>
// T modular_addition(T a, T b, T m) {
//    static_assert( !(ut_numeric_limits<T>::is_signed), "" );
//    T tmp = b-m;
//    return (a+tmp >= a) ? a+b : a+tmp;
// }


}}  // end namespace

#endif
