// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATION_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <type_traits>
#if defined(_MSC_VER)
#  include <immintrin.h>
#  include <intrin.h>
#endif

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#endif

// Please consider using the montgomery_arithmetic MontgomeryForm's multiply
// function instead of the standard modular multiplication of this file, if you
// are heavily using modular multiplication in your code, or if you are in
// danger of (or certain of) getting the slow_modular_multiplication function
// below.  In such cases, there's a good chance that montgomery multiplication
// will greatly improve your code's performance.

// By default, if an inline asm modmult function is available we use it unless
// explicitly disallowed via HURCHALLA_DISALLOW_INLINE_ASM_MODMUL.  We do this
// by default because at least for x86, the asm version is many times faster
// than the non-asm version, and because the inline asm is relatively simple
// (though all inline asm is very difficult to verify as correct).
// Presumably if all inline asm is enabled, asm modmul should not be disallowed.
#ifdef HURCHALLA_ALLOW_INLINE_ASM_ALL
#  undef HURCHALLA_DISALLOW_INLINE_ASM_MODMUL
#endif

namespace hurchalla { namespace detail {


// Slow implementation that works for all compilers and architectures.
// Ideally for best performance, call with a >= b.
// Credit: this code was adapted from mulmod() at
//   http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=primalityTesting
// Note: uses a static member function to disallow ADL.
struct slow_modular_multiplication {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);
    HPBC_PRECONDITION2(b<modulus);

    T result = 0;
    while (b > 0) {
        namespace hc = ::hurchalla;
        T tmp = hc::modular_addition_prereduced_inputs(result, a, modulus);
          // result = (b&1) ? tmp : result
        result = hc::conditional_select((b & 1), tmp, result);
        a = hc::modular_addition_prereduced_inputs(a, a, modulus);
        b = static_cast<T>(b >> 1);
    }
    return result;
  }
};


// Note that there are fast specializations of impl_modular_multiplication
// available below in this file whenever possible. They are platform specific,
// though the uint8_t version is available for any 16 bit or greater target;
// likewise uint16_t, uint32_t, and uint64_t versions are available for any
// platform target with bit width that is at least twice the size of the
// prospective uint_t's bit width.
// Additionally for the x86/x64 ISAs I've made versions using assembly language
// for uint32_t and uint64_t.
// From preliminary investigation I did on ARM, I don't believe ARM would get
// any significant gain by writing assembly language versions (though like any
// ISA it benefits from the non-asm specializations, if available for the type
// T). At the time of this writing none of the ARM ISAs appear to provide an
// instruction for division of a 128 bit dividend by a 64 bit divisor (with a 64
// bit quotient); ARM also doesn't seem to have any instruction to divide a 64
// bit dividend by a 32 bit divisor. Not all ARM ISAs even have the standard
// division instructions (32bit by 32bit, or 64bit by 64 bit). ARM does have a
// UMULL instruction though, which does 32bit by 32bit multiplication for
// (effectively) a 64bit result. ARM64 also has UMULH for the high 64 bits of
// a 64bit x 64bit -> 128 bit multiply (for MSVC this is instrinsic __umulh).


// primary template
// In most cases a type will match and resolve to one of the specializations
// below, instead of this primary template.
template <typename T>
struct impl_modular_multiplication {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return true; }
#ifdef HURCHALLA_COMPILE_ERROR_ON_SLOW_MATH
  // cause compile error instead of falling back to slow modular multiplication.
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus) = delete;
#else
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    return slow_modular_multiplication::call(a, b, modulus);
  }
#endif
};


#if !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                   HURCHALLA_TARGET_BIT_WIDTH >= 16
template <> struct impl_modular_multiplication<std::uint8_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static
  std::uint8_t call(std::uint8_t a, std::uint8_t b, std::uint8_t modulus)
  {
    // Calculate (a*b)%modulus, guaranteeing no overflow on a*b.
    // We use safely_promote_unsigned<T> to avoid undefined behavior.  See
    // https://jeffhurchalla.com/2019/01/16/c-c-surprises-and-undefined-behavior-due-to-unsigned-integer-promotion/
    using P = safely_promote_unsigned<std::uint16_t>::type;
    return (std::uint8_t)((P)a*(P)b % (P)modulus);
  }
};
#endif

#if !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                   HURCHALLA_TARGET_BIT_WIDTH >= 32
template <> struct impl_modular_multiplication<std::uint16_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static
  std::uint16_t call(std::uint16_t a, std::uint16_t b, std::uint16_t modulus)
  {
    using P = safely_promote_unsigned<std::uint32_t>::type;
    return (std::uint16_t)((P)a*(P)b % (P)modulus);
  }
};
#endif


// Note: for HURCHALLA_TARGET_ISA_X86_32 and HURCHALLA_TARGET_ISA_X86_64, 32bit
// mul and div is faster (on current and past intel/amd cpus) than 64bit mul
// and div.  To use 32bit mul and div here, we need access to the two-register
// wide product and dividend, which requires assembly language or intrinsics.

#if defined(_MSC_VER) && (_MSC_VER >= 1920) && \
  (defined(HURCHALLA_TARGET_ISA_X86_64) || defined(HURCHALLA_TARGET_ISA_X86_32))
// _MSC_VER >= 1920 indicates Visual Studio 2019 or higher. VS2019 (for x86/x64)
// is the first version to support _udiv64 used below.
template <> struct impl_modular_multiplication<std::uint32_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint32_t call(
                        std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
  {
    using P = safely_promote_unsigned<std::uint64_t>::type;
    std::uint32_t result;
    // Note: at the time of this writing (April 2020), the MS documentation of
    // udiv64 is likely incomplete.  udiv64 is almost certainly a direct
    // translation of the x86/x64 assembly "div" instruction (which we want).
    // See  https://developercommunity.visualstudio.com/content/problem/896815/-udiv128-causes-integer-overflow-and-doesnt-have-a.html
    _udiv64(__emulu(a, b), modulus, &result);
    HPBC_POSTCONDITION2((P)result == (P)a*(P)b % (P)modulus);
    return result;
  }
};
#elif !defined(HURCHALLA_DISALLOW_INLINE_ASM_MODMUL) && defined(_MSC_VER) && \
      defined(HURCHALLA_TARGET_ISA_X86_32)
// Since this is x86 msvc and we will use inline asm, we must ensure this
// function doesn't use __fastcall or __vectorcall (see
// https://docs.microsoft.com/en-us/cpp/assembler/inline/using-and-preserving-registers-in-inline-assembly ).
// To do this, we declare the function with the __cdecl modifier, which forces
// the function to use cdecl calling convention.  This overrides any potential
// compiler flag that might specify __fastcall or __vectorcall.
template <> struct impl_modular_multiplication<std::uint32_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint32_t __cdecl call(
                        std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
  {
    using P = safely_promote_unsigned<std::uint64_t>::type;
    std::uint32_t result;
    __asm {
        mov eax, a
        mul b           ; EDX:EAX = EAX*b; high-order bits of the product in EDX
        div modulus     ; (quotient EAX, remainder EDX) = EDX:EAX/modulus
        mov result, edx ; save the remainder
    }
    HPBC_POSTCONDITION2((P)result == (P)a*(P)b % (P)modulus);
    return result;
  }
};
#elif !defined(HURCHALLA_DISALLOW_INLINE_ASM_MODMUL) && !defined(_MSC_VER) && \
      ( defined(HURCHALLA_TARGET_ISA_X86_64) || \
        defined(HURCHALLA_TARGET_ISA_X86_32) )
template <> struct impl_modular_multiplication<std::uint32_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint32_t call(
                        std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
  {
    using P = safely_promote_unsigned<std::uint64_t>::type;
    // gnu/AT&T style inline asm
    // [inout operands]:  EAX = tmp  (inout lets us safely overwrite EAX)
    // mull %[b]:         EDX:EAX = EAX*b; high-order bits of the product in EDX
    // divl %[m]:         (quotient EAX, remainder EDX) = EDX:EAX/modulus
    // [output operands]: result = EDX
    // [clobber list]:    both mull and divl clobber the FLAGS register ["cc"]
    std::uint32_t result;
    std::uint32_t tmp = a;  // we prefer not to overwrite an input (a)
    __asm__ ("mull %[b] \n\t"
             "divl %[m] \n\t"
             : "=&d"(result), "+&a"(tmp)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [b]"r"(b), [m]"r"(modulus)
# else
             : [b]"rm"(b), [m]"rm"(modulus)
# endif
             : "cc");
    HPBC_POSTCONDITION2((P)result == (P)a*(P)b % (P)modulus);
    return result;
  }
};
#elif !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                     HURCHALLA_TARGET_BIT_WIDTH >= 64
template <> struct impl_modular_multiplication<std::uint32_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint32_t call(
                        std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
  {
    using P = safely_promote_unsigned<std::uint64_t>::type;
    return (std::uint32_t)((P)a*(P)b % (P)modulus);
  }
};
#endif


#if defined(_MSC_VER) && _MSC_VER >= 1920 && \
          defined(HURCHALLA_TARGET_ISA_X86_64)
// _MSC_VER >= 1920 indicates Visual Studio 2019 or higher. VS2019 (for x64)
// is the first version to support _udiv128 used below.
template <> struct impl_modular_multiplication<std::uint64_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint64_t call(
                        std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    // Note: at the time of this writing (April 2020), the MS documentation of
    // udiv128 is likely incomplete.  udiv128 is almost certainly a direct
    // translation of the x64 assembly "div" instruction (which we want!).
    // See  https://developercommunity.visualstudio.com/content/problem/896815/-udiv128-causes-integer-overflow-and-doesnt-have-a.html
    std::uint64_t productHigh, result;
    std::uint64_t productLow = _umul128(a, b, &productHigh);
    _udiv128(productHigh, productLow, modulus, &result);
    HPBC_POSTCONDITION3(result ==
                              slow_modular_multiplication::call(a, b, modulus));
    return result;
  }
};
#elif defined(_MSC_VER) && defined(HURCHALLA_TARGET_ISA_X86_64)
extern "C" std::uint64_t modular_multiply_uint64_asm_UID7b5f83fc983(
              std::uint64_t a, std::uint64_t b, std::uint64_t modulus) noexcept;
template <> struct impl_modular_multiplication<std::uint64_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint64_t call(
                        std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    using std::uint64_t;
    // The older versions of MSVC don't have the _udiv128 intrinsic.  Since
    // MSVC doesn't support inline asm for 64 bit targets, use an asm function
    uint64_t result = modular_multiply_uint64_asm_UID7b5f83fc983(a, b, modulus);
    HPBC_POSTCONDITION3(result ==
                              slow_modular_multiplication::call(a, b, modulus));
    return result;
  }
};
#elif !defined(HURCHALLA_DISALLOW_INLINE_ASM_MODMUL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64)    // inline asm with gnu/AT&T syntax
template <> struct impl_modular_multiplication<std::uint64_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint64_t call(
                        std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    // [inout operands]:  RAX = tmp  (inout lets us safely overwrite RAX)
    // mulq %[b]:         RDX:RAX = RAX*b; high-order bits of the product in RDX
    // divq %[m]:         (quotient RAX, remainder RDX) = RDX:RAX/modulus
    // [output operands]: result = RDX
    // [clobber list]:    both mulq and divq clobber the FLAGS register ["cc"]
    std::uint64_t result;
    std::uint64_t tmp = a;  // we prefer not to overwrite an input (a)
    __asm__ ("mulq %[b] \n\t"
             "divq %[m] \n\t"
             : "=&d"(result), "+&a"(tmp)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [b]"r"(b), [m]"r"(modulus)
# else
             : [b]"rm"(b), [m]"rm"(modulus)
# endif
             : "cc");
    HPBC_POSTCONDITION3(result ==
                              slow_modular_multiplication::call(a, b, modulus));
    return result;
  }
};
#elif !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                    HURCHALLA_TARGET_BIT_WIDTH >= 128
// this is speculative since I don't know of any 128 bit ALUs.
template <> struct impl_modular_multiplication<std::uint64_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint64_t call(
                        std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    using P = safely_promote_unsigned<std::uint128_t>::type;
    return (std::uint64_t)((P)a*(P)b % (P)modulus);
  }
};
#elif !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                                            (HURCHALLA_COMPILER_HAS_UINT128_T())
// It's uncertain that division using __uint128_t on a 64bit system would be any
// better than letting ourselves fall back to the primary template
// impl_modular_multiplication that uses slow_modular_multiplication.  The code
// below should be correct as-is. If you wish to try it, you can optionally
// uncomment the following section.
/*
template <> struct impl_modular_multiplication<std::uint64_t> {
  HURCHALLA_FORCE_INLINE static constexpr bool has_slow_perf() { return false; }
  HURCHALLA_FORCE_INLINE static std::uint64_t call(
                        std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    using P = safely_promote_unsigned<__uint128_t>::type;
    return (std::uint64_t)((P)a*(P)b % (P)modulus);
  }
};
*/
#endif


}}  // end namespace

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif


#endif  // include guard
