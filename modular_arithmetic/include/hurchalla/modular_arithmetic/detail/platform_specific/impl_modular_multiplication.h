// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATION_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
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

namespace hurchalla { namespace modular_arithmetic {


/*  Generic (non-platform specific) implementation for
T modular_multiplication_prereduced_inputs(T a, T b, T modulus).
Ideally for best performance, call with a >= b.
Notes:
   This code was adapted from mulmod() at
 http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=primalityTesting
   I've made faster non-template function overloads later in this file when
   possible. They are platform specific, though the uint8_t version is available
   for any 16 bit or greater target; likewise uint16_t, uint32_t, and uint64_t
   versions are available for any platform target with bit width that is at
   least twice the size of the prospective uint_t's bit width.
   Additionally for the x86/x64 ISAs I've made overloads using assembly language
   for uint32_t and uint64_t, later in this file.
   From preliminary investigation I did on ARM, I don't believe ARM would get
   any significant gain by writing assembly language overloads of this function.
   At the time of this writing none of the ARM ISAs appear to provide an
   instruction for division of a 128 bit dividend by a 64 bit divisor (with a 64
   bit quotient); ARM also doesn't seem to have any instruction to divide a 64
   bit dividend by a 32 bit divisor. Not all ARM ISAs even have the standard
   division instructions (32bit by 32bit, or 64bit by 64 bit). ARM does have a
   UMULL instruction though, which does 32bit by 32bit multiplication for
   (effectively) a 64bit result.
Code review notes: Everything appears correct.
*/
template <typename T>
T slow_modular_multiplication(T a, T b, T modulus)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a>=0 && a<modulus);
    HPBC_PRECONDITION2(b>=0 && b<modulus);

    T result = 0;
    while (b > 0) {
        if (b & 1)
            result = modular_addition_prereduced_inputs(result, a, modulus);
        a = modular_addition_prereduced_inputs(a, a, modulus);
        b = static_cast<T>(b >> 1);
    }
    return result;
}




// ---- TEMPLATE versions of impl_modular_multiplication_prereduced_inputs ----


// True for signed integral types *KNOWN* to std::type_traits, otherwise false.
// I.e. true for the integer types that are both signed and native (int, short,
// long, int64_t, etc).
template <typename T>
struct IsSignedAndNativeInteger : std::integral_constant<
            bool,
            std::is_integral<T>::value && std::is_signed<T>::value
        > {};


// Enabled for matching to any T that is *not* a signed-and-native integer type.
// Note that typically, unsigned integer types resolve to one of the platform
// specific nontemplate overloads below instead of to this template (nontemplate
// functions get first priority for potential matches).
template <typename T>
typename std::enable_if<!(IsSignedAndNativeInteger<T>::value), T>::type
#ifdef HURCHALLA_COMPILE_ERROR_ON_SLOW_MATH
  // cause a compile error instead of falling back to the slow template function
  impl_modular_multiplication_prereduced_inputs(T a, T b, T modulus) = delete;
#else
  impl_modular_multiplication_prereduced_inputs(T a, T b, T modulus)
  {
      return slow_modular_multiplication(a, b, modulus);
  }
#endif

// Enabled for matching to any signed-and-native integer type T.  The internal
// call typically resolves to one of the platform specific overloads below.  If
// not, it resolves to the template just above enabled for not-signed-and-native
// types.
template <typename T>
typename std::enable_if<IsSignedAndNativeInteger<T>::value, T>::type
impl_modular_multiplication_prereduced_inputs(T a, T b, T modulus)
{
    HPBC_PRECONDITION2(a>=0 && b>=0 && modulus>0 && a<modulus && b<modulus);
    using U = typename std::make_unsigned<T>::type;
    static_assert(!IsSignedAndNativeInteger<U>::value, "");
    return (T)impl_modular_multiplication_prereduced_inputs((U)a, (U)b,
                                                                (U)modulus);
}




// -------- PLATFORM SPECIFIC nontemplate overloads ----------

// Note: these fast nontemplate function overloads get first priority for being
// called (see http://www.gotw.ca/publications/mill17.htm ), when both one of
// these nontemplate functions and the generic template function match the
// caller's provided argument type(s).


#if !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                   HURCHALLA_TARGET_BIT_WIDTH >= 16
inline uint8_t impl_modular_multiplication_prereduced_inputs(uint8_t a,
                                            uint8_t b, uint8_t modulus)
{
    // Calculate (a*b)%modulus, guaranteeing no overflow on a*b.
    // Note: the C++ standard mandates 'unsigned int' is at least 16 bits wide
    return (uint8_t)((unsigned int)a*(unsigned int)b % (unsigned int)modulus);
}
#endif


#if !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                   HURCHALLA_TARGET_BIT_WIDTH >= 32
inline uint16_t impl_modular_multiplication_prereduced_inputs(uint16_t a,
                                            uint16_t b, uint16_t modulus)
{
    return (uint16_t)((uint32_t)a*(uint32_t)b % (uint32_t)modulus);
}
#endif




// Note: for HURCHALLA_TARGET_ISA_X86_32 and HURCHALLA_TARGET_ISA_X86_64, 32bit
// mul and div is faster (on current and past intel/amd cpus) than 64bit mul
// and div.  To use 32bit mul and div here, we need access to the two-register
// wide product and dividend, which requires assembly language or intrinsics.

#if defined(_MSC_VER) && (_MSC_VER >= 1920) && \
  (defined(HURCHALLA_TARGET_ISA_X86_64) || defined(HURCHALLA_TARGET_ISA_X86_32))
// _MSC_VER >= 1920 indicates Visual Studio 2019 or higher. VS2019 (for x86/x64)
// is the first version to support _udiv64 used below.
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
    // Note: at the time of this writing (April 2020), the MS documentation of
    // udiv64 is likely incomplete.  udiv64 is almost certainly a direct
    // translation of the x86/x64 assembly "div" instruction (which we want).
    // See  https://developercommunity.visualstudio.com/content/problem/896815/-udiv128-causes-integer-overflow-and-doesnt-have-a.html
    uint32_t result;
    _udiv64(__emulu(a, b), modulus, &result);
    HPBC_POSTCONDITION2((uint64_t)result == 
                        (uint64_t)a*(uint64_t)b % (uint64_t)modulus);
    return result;
}
#elif defined(HURCHALLA_ALLOW_INLINE_ASM_MODMUL) && defined(_MSC_VER) && \
      defined(HURCHALLA_TARGET_ISA_X86_32)
// Since this is x86 msvc and we will use inline asm, we must ensure this
// function doesn't use __fastcall or __vectorcall (see
// https://docs.microsoft.com/en-us/cpp/assembler/inline/using-and-preserving-registers-in-inline-assembly ).
// To do this, we declare this function with the __cdecl modifier, which forces
// the function to use cdecl calling convention.  This overrides any potential
// compiler flag that might specify __fastcall or __vectorcall.
inline uint32_t __cdecl impl_modular_multiplication_prereduced_inputs(
                                       uint32_t a, uint32_t b, uint32_t modulus)
{
    uint32_t result;
    __asm {
        mov eax, a
        mul b           ; EDX:EAX = EAX*b; high-order bits of the product in EDX
        div modulus     ; (quotient EAX, remainder EDX) = EDX:EAX/modulus
        mov result, edx ; save the remainder
    }
    HPBC_POSTCONDITION2((uint64_t)result ==
                        (uint64_t)a*(uint64_t)b % (uint64_t)modulus);
    return result;
}
#elif defined(HURCHALLA_ALLOW_INLINE_ASM_MODMUL) && !defined(_MSC_VER) && \
      ( defined(HURCHALLA_TARGET_ISA_X86_64) || \
        defined(HURCHALLA_TARGET_ISA_X86_32) )
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
    // gnu/AT&T style inline asm
    // [input operands]:  EAX = a
    // mull %3:           EDX:EAX = EAX*b; high-order bits of the product in EDX
    // divl %4:           (quotient EAX, remainder EDX) = EDX:EAX/modulus
    // [output operands]: result = EDX
    // [clobber list]:    both mull and divl clobber the FLAGS register ["cc"]
    uint32_t result, dummy;
    __asm__ ("mull %3\n\t"
             "divl %4"
             : "=&d"(result), "=&a"(dummy)
             : "1"(a), "r"(b), "r"(modulus)
             : "cc");
    HPBC_POSTCONDITION2((uint64_t)result ==
                        (uint64_t)a*(uint64_t)b % (uint64_t)modulus);
    return result;
}
#elif !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                     HURCHALLA_TARGET_BIT_WIDTH >= 64
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
    return (uint32_t)((uint64_t)a*(uint64_t)b % (uint64_t)modulus);
}
#endif




#if defined(_MSC_VER) && _MSC_VER >= 1920 && \
          defined(HURCHALLA_TARGET_ISA_X86_64)
// _MSC_VER >= 1920 indicates Visual Studio 2019 or higher. VS2019 (for x64)
// is the first version to support _udiv128 used below.
inline uint64_t impl_modular_multiplication_prereduced_inputs(uint64_t a,
                                            uint64_t b, uint64_t modulus)
{
    // Note: at the time of this writing (April 2020), the MS documentation of
    // udiv128 is likely incomplete.  udiv128 is almost certainly a direct
    // translation of the x64 assembly "div" instruction (which we want!).
    // See  https://developercommunity.visualstudio.com/content/problem/896815/-udiv128-causes-integer-overflow-and-doesnt-have-a.html
    uint64_t productHigh, result;
    uint64_t productLow = _umul128(a, b, &productHigh);
    _udiv128(productHigh, productLow, modulus, &result);
    HPBC_POSTCONDITION3(result == slow_modular_multiplication(a, b, modulus));
    return result;
}
#elif defined(_MSC_VER) && defined(HURCHALLA_TARGET_ISA_X86_64)
extern "C" uint64_t modular_multiply_uint64_asm_UID7b5f83fc983(uint64_t a,
                                         uint64_t b, uint64_t modulus) noexcept;
inline uint64_t impl_modular_multiplication_prereduced_inputs(uint64_t a,
                                            uint64_t b, uint64_t modulus)
{
    // The older versions of MSVC don't have the _udiv128 intrinsic.  Since
    // MSVC doesn't support inline asm for 64 bit targets, use an asm function
    uint64_t result = modular_multiply_uint64_asm_UID7b5f83fc983(a, b, modulus);
    HPBC_POSTCONDITION3(result == slow_modular_multiplication(a, b, modulus));
    return result;
}
#elif defined(HURCHALLA_ALLOW_INLINE_ASM_MODMUL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64)    // inline asm with gnu/AT&T syntax
inline uint64_t impl_modular_multiplication_prereduced_inputs(uint64_t a,
                                            uint64_t b, uint64_t modulus)
{
    // [input operands]:  RAX = a
    // mulq %3:           RDX:RAX = RAX*b; high-order bits of the product in RDX
    // divq %4:           (quotient RAX, remainder RDX) = RDX:RAX/modulus
    // [output operands]: result = RDX
    // [clobber list]:    both mulq and divq clobber the FLAGS register ["cc"]
    uint64_t result, dummy;
    __asm__ ("mulq %3\n\t"
             "divq %4"
             : "=&d"(result), "=&a"(dummy)
             : "1"(a), "r"(b), "r"(modulus)
             : "cc");
    HPBC_POSTCONDITION3(result == slow_modular_multiplication(a, b, modulus));
    return result;
}
#elif !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                    HURCHALLA_TARGET_BIT_WIDTH >= 128
// this is speculative since I don't know of any 128 bit ALUs.
inline uint64_t impl_modular_multiplication_prereduced_inputs(uint64_t a,
                                            uint64_t b, uint64_t modulus)
{
    return (uint64_t)((uint128_t)a*(uint128_t)b % (uint128_t)modulus);
}
/*
// For the next #elif, it's uncertain that division using __uint128_t on a 64bit
// system would be any better than the generic template version of this
// function.
// The code below should be correct as-is. If you wish to try it, you can
// optionally uncomment this section to enable it.
//
#elif !defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) && \
                                            (HURCHALLA_COMPILER_HAS_UINT128_T())
inline uint64_t impl_modular_multiplication_prereduced_inputs(uint64_t a,
                                            uint64_t b, uint64_t modulus)
{
    using U = __uint128_t;
    return (uint64_t)((U)a*(U)b % (U)modulus);
}
*/
#endif



}}  // end namespace

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif


#endif  // include guard
