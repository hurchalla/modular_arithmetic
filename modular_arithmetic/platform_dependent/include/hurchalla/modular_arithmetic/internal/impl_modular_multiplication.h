
#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATION_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATION_H__INCLUDED


#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

#include <cstdint>
#include <limits>

namespace hurchalla { namespace modular_arithmetic {


/*  Generic (non-platform specific) implementation of the contract for
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
Code review/testing notes:
   Analytical correctness: Appears perfect.
   Empirical correctness: Impossible to test exhaustively, but passed all tests.
*/
#ifndef COMPILE_ERROR_ON_SLOW_MATH
template <typename T>
T impl_modular_multiplication_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(std::numeric_limits<T>::is_integer &&
                 !(std::numeric_limits<T>::is_signed), "");
    precondition(modulus>0);
    precondition(a<modulus);
    precondition(b<modulus);

    T result = 0;
    while (b > 0) {
        namespace ma = ::hurchalla::modular_arithmetic;
        if (b & 1)
            result = ma::modular_addition_prereduced_inputs(result, a, modulus);
        a = ma::modular_addition_prereduced_inputs(a, a, modulus);
        b = b >> 1;
    }
    return result;
}
#else
// cause a compile error if instantiating this (slow) template function
template <typename T>
T impl_modular_multiplication_prereduced_inputs(T a, T b, T modulus) = delete;
#endif // #ifndef COMPILE_ERROR_ON_SLOW_MATH




// -------- PLATFORM SPECIFIC overloads ----------

// Note: these fast nontemplate function overloads get first priority for being
// called (see http://www.gotw.ca/publications/mill17.htm ), when both a
// function overload and the generic template function match the argument type.



#if defined(TARGET_ISA_HAS_DIVIDE) && TARGET_BIT_WIDTH >= 16
inline uint8_t impl_modular_multiplication_prereduced_inputs(uint8_t a,
                                            uint8_t b, uint8_t modulus)
{
    // calculate (a*b)%modulus, guaranteeing no overflow on a*b
    return (uint8_t)((uint16_t)a*(uint16_t)b % (uint16_t)modulus);
}
#endif


#if defined(TARGET_ISA_HAS_DIVIDE) && TARGET_BIT_WIDTH >= 32
inline uint16_t impl_modular_multiplication_prereduced_inputs(uint16_t a,
                                            uint16_t b, uint16_t modulus)
{
    return (uint16_t)((uint32_t)a*(uint32_t)b % (uint32_t)modulus);
}
#endif




// Note: for TARGET_ISA_X86_64 (and TARGET_ISA_X86_32), 32bit mul and div should
// be faster than 64bit mul and div.  For 32bit mul and div, we need access to
// the two-register wide product and dividend, which requires assembly language.

#if defined(_MSC_VER) && _MSC_VER >= 1920 && (defined(TARGET_ISA_X86_64) || defined(TARGET_ISA_X86_32))
// _MSC_VER >= 1920 indicates Visual Studio 2019 or higher. VS2019 (for x86/x64)
// is the first version to support _udiv64 used below.
#  include <immintrin.h>
#  include <intrin.h>
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
    // Note: at the time of this writing (April 2020), the MS documentation of
    // udiv64 is likely incomplete.  udiv64 is almost certainly a direct
    // translation of the x86/x64 assembly "div" instruction (which we want!).
    // See  https://developercommunity.visualstudio.com/content/problem/896815/-udiv128-causes-integer-overflow-and-doesnt-have-a.html
    uint32_t result;
    _udiv64(__emulu(a, b), modulus, &result);
    postcondition2((uint64_t)result==(uint64_t)a*(uint64_t)b%(uint64_t)modulus);
    return result;
}
#elif defined(_MSC_VER) && defined(TARGET_ISA_X86_64)
extern "C" uint32_t modular_multiply_uint32_asm_UID7b5f83fc983(uint32_t a,
                                         uint32_t b, uint32_t modulus) noexcept;
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
    // MSVC doesn't support inline asm for 64 bit code - using my assembly code
    // will likely require an unavoidable function call.  With modern intel CPUs
    // that now have relatively fast DIV instructions, it's unclear that setting
    // up the stack and unwinding the stack and branching to the asm function,
    // all just to use a single 32 bit asm MUL and DIV, will be any faster than
    // simply using 64 bit multiply and divide in standard C++.
    // If getting every bit of performance matters to you, for this application
    // consider using a different compiler which supports inline asm.
    //
    // Absent performance measurements, I've chosen to use simple C++ 64 bit
    // multiply and divide.  If you want asm instead, change to #if 0, below.
  #if 1
    return (uint32_t)((uint64_t)a*(uint64_t)b % (uint64_t)modulus);
  #else
    // The older versions of MSVC don't have the _udiv64 intrinsic.  Since
    // MSVC doesn't support inline asm for 64 bit targets, use an asm function
    uint32_t result = modular_multiply_uint32_asm_UID7b5f83fc983(a, b, modulus);
    postcondition2((uint64_t)result==(uint64_t)a*(uint64_t)b%(uint64_t)modulus);
    return result;
  #endif
}
#elif defined(_MSC_VER) && defined(TARGET_ISA_X86_32)   // inline asm, MS syntax
// Since this is x86 msvc and will use inline asm, we must ensure this function
// isn't using __fastcall or __vectorcall (see
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
    postcondition2((uint64_t)result==(uint64_t)a*(uint64_t)b%(uint64_t)modulus);
    return result;
}
#elif defined(TARGET_ISA_X86_64) || defined(TARGET_ISA_X86_32) // gnu/AT&T style
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
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
    postcondition2((uint64_t)result==(uint64_t)a*(uint64_t)b%(uint64_t)modulus);
    return result;
}
#elif defined(TARGET_ISA_HAS_DIVIDE) && TARGET_BIT_WIDTH >= 64
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
    return (uint32_t)((uint64_t)a*(uint64_t)b % (uint64_t)modulus);
}
#endif




#if defined(_MSC_VER) && _MSC_VER >= 1920 && defined(TARGET_ISA_X86_64)
// _MSC_VER >= 1920 indicates Visual Studio 2019 or higher. VS2019 (for x64)
// is the first version to support _udiv128 used below.
#  include <immintrin.h>
#  include <intrin.h>
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
    postcondition3(result == impl_modular_multiplication_prereduced_inputs
                                               <uint64_t>(a, b, modulus));
    return result;
}
#elif defined(_MSC_VER) && defined(TARGET_ISA_X86_64)
extern "C" uint64_t modular_multiply_uint64_asm_UID7b5f83fc983(uint64_t a,
                                         uint64_t b, uint64_t modulus) noexcept;
inline uint64_t impl_modular_multiplication_prereduced_inputs(uint64_t a,
                                            uint64_t b, uint64_t modulus)
{
    // The older versions of MSVC don't have the _udiv128 intrinsic.  Since
    // MSVC doesn't support inline asm for 64 bit targets, use an asm function
    uint64_t result = modular_multiply_uint64_asm_UID7b5f83fc983(a, b, modulus);
    postcondition3(result == impl_modular_multiplication_prereduced_inputs
                                               <uint64_t>(a, b, modulus));
    return result;
}
#elif defined(TARGET_ISA_X86_64)    // use inline asm with gnu/AT&T syntax
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
    postcondition3(result == impl_modular_multiplication_prereduced_inputs
                                               <uint64_t>(a, b, modulus));
    return result;
}
#elif defined(TARGET_ISA_HAS_DIVIDE) && TARGET_BIT_WIDTH >= 128
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
#elif defined(TARGET_ISA_HAS_DIVIDE) && defined(__SIZEOF_INT128__)
// The macro __SIZEOF_INT128__ indicates if __int128 is supported.  See
// https://stackoverflow.com/questions/16088282/is-there-a-128-bit-integer-in-gcc
inline uint64_t impl_modular_multiplication_prereduced_inputs(uint64_t a,
                                            uint64_t b, uint64_t modulus)
{
    using U = __uint128_t;
    return (uint64_t)((U)a*(U)b % (U)modulus);
}
*/
#endif



}}  // end namespace

#endif  // include guard