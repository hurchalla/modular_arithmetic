
#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATION_H__INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_MULTIPLICATION_H__INCLUDED


#include "hurchalla/platform_specific/defines_target_platform.h"
#include "hurchalla/platform_specific/defines_compiler.h"

#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"

#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


/*  Generic (non-platform specific) implementation of the contract for
modular_multiplication_prereduced_inputs(T a, T b, T modulus).
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
   From preliminary investigation I did on ARM, I don't believe there would be
   any significant gain by writing assembly language versions of this function.
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
template <typename T>
inline T impl_modular_multiplication_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(std::is_unsigned<T>::value);  //T unsigned integral type
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
// be faster than 64bit mul and div.  The 32bit mul and div need a two-register
// wide product and dividend, which requires assembly language.
#if defined(_MSC_VER) && defined(TARGET_ISA_X86_64)
extern "C" uint32_t modular_multiply_uint32_asm_UID7b5f83fc983(uint32_t a,
                                             uint32_t b, uint32_t modulus);
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
    // MSVC doesn't support inline asm for 64 bit targets, so use asm function
    uint32_t result = modular_multiply_uint32_asm_UID7b5f83fc983(a, b, modulus);
    postcondition3((uint64_t)result==(uint64_t)a*(uint64_t)b%(uint64_t)modulus);
    return result;
}
#elif defined(_MSC_VER) && defined(TARGET_ISA_X86_32)   // inline asm, MS syntax
// Since this is x86 msvc and will use inline asm, we must ensure this function
// isn't using __fastcall or __vectorcall (see https://docs.microsoft.com/en-us/cpp/assembler/inline/using-and-preserving-registers-in-inline-assembly ).
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
    postcondition3((uint64_t)result==(uint64_t)a*(uint64_t)b%(uint64_t)modulus);
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
    postcondition3((uint64_t)result==(uint64_t)a*(uint64_t)b%(uint64_t)modulus);
    return result;
}
#elif defined(TARGET_ISA_HAS_DIVIDE) && TARGET_BIT_WIDTH >= 64
inline uint32_t impl_modular_multiplication_prereduced_inputs(uint32_t a,
                                            uint32_t b, uint32_t modulus)
{
    return (uint32_t)((uint64_t)a*(uint64_t)b % (uint64_t)modulus);
}
#endif



#if defined(_MSC_VER) && defined(TARGET_ISA_X86_64)
extern "C" uint64_t modular_multiply_uint64_asm_UID7b5f83fc983(uint64_t a,
                                             uint64_t b, uint64_t modulus);
inline uint64_t impl_modular_multiplication_prereduced_inputs(uint64_t a,
                                            uint64_t b, uint64_t modulus)
{
    // MSVC doesn't support inline asm for 64 bit targets, so use asm function
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
#endif


}}  // end namespace

#endif  // include guard
