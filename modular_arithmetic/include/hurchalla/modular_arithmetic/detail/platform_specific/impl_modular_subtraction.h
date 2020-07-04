// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace modular_arithmetic {


// MSVC doesn't support inline asm, and clang compiles the template to optimal
// asm (with more flexibility), so we skip MSVC and only compile this for clang
// when we want to test inline asm (especially for clang compilation warnings).
#if defined(HURCHALLA_ALLOW_ALL_INLINE_ASM) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER) && \
    ( !defined(__clang__) || defined(HURCHALLA_TEST_INLINE_ASM) )
inline uint32_t impl_modular_subtraction_prereduced_inputs(uint32_t a,
                                                   uint32_t b, uint32_t modulus)
{
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.
    // Note: we want to make sure the LEA instruction doesn't use RBP/EBP or R13
    // for the base register, since that would necessitate a slower form of LEA
    // that has an extra 2 cycles latency and half the throughput of the fast
    // form.  We prevent this by using the "U" constraint which allows only RAX,
    // RCX, RDX, R8, R9, R10, R11 for the Microsoft x64 calling convention, and
    // allows only RAX, RCX, RDX, RSI, RDI, R8, R9, R10, R11 for the System V
    // AMD64 calling convention.  ICC (intel compiler) and clang don't support
    // "U" so we use "Q" for them (this allows rax, rbx, rcx, rdx).
    uint32_t result, dummy;
    __asm__ ("subl %3, %2 \n\t"
             "leal (%q0, %q4), %1 \n\t"
             "cmovbl %1, %0 \n\t"
#if defined(__INTEL_COMPILER) || defined(__clang__)
                 : "=&Q"(result), "=r"(dummy)
#else
                 : "=&U"(result), "=r"(dummy)
#endif
             : "0"(a), "r"(b), "r"(modulus)
             : "cc");
    return result;
}
inline uint64_t impl_modular_subtraction_prereduced_inputs(uint64_t a,
                                                   uint64_t b, uint64_t modulus)
{
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.
    // Note: the issues and solutions with LEA and RBP/EBP/R13 are the same
    // here as for the uint32_t version of this function above.
    uint64_t result, dummy;
    __asm__ ("subq %3, %2 \n\t"
             "leaq (%0, %4), %1 \n\t"
             "cmovbq %1, %0 \n\t"
#if defined(__INTEL_COMPILER) || defined(__clang__)
                 : "=&Q"(result), "=r"(dummy)
#else
                 : "=&U"(result), "=r"(dummy)
#endif
             : "0"(a), "r"(b), "r"(modulus)
             : "cc");
    return result;
}
#endif



// -----------------------------------------------------------------------------
// Template function version.  By C++ rules, this function has second priority
// for function argument matching - the non-template versions above have first
// priority.
// -----------------------------------------------------------------------------

template <typename T>
T impl_modular_subtraction_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a>=0 && a<modulus);  // i.e. the input must be prereduced
    HPBC_PRECONDITION2(b>=0 && b<modulus);  // i.e. the input must be prereduced

    // POSTCONDITION:
    // Let a conceptual "%%" operator represent a modulo operator that always
    // returns a non-negative remainder.
    // This function returns (a-b) %% modulus, performed as if a and b are
    // infinite precision signed ints (and thus as if it is impossible for the
    // subtraction (a-b) to overflow).

    // We want essentially-  result = (a-b < 0) ? a-b+modulus : a-b
    //    But for unsigned type T, (a-b < 0) is always false.  So instead we use
    T tmp = static_cast<T>(a-b);
    T result = (a<b) ? static_cast<T>(modulus+tmp) : tmp;
    return result;
}


}}  // end namespace

#endif
