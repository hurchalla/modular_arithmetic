// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace modular_arithmetic {


// MSVC doesn't support inline asm so we skip it.
#if defined(HURCHALLA_ALLOW_INLINE_ASM_MODSUB) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
inline uint32_t impl_modular_subtraction_prereduced_inputs(uint32_t a,
                                                   uint32_t b, uint32_t modulus)
{
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    // By calculating diff outside of the __asm__, we allow the compiler to loop
    // hoist diff, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint32_t diff = modulus - b;
    uint32_t tmp = diff + a;

    uint32_t result;
    __asm__ ("subl %[b], %0 \n\t"      /* result = a - b */
             "cmovbl %[tmp], %0 \n\t"  /* result = (result>=a) ? sum : result */
             : "=&r"(result)
             : "0"(a), [b]"r"(b), [tmp]"r"(tmp)
             : "cc");

    HPBC_POSTCONDITION2(result<modulus);  // uint32_t guarantees result>=0.
    return result;
}
inline uint64_t impl_modular_subtraction_prereduced_inputs(uint64_t a,
                                                   uint64_t b, uint64_t modulus)
{
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    // By calculating diff outside of the __asm__, we allow the compiler to loop
    // hoist diff, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint64_t diff = modulus - b;
    uint64_t tmp = diff + a;

    uint64_t result;
    __asm__ ("subq %[b], %0 \n\t"      /* result = a - b */
             "cmovbq %[tmp], %0 \n\t"  /* result = (result>=a) ? sum : result */
             : "=&r"(result)
             : "0"(a), [b]"r"(b), [tmp]"r"(tmp)
             : "cc");

    HPBC_POSTCONDITION2(result<modulus);  // uint64_t guarantees result>=0.
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

    HPBC_POSTCONDITION2(0<=result && result<modulus);
    return result;
}


}}  // end namespace

#endif
