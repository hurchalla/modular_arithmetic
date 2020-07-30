// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace modular_arithmetic {


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
    //    But (a-b) overflows whenever b>a, so instead of testing if (a-b < 0),
    //   we test the alternative predicate (a < b).  This gives us our desired
    //   result without any problem of overflow.  So we can and should use:
    //   result = (a < b) ? a-b+modulus : a-b

    // we calculate diff here to encourage the compiler to hoist it out of a
    // loop (assuming this function is inlined inside a loop).
    T diff = static_cast<T>(modulus - b);
    T result = (a < b) ? static_cast<T>(a + diff) : static_cast<T>(a - b);

    HPBC_POSTCONDITION2(0<=result && result<modulus);
    return result;
}


// ----------------------------------------------------------------------------
// Non-template function overloads.  By C++ rules, these overloads have first
// priority for function argument matching - the corresponding template version
// above has second priority.
// ----------------------------------------------------------------------------


// MSVC doesn't support inline asm so we skip it.
#if defined(HURCHALLA_ALLOW_INLINE_ASM_MODSUB) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
inline std::uint32_t impl_modular_subtraction_prereduced_inputs(std::uint32_t a,
                                         std::uint32_t b, std::uint32_t modulus)
{
    using std::uint32_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    // By calculating diff outside of the __asm__, we allow the compiler to loop
    // hoist diff, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint32_t diff = modulus - b;
    uint32_t tmp = a + diff;

    uint32_t result;
    __asm__ ("subl %[b], %0 \n\t"      /* result = a - b */
             "cmovbl %[tmp], %0 \n\t"  /* result = (a < b) ? tmp : result */
             : "=&r"(result)
             : "0"(a), [b]"r"(b), [tmp]"r"(tmp)
             : "cc");

    HPBC_POSTCONDITION2(result<modulus);  // uint32_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
           impl_modular_subtraction_prereduced_inputs<uint32_t>(a, b, modulus));
    return result;
}

inline std::uint64_t impl_modular_subtraction_prereduced_inputs(std::uint64_t a,
                                         std::uint64_t b, std::uint64_t modulus)
{
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    // By calculating diff outside of the __asm__, we allow the compiler to loop
    // hoist diff, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint64_t diff = modulus - b;
    uint64_t tmp = a + diff;

    uint64_t result;
    __asm__ ("subq %[b], %0 \n\t"      /* result = a - b */
             "cmovbq %[tmp], %0 \n\t"  /* result = (a < b) ? tmp : result */
             : "=&r"(result)
             : "0"(a), [b]"r"(b), [tmp]"r"(tmp)
             : "cc");

    HPBC_POSTCONDITION2(result<modulus);  // uint64_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
           impl_modular_subtraction_prereduced_inputs<uint64_t>(a, b, modulus));
    return result;
}
#endif


}}  // end namespace

#endif
