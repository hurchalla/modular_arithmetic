// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace modular_arithmetic {


template <typename T>
T impl_absolute_value_difference(T a, T b)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION(a >= 0);
    HPBC_PRECONDITION(b >= 0);

    // POSTCONDITION:
    // This function returns absolute_value(a-b).
    T result = (a > b) ? static_cast<T>(a - b) : static_cast<T>(b - a);
    HPBC_POSTCONDITION(result<=a || result<=b);
    return result;
}


// ----------------------------------------------------------------------------
// Non-template function overloads.  By C++ rules, these overloads have first
// priority for function argument matching - the corresponding template version
// above has second priority.
// ----------------------------------------------------------------------------


// MSVC doesn't support inline asm so we skip it.
#if defined(HURCHALLA_ALLOW_INLINE_ASM_ABSDIFF) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
inline uint32_t impl_absolute_value_difference(uint32_t a, uint32_t b)
{
    // Type uint32_t guarantees a>=0 and b>=0.

    uint32_t diff = b - a;

    uint32_t result;
    __asm__ ("subl %[b], %0 \n\t"       /* result = a - b */
             "cmovbl %[diff], %0 \n\t"  /* result = (a < b) ? diff : result */
             : "=&r"(result)
             : "0"(a), [b]"r"(b), [diff]"r"(diff)
             : "cc");

    HPBC_POSTCONDITION2(result<=a || result<=b);
    HPBC_POSTCONDITION2(result== impl_absolute_value_difference<uint32_t>(a,b));
    return result;
}

inline uint64_t impl_absolute_value_difference(uint64_t a, uint64_t b)
{
    // Type uint64_t guarantees a>=0 and b>=0.

    uint64_t diff = b - a;

    uint64_t result;
    __asm__ ("subq %[b], %0 \n\t"       /* result = a - b */
             "cmovbq %[diff], %0 \n\t"  /* result = (a < b) ? diff : result */
             : "=&r"(result)
             : "0"(a), [b]"r"(b), [diff]"r"(diff)
             : "cc");

    HPBC_POSTCONDITION2(result<=a || result<=b);
    HPBC_POSTCONDITION2(result== impl_absolute_value_difference<uint64_t>(a,b));
    return result;
}
#endif


}}  // end namespace

#endif
