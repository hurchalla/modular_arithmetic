// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace modular_arithmetic { namespace detail {


template <typename T>
T impl_absolute_value_difference(T a, T b)
{
    namespace ut = hurchalla::util;
    static_assert(ut::ut_numeric_limits<T>::is_integer, "");
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
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_ABSDIFF)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
inline
std::uint32_t impl_absolute_value_difference(std::uint32_t a, std::uint32_t b)
{
    using std::uint32_t;
    // Type uint32_t guarantees a>=0 and b>=0.

    uint32_t diff = b - a;
    uint32_t tmp = a;  // in C++ we prefer not to overwrite an input (a)
    __asm__ ("subl %[b], %[tmp] \n\t"       /* tmp = a - b */
             "cmovbl %[diff], %[tmp] \n\t"  /* tmp = (a < b) ? diff : tmp */
             : [tmp]"+&r"(tmp)
             : [b]"r"(b), [diff]"r"(diff)
             : "cc");
    uint32_t result = tmp;

    HPBC_POSTCONDITION2(result<=a || result<=b);
    HPBC_POSTCONDITION2(result== impl_absolute_value_difference<uint32_t>(a,b));
    return result;
}

inline
std::uint64_t impl_absolute_value_difference(std::uint64_t a, std::uint64_t b)
{
    using std::uint64_t;
    // Type uint64_t guarantees a>=0 and b>=0.

    uint64_t diff = b - a;
    uint64_t tmp = a;  // in C++ we prefer not to overwrite an input (a)
    __asm__ ("subq %[b], %[tmp] \n\t"       /* tmp = a - b */
             "cmovbq %[diff], %[tmp] \n\t"  /* tmp = (a < b) ? diff : tmp */
             : [tmp]"+&r"(tmp)
             : [b]"r"(b), [diff]"r"(diff)
             : "cc");
    uint64_t result = tmp;

    HPBC_POSTCONDITION2(result<=a || result<=b);
    HPBC_POSTCONDITION2(result== impl_absolute_value_difference<uint64_t>(a,b));
    return result;
}
#endif


}}}  // end namespace

#endif
