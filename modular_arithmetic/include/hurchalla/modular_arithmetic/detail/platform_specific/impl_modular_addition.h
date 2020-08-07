// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


template <typename T>
T impl_modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(0<=a && a<modulus);  // i.e. the input must be prereduced
    HPBC_PRECONDITION2(0<=b && b<modulus);  // i.e. the input must be prereduced

    // We want essentially-  result = (a+b < modulus) ? a+b : a+b-modulus
    //   But due to the potential for overflow on a+b, we need to instead test
    //   the alternative predicate (a < modulus-b), which gives us our desired
    //   result without any problem of overflow.  So we can and should use:
    //   result = (a < modulus-b) ? a+b : a+b-modulus
    T tmp = static_cast<T>(modulus - b);
    T result = (a < tmp) ? static_cast<T>(a + b) : static_cast<T>(a - tmp);

    HPBC_POSTCONDITION2(static_cast<T>(0)<=result && result<modulus);
    return result;
}


// ----------------------------------------------------------------------------
// Non-template function overloads.  By C++ rules, these overloads have first
// priority for function argument matching - the corresponding template version
// above has second priority.
// ----------------------------------------------------------------------------


// MSVC doesn't support inline asm, so we skip it.
#if defined(HURCHALLA_ALLOW_INLINE_ASM_MODADD) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
inline std::uint32_t impl_modular_addition_prereduced_inputs(std::uint32_t a,
                                         std::uint32_t b, std::uint32_t modulus)
{
    using std::uint32_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    // By calculating tmp outside of the __asm__, we allow the compiler to
    // potentially loop hoist tmp, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint32_t tmp = modulus - b;
    uint32_t sum = a + b;
    uint32_t tmp2 = a;  // in C++ we prefer not to overwrite an input (a)
    __asm__ ("subl %[tmp], %[tmp2] \n\t"    /* tmp2 = a - tmp */
             "cmovbl %[sum], %[tmp2] \n\t"  /* tmp2 = (a<tmp) ? sum : tmp2 */
             : [tmp2]"+&r"(tmp2)
             : [tmp]"rm"(tmp), [sum]"r"(sum)
             : "cc");
    uint32_t result = tmp2;

    HPBC_POSTCONDITION2(result<modulus);  // uint32_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
              impl_modular_addition_prereduced_inputs<uint32_t>(a, b, modulus));
    return result;
}

inline std::uint64_t impl_modular_addition_prereduced_inputs(std::uint64_t a,
                                         std::uint64_t b, std::uint64_t modulus)
{
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    // By calculating tmp outside of the __asm__, we allow the compiler to
    // potentially loop hoist tmp, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint64_t tmp = modulus - b;
    uint64_t sum = a + b;
    uint64_t tmp2 = a;  // in C++ we prefer not to overwrite an input (a)
    __asm__ ("subq %[tmp], %[tmp2] \n\t"    /* tmp2 = a - tmp */
             "cmovbq %[sum], %[tmp2] \n\t"  /* tmp2 = (a<tmp) ? sum : tmp2 */
             : [tmp2]"+&r"(tmp2)
             : [tmp]"rm"(tmp), [sum]"r"(sum)
             : "cc");
    uint64_t result = tmp2;

    HPBC_POSTCONDITION2(result<modulus);  // uint64_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
              impl_modular_addition_prereduced_inputs<uint64_t>(a, b, modulus));
    return result;
}
#endif


}}  // end namespace

#endif
