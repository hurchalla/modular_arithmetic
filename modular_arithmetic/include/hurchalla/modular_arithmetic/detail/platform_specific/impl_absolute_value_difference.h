// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_ABSOLUTE_VALUE_DIFFERENCE_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// note: uses a static member function to disallow ADL.
struct default_impl_absdiff {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
#if 0
    T result = (a > b) ? static_cast<T>(a - b) : static_cast<T>(b - a);
#else
    T result = static_cast<T>(b - a);
    HURCHALLA_CMOV(a > b, result, static_cast<T>(a - b));
#endif
    // POSTCONDITION:
    // This function returns absolute_value(a-b).
    HPBC_POSTCONDITION(result<=a || result<=b);
    return result;
  }
};


// primary template
template <typename T>
struct impl_absolute_value_difference {
  HURCHALLA_FORCE_INLINE static T call(T a, T b)
  {
    return default_impl_absdiff::call(a, b);
  }
};


// MSVC doesn't support inline asm so we skip it.
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_ABSDIFF)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

template <>
struct impl_absolute_value_difference<std::uint32_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint32_t call(std::uint32_t a, std::uint32_t b)
  {
    using std::uint32_t;
    uint32_t diff = b - a;
    uint32_t tmp = a;  // we prefer not to overwrite an input (a)
    __asm__ ("subl %[b], %[tmp] \n\t"       /* tmp = a - b */
             "cmovbl %[diff], %[tmp] \n\t"  /* tmp = (a < b) ? diff : tmp */
             : [tmp]"+&r"(tmp)
             : [b]"r"(b), [diff]"r"(diff)
             : "cc");
    uint32_t result = tmp;

    HPBC_POSTCONDITION2(result<=a || result<=b);
    HPBC_POSTCONDITION2(result == default_impl_absdiff::call(a, b));
    return result;
  }
};

template <>
struct impl_absolute_value_difference<std::uint64_t> {
  HURCHALLA_FORCE_INLINE
  static std::uint64_t call(std::uint64_t a, std::uint64_t b)
  {
    using std::uint64_t;
    uint64_t diff = b - a;
    uint64_t tmp = a;  // we prefer not to overwrite an input (a)
    __asm__ ("subq %[b], %[tmp] \n\t"       /* tmp = a - b */
             "cmovbq %[diff], %[tmp] \n\t"  /* tmp = (a < b) ? diff : tmp */
             : [tmp]"+&r"(tmp)
             : [b]"r"(b), [diff]"r"(diff)
             : "cc");
    uint64_t result = tmp;

    HPBC_POSTCONDITION2(result<=a || result<=b);
    HPBC_POSTCONDITION2(result == default_impl_absdiff::call(a, b));
    return result;
  }
};

#endif


}}  // end namespace

#endif
