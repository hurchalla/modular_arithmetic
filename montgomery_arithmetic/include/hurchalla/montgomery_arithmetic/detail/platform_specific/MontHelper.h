// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_HELPER_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONT_HELPER_H_INCLUDED


#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <algorithm>

namespace hurchalla { namespace montgomery_arithmetic { namespace detail {



namespace detail_mh {

template <typename T>
HURCHALLA_FORCE_INLINE T default_modsub_canonical_subtrahend(T x, T y, T n)
{
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical

    T diff = static_cast<T>(x - y);
    T result = static_cast<T>(diff + n);
    // encourage compiler to use conditional move via ternary operator
    result = (x>=y) ? diff : result;

    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    return result;
}

template <typename T>
HURCHALLA_FORCE_INLINE T default_modadd_canonical_second_addend(T x, T y, T n)
{
    HPBC_PRECONDITION2(y < n);  // the second addend must be canonical

    // Naively, we would like to set  result = (x+y >= n) ? (x+y-n) : x+y.
    // But x+y possibly could overflow, so instead we will use the equivalent
    // conditional (x >= n-y).  This is safe since due to the precondition y<n,
    // n-y will never overflow.  So we set
    // result = (x >= n-y) ? (x-(n-y)) : x+y
    T tmp = static_cast<T>(n - y);
    T sum = static_cast<T>(x + y);
    T tmp2 = static_cast<T>(x - tmp);
    // encourage compiler to use conditional move via ternary operator
    T result = (x>=tmp) ? tmp2 : sum;

    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    return result;
}

} // end namespace detail_mh



// primary template
template <typename T>
struct MontHelper
{
  // modsub_canonical_subtrahend()  returns x-y (mod n).
  // y must be canonical (meaning: 0 <= y < n).
  // The return value is not necessarily canonical, but it is less than or equal
  // to max(x, n-1).
  static HURCHALLA_FORCE_INLINE T modsub_canonical_subtrahend(T x, T y, T n)
  {
    return detail_mh::default_modsub_canonical_subtrahend(x, y, n);
  }
  // modadd_canonical_second_addend()  returns x+y (mod n).
  // y must be canonical (meaning: 0 <= y < n).
  // The return value is not necessarily canonical, but it is less than or equal
  // to max(x, n-1).
  static HURCHALLA_FORCE_INLINE T modadd_canonical_second_addend(T x, T y, T n)
  {
    return detail_mh::default_modadd_canonical_second_addend(x, y, n);
  }
};


#if defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) && \
      defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

// specialization for uint64_t (for x86_64)
template <>
struct MontHelper<std::uint64_t>
{
  using T = std::uint64_t;
  // Technical note: these functions are out-of-place.  I.e. at the assembly
  // level, the result register is different from any of the input registers.

  // modsub_canonical_subtrahend()  returns x-y (mod n).
  // y must be canonical (meaning: 0 <= y < n).
  // The return value is not necessarily canonical, but it is less than or equal
  // to max(x, n-1).
  static HURCHALLA_FORCE_INLINE T modsub_canonical_subtrahend(T x, T y, T n)
  {
    HPBC_PRECONDITION2(y < n);  // the subtrahend must be canonical

    // We use the "UabcdSD" constraint below so that the LEA instruction doesn't
    // use RBP/EBP or R13 for the base register (which results in slow lea).
    // See impl_modular_subtraction.h for more info.
    T tmp = x;
    T result;
    __asm__ ("subq %[y], %[tmp] \n\t"            /* tmp = x - y */
             "leaq (%[tmp], %[n]), %[res] \n\t"  /* res = tmp + n */
             "cmovaeq %[tmp], %[res] \n\t"       /* res = (x>=y) ? tmp : res */
#  if defined(__INTEL_COMPILER) || defined(__clang__)
             : [tmp]"+&abcdSD"(tmp), [res]"=r"(result)
#  else
             : [tmp]"+&UabcdSD"(tmp), [res]"=r"(result)
#  endif
             : [n]"r"(n), [y]"rm"(y)
             : "cc");
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    HPBC_POSTCONDITION2(result ==
                       detail_mh::default_modsub_canonical_subtrahend(x, y, n));
    return result;
  }

  // modadd_canonical_second_addend()  returns x+y (mod n).
  // y must be canonical (meaning: 0 <= y < n).
  // The return value is not necessarily canonical, but it is less than or equal
  // to max(x, n-1).
  static HURCHALLA_FORCE_INLINE T modadd_canonical_second_addend(T x, T y, T n)
  {
    HPBC_PRECONDITION2(y < n);  // the second addend must be canonical
    T tmp = static_cast<T>(n - y);
    T sum = static_cast<T>(x + y);
    T tmp2 = x;
    __asm__ ("subq %[tmp], %[tmp2] \n\t"     /* tmp2 = x - tmp */
             "cmovaeq %[tmp2], %[sum] \n\t"  /* sum = (x>=tmp) ? tmp2 : sum */
             : [tmp2]"+&r"(tmp2), [sum]"+r"(sum)
             : [tmp]"rm"(tmp)
             : "cc");
    T result = sum;
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    HPBC_POSTCONDITION2(result ==
                    detail_mh::default_modadd_canonical_second_addend(x, y, n));
    return result;
  }
};


// specialization for uint32_t (for x86_64)
// See MontHelper<std::uint64_t> for comments.  Essentially everything written
// there applies here as well.
template <>
struct MontHelper<std::uint32_t>
{
  using T = std::uint32_t;

  static HURCHALLA_FORCE_INLINE T modsub_canonical_subtrahend(T x, T y, T n)
  {
    HPBC_PRECONDITION2(y < n);

    // For info on the "UabcdSD" constraint, see the the uint64_t specialization
    T tmp = x;
    T result;
    __asm__ ("subl %[y], %[tmp] \n\t"             /* tmp = x - y */
             "leal (%q[tmp], %q[n]), %[res] \n\t" /* res = tmp + n */
             "cmovael %[tmp], %[res] \n\t"        /* res = (x>=y) ? tmp : res */
#  if defined(__INTEL_COMPILER) || defined(__clang__)
             : [tmp]"+&abcdSD"(tmp), [res]"=r"(result)
#  else
             : [tmp]"+&UabcdSD"(tmp), [res]"=r"(result)
#  endif
             : [n]"r"(n), [y]"rm"(y)
             : "cc");
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    HPBC_POSTCONDITION2(result ==
                       detail_mh::default_modsub_canonical_subtrahend(x, y, n));
    return result;
  }

  static HURCHALLA_FORCE_INLINE T modadd_canonical_second_addend(T x, T y, T n)
  {
    HPBC_PRECONDITION2(y < n);
    T tmp = static_cast<T>(n - y);
    T sum = static_cast<T>(x + y);
    T tmp2 = x;
    __asm__ ("subl %[tmp], %[tmp2] \n\t"      /* tmp2 = x - tmp */
             "cmovael %[tmp2], %[sum] \n\t"   /* sum = (x>=tmp) ? tmp2 : sum */
             : [tmp2]"+&r"(tmp2), [sum]"+r"(sum)
             : [tmp]"rm"(tmp)
             : "cc");
    T result = sum;
    HPBC_POSTCONDITION2(result <= std::max(x, static_cast<T>(n-1)));
    HPBC_POSTCONDITION2(result ==
                    detail_mh::default_modadd_canonical_second_addend(x, y, n));
    return result;
  }
};

#endif   // defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) &&
         // defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)


}}} // end namespace

#endif
