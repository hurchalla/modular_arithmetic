// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// note: uses a static member function to disallow ADL.
struct default_modsub {
  template <typename T>
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // i.e. the input must be prereduced
    HPBC_PRECONDITION2(b<modulus);  // i.e. the input must be prereduced

    // POSTCONDITION:
    // Let a conceptual "%%" operator represent a modulo operator that always
    // returns a non-negative remainder.
    // This function returns (a-b) %% modulus, performed as if a and b are
    // infinite precision signed ints (and thus as if it is impossible for the
    // subtraction (a-b) to overflow).

    // We want essentially-  result = (a-b >= 0) ? a-b : a-b+modulus
    //    But (a-b) overflows whenever b>a, so instead of testing if (a-b >= 0),
    //   we test the alternative predicate (a >= b).  This gives us our desired
    //   result without any problem of overflow.  So we can and should use:
    //   result = (a>=b) ? a-b : a-b+modulus

    // This implementation is designed for low uop count and low register use.
    // An implementation is possible with expected lower latency and higher uop
    // count and higher register use, but it's not preferred (see the comments
    // below for that alternative).
    T diff = static_cast<T>(a - b);
    T result = static_cast<T>(diff + modulus);
# if 0
    result = (a >= b) ? diff : result;
# else
    HURCHALLA_CMOV(a >= b, result, diff);
# endif

    HPBC_POSTCONDITION2(0<=result && result<modulus);
    return result;
  }
};
// Implementation note: if the compiler generates optimal assembly/machine code,
// then the above function will have 3 cycles latency (on x86).  One cycle for
// the subtraction, one for the addition, and finally one for a cmov.
// In principle we could create an alternative modular_subtraction function
// which has only 2 cycles latency in ideal situations (and 3 cycles in non-
// ideal situations) but that requires as many as two extra uops.  This
// alternative and sometimes lower latency code looks like this:
//    T diff = n - b;
//    T tmp = a + diff;           // ideally a LEA instruction
//    a = a - b;
//    tmp = (a >= b) ? a : tmp;   // ideally a CMOVAE instruction
// If the modular_subtraction function is called from within a loop, and b
// remains unchanged throughout the loop, then the compiler will ideally loop-
// hoist the "diff = n - b" outside of the loop, if we assume that the modulus
// n is also constant throughout the loop (which is extremely likely).
// The  a+diff  and  a-b  can run in parallel in the same cycle, and if the
// ternary operator is implemented via CMOVAE, then this means the
// modular_subtraction function would have 2 cycles total latency.  But this
// requires that the compiler is able to loop hoist the  diff = n-b,  so this is
// the ideal situation mentioned above in this note.  If the compiler can't loop
// hoist the calculation of diff (presumably because b does not remain constant
// throughout a loop), then it would have 3 cycles latency and would need 5 uops
// (including one uop to copy n to another register prior to the calculation of
// diff).  In contrast the default_modsub::call function has 3 cycles latency
// and needs 3 uops, but can not take advantage of potential loop hoisting to
// reduce the latency to 2.
// I do not believe it is worth the added complexity to offer the alternative
// version of modular_subtraction discussed in this note.  If someone is writing
// code that needs a low latency modular_subtraction within a loop where loop
// hoisting would be applicable, then an option usable right now is to do a
// modular negate of b outside the loop, and then use modular_addition inside
// the loop, because modular_addition *does* take advantage of the potential of
// loop hoisting to lower the total latency to 2 cycles.


// primary template
template <typename T>
struct impl_modular_subtraction {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    return default_modsub::call(a, b, modulus);
  }
};


// MSVC doesn't support inline asm so we skip it.
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MODSUB)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
template <>
struct impl_modular_subtraction<std::uint32_t> {
  HURCHALLA_FORCE_INLINE static
  std::uint32_t call(std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
  {
    using std::uint32_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    // Note: we want to make sure the LEA instruction doesn't use RBP/EBP or R13
    // for the base register, since that would necessitate a slower form of LEA
    // that has an extra 2 cycles latency and half the throughput of the fast
    // form.  We prevent this by using the "U" constraint, which allows  RAX,
    // RCX, RDX, RSI, RDI, R8, R9, R10, R11  for the System V AMD64 calling
    // convention, and  RAX, RCX, RDX, R8, R9, R10, R11  for the Microsoft x64
    // convention.  ICC (intel compiler) and clang don't support "U", so we use
    // "abcdSD" for them (allowing rax, rbx, rcx, rdx, rsi, rdi).
    uint32_t tmp = a;  // we prefer not to overwrite an input (a)
    uint32_t result;
    __asm__ ("subl %[b], %[tmp] \n\t"             /* tmp = a - b */
             "leal (%q[tmp], %q[m]), %[res] \n\t" /* res = tmp + modulus */
             "cmovael %[tmp], %[res] \n\t"        /* res = (a>=b) ? tmp : res */

# if defined(__INTEL_COMPILER)
             : [tmp]"+&abcdSD"(tmp), [res]"=r"(result)
# elif defined(__clang__)    /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
                             /* clang seems to use the first register listed. */
                             /* rcx is probably a good first choice. */
             : [tmp]"+&cabdSD"(tmp), [res]"=r"(result)
# else
             : [tmp]"+&UabcdSD"(tmp), [res]"=r"(result)
# endif

# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [m]"r"(modulus), [b]"r"(b)
# else
             : [m]"r"(modulus), [b]"rm"(b)
# endif

             : "cc");

    HPBC_POSTCONDITION2(result < modulus);  // uint32_t guarantees result>=0.
    HPBC_POSTCONDITION2(result == default_modsub::call(a, b, modulus));
    return result;
  }
};

template <>
struct impl_modular_subtraction<std::uint64_t> {
  HURCHALLA_FORCE_INLINE static
  std::uint64_t call(std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    // Note: the issues and solutions with LEA and RBP/EBP/R13 are the same here
    // as for the uint32_t version of this function above.
    uint64_t tmp = a;  // we prefer not to overwrite an input (a)
    uint64_t result;
    __asm__ ("subq %[b], %[tmp] \n\t"            /* tmp = a - b */
             "leaq (%[tmp], %[m]), %[res] \n\t"  /* res = tmp + modulus */
             "cmovaeq %[tmp], %[res] \n\t"       /* res = (a>=b) ? tmp : res */

# if defined(__INTEL_COMPILER)
             : [tmp]"+&abcdSD"(tmp), [res]"=r"(result)
# elif defined(__clang__)    /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
                             /* clang seems to use the first register listed. */
                             /* rcx is probably a good first choice. */
             : [tmp]"+&cabdSD"(tmp), [res]"=r"(result)
# else
             : [tmp]"+&UabcdSD"(tmp), [res]"=r"(result)
# endif

# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [m]"r"(modulus), [b]"r"(b)
# else
             : [m]"r"(modulus), [b]"rm"(b)
# endif

             : "cc");

    HPBC_POSTCONDITION2(result < modulus);  // uint64_t guarantees result>=0.
    HPBC_POSTCONDITION2(result == default_modsub::call(a, b, modulus));
    return result;
  }
};
#endif


}}  // end namespace

#endif
