// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_SUBTRACTION_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// With regard to the LowlatencyTag vs. the LowuopsTag versions here:
// If 'modulus' was not set/modified recently before the call of this modular
// subtraction function, AND if either 'a' or 'b' was not set/modified recently
// before the call, then the LowlatencyTag versions should usually provide lower
// latency than the LowuopsTag versions.
// Note that the LowlatencyTag versions will typically use more uops and create
// more pressure on the ALU than LowuopsTag, unless the compiler can loop hoist
// the extra instruction(s) that the LowlatencyTag versions use that involve 'a'
// (or 'b') and 'modulus'.

// Fyi: the purpose of having structs with static member functions is to
// disallow ADL and to make specializations simple and easy.

// primary template for default implementation
template <class PTAG>
struct default_modsub_unsigned {
};

// LowuopsTag, for low uops and low ALU use.
template <>
struct default_modsub_unsigned<LowuopsTag> {
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

    // These lines minimize uop count, register use, and ALU pressure.
    T diff = static_cast<T>(a - b);
    T result = static_cast<T>(diff + modulus);
      // result = (a >= b) ? diff : result
    result = ::hurchalla::conditional_select(a >= b, diff, result);

    HPBC_POSTCONDITION2(0<=result && result<modulus);
    return result;
  }
};

// LowlatencyTag.  If modulus and 'a'(or 'b') were not set/modified recently
// before the call, then the below function would typically have lowest possible
// latency.  (Specifically, diff = b - modulus  normally could either be loop
// hoisted by the compiler, or computed at the same time as earlier work by the
// CPU, costing zero latency.)
template <>
struct default_modsub_unsigned<LowlatencyTag> {
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

    // Naively, to achieve low latency, we might assume we need one version of
    // this function for when 'a' has remained mostly constant, and a different
    // version for when 'b' has remained constant.  It's fairly straightforward
    // to understand that when 'b' hasn't changed in the lines of code preceding
    // this function call, the compiler (due to function inlining) will be able
    // to schedule the line  T diff = b - modulus  either outside of a loop, or
    // at an earlier point where it can overlap with other CPU work at no extra
    // latency (via CPU pipelining or superscalar execution).
    // However, the important thing to point out is that in practice, this
    // function will probably also always work at achieving low latency when it
    // is 'a' that hasn't changed in the lines preceding the function call.
    // If/when the compiler sees that this is the case, it will almost always be
    // smart enough to transform the two lines below:
    // diff = b - modulus
    // tmp = a - diff
    // effectively (in its produced assembly) into these following two lines:
    // sum = a + modulus
    // tmp = sum - b
    // This transformation does not change the result of the function in any
    // way, but when the compiler sees that neither 'a' nor 'modulus' have
    // changed recently, the compiler now will be able to loop hoist 'sum' or
    // schedule it where it is sure to be calculated (by the CPU) in parallel
    // with earlier instructions, at a cost of zero additional latency.
    // This is why we only need one version of this function for low latency.
    // [Note: inside an inline asm section, we need to calculate  "a - b"  both
    // to have the result available and in order to set the carry flag for the
    // cmov that follows (note, on x86 the subtraction will overwrite 'a', or a
    // copy of it if the compiler made a copy before entering the asm section).
    // I don't think we can avoid doing this in the inline asm.  So if we were
    // for some reason to make two versions of the function, as best I can tell
    // we would want/need to use this same inline asm for both versions.]

    T diff = b - modulus;
    // note: the next two lines can begin on the same clock cycle
    T tmp = a - diff;
    T result = a - b;
      // result = (a < b) ? tmp : result;  //on x86, ideally a CMOVB instruction
    result = ::hurchalla::conditional_select(a < b, tmp, result);

    HPBC_POSTCONDITION2(0<=result && result<modulus);
    return result;
  }
};



// primary template
template <typename T, class PTAG>
struct impl_modular_subtraction_unsigned {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    return default_modsub_unsigned<PTAG>::call(a, b, modulus);
  }
};


// MSVC doesn't support inline asm so we skip it.
#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MODSUB)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)

template <>
struct impl_modular_subtraction_unsigned<std::uint32_t, LowuopsTag> {
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
    HPBC_POSTCONDITION2(result ==
                      default_modsub_unsigned<LowuopsTag>::call(a, b, modulus));
    return result;
  }
};

template <>
struct impl_modular_subtraction_unsigned<std::uint64_t, LowuopsTag> {
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
    HPBC_POSTCONDITION2(result ==
                      default_modsub_unsigned<LowuopsTag>::call(a, b, modulus));
    return result;
  }
};

#ifdef HURCHALLA_ENABLE_INLINE_ASM_128_BIT
template <>
struct impl_modular_subtraction_unsigned<__uint128_t, LowuopsTag> {
  HURCHALLA_FORCE_INLINE static
  __uint128_t call(__uint128_t a, __uint128_t b, __uint128_t modulus)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // __uint128_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // __uint128_t guarantees b>=0.

// We can't use LEA here, since our 128 bit operands would necessitate an add
// with carry to calculate a high 64 bit part, and LEA can neither produce nor
// consume a carry.  Therefore we'll implement this alternative in asm:
//  __uint128_t zero = 0;
//  __uint128_t diff = a - b;
//  __uint128_t modulus_or_zero = (a >= b) ? zero : modulus;
//  __uint128_t result = diff + modulus_or_zero;
// Note: since we aren't using LEA, we have no concern with RBP/EBP/R13
// like we did for the 32bit/64bit functions above.

    uint64_t reg = 0;

    uint64_t alo = static_cast<uint64_t>(a);
    uint64_t ahi = static_cast<uint64_t>(a >> 64);
    uint64_t blo = static_cast<uint64_t>(b);
    uint64_t bhi = static_cast<uint64_t>(b >> 64);
    uint64_t mlo = static_cast<uint64_t>(modulus);
    uint64_t mhi = static_cast<uint64_t>(modulus >> 64);
    __asm__ ("subq %[blo], %[alo] \n\t"         /* diff = a - b */
             "sbbq %[bhi], %[ahi] \n\t"
             "cmovaeq %[reg], %[mlo] \n\t"      /* mozlo = (a>=b) ? 0 : mlo */
             "cmovbq %[mhi], %[reg] \n\t"       /* mozhi = (a<b)  ? mhi : 0 */
             : [alo]"+&r"(alo), [ahi]"+&r"(ahi), [mlo]"+&r"(mlo), [reg]"+&r"(reg)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [blo]"r"(blo), [bhi]"r"(bhi), [mhi]"r"(mhi)
# else
             : [blo]"rm"(blo), [bhi]"rm"(bhi), [mhi]"rm"(mhi)
# endif
             : "cc");
    uint64_t difflo = alo;
    uint64_t diffhi = ahi;
    uint64_t mozlo = mlo;
    uint64_t mozhi = reg;
    __uint128_t diff = (static_cast<__uint128_t>(diffhi) << 64) | difflo;
    __uint128_t moz = (static_cast<__uint128_t>(mozhi) << 64) | mozlo;

    __uint128_t result = diff + moz;

    HPBC_POSTCONDITION2(result < modulus);  // __uint128_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
                      default_modsub_unsigned<LowuopsTag>::call(a, b, modulus));
    return result;
  }
};
#endif


// See explanation inside  default_modsub_unsigned<LowlatencyTag>  for why we
// don't need two different low latency functions for taking advantage of when
// 'b' was recently unchanged vs. when 'a' was recently unchanged.
template <>
struct impl_modular_subtraction_unsigned<std::uint32_t, LowlatencyTag> {
  HURCHALLA_FORCE_INLINE static
  std::uint32_t call(std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
  {
    using std::uint32_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    uint32_t diff = b - modulus;
    uint32_t tmp = a - diff;

    uint32_t a2 = a;   // we prefer not to overwrite an input
    // Note we don't use LEA, so we don't worry about RBP/EBP or R13
    __asm__ ("subl %[b], %[a2] \n\t"            /* res = a - b */
             "cmovbl %[tmp], %[a2] \n\t"        /* res = (a<b) ? tmp : res */
             : [a2]"+&r"(a2)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [b]"r"(b), [tmp]"r"(tmp)
# else
             : [b]"rm"(b), [tmp]"rm"(tmp)
# endif
             : "cc");
    uint32_t result = a2;

    HPBC_POSTCONDITION2(result < modulus);  // uint32_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
                   default_modsub_unsigned<LowlatencyTag>::call(a, b, modulus));
    return result;
  }
};

template <>
struct impl_modular_subtraction_unsigned<std::uint64_t, LowlatencyTag> {
  HURCHALLA_FORCE_INLINE static
  std::uint64_t call(std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    uint64_t diff = b - modulus;
    uint64_t tmp = a - diff;

    uint64_t a2 = a;   // we prefer not to overwrite an input
    // Note we don't use LEA, so we don't worry about RBP/EBP or R13
    __asm__ ("subq %[b], %[a2] \n\t"            /* res = a - b */
             "cmovbq %[tmp], %[a2] \n\t"        /* res = (a<b) ? tmp : res */
             : [a2]"+&r"(a2)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [b]"r"(b), [tmp]"r"(tmp)
# else
             : [b]"rm"(b), [tmp]"rm"(tmp)
# endif
             : "cc");
    uint64_t result = a2;

    HPBC_POSTCONDITION2(result < modulus);  // uint64_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
                   default_modsub_unsigned<LowlatencyTag>::call(a, b, modulus));
    return result;
  }
};

#ifdef HURCHALLA_ENABLE_INLINE_ASM_128_BIT
template <>
struct impl_modular_subtraction_unsigned<__uint128_t, LowlatencyTag> {
  HURCHALLA_FORCE_INLINE static
  __uint128_t call(__uint128_t a, __uint128_t b, __uint128_t modulus)
  {
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // __uint128_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // __uint128_t guarantees b>=0.

    __uint128_t diff = b - modulus;
    __uint128_t tmp = a - diff;

    uint64_t tmplo = static_cast<uint64_t>(tmp);
    uint64_t tmphi = static_cast<uint64_t>(tmp >> 64);
    uint64_t alo = static_cast<uint64_t>(a);
    uint64_t ahi = static_cast<uint64_t>(a >> 64);
    uint64_t blo = static_cast<uint64_t>(b);
    uint64_t bhi = static_cast<uint64_t>(b >> 64);
    // Note we don't use LEA, so we don't worry about RBP/EBP or R13
    __asm__ ("subq %[blo], %[alo] \n\t"         /* res = a - b */
             "sbbq %[bhi], %[ahi] \n\t"
             "cmovbq %[tmplo], %[alo] \n\t"     /* res = (a<b) ? tmp : res */
             "cmovbq %[tmphi], %[ahi] \n\t"
             : [alo]"+&r"(alo), [ahi]"+&r"(ahi)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [blo]"r"(blo), [bhi]"r"(bhi), [tmplo]"r"(tmplo), [tmphi]"r"(tmphi)
# else
             : [blo]"rm"(blo), [bhi]"rm"(bhi), [tmplo]"rm"(tmplo), [tmphi]"rm"(tmphi)
# endif
             : "cc");
    __uint128_t result = (static_cast<__uint128_t>(ahi) << 64) | alo;

    HPBC_POSTCONDITION2(result < modulus);  // __uint128_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
                   default_modsub_unsigned<LowlatencyTag>::call(a, b, modulus));
    return result;
  }
};
#endif

// end of inline asm functions for x86_64
#endif



// version for unsigned T
template <typename T, class PTAG, bool = ut_numeric_limits<T>::is_signed>
struct impl_modular_subtraction {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    return impl_modular_subtraction_unsigned<T,PTAG>::call(a, b, modulus);
  }
};

// version for signed T
template <typename T, class PTAG>
struct impl_modular_subtraction<T, PTAG, true> {
  HURCHALLA_FORCE_INLINE static T call(T a, T b, T modulus)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::is_signed, "");
    static_assert(static_cast<T>(-1) == ~(static_cast<T>(0)),
                                  "T must use two's complement representation");
    using U = typename extensible_make_unsigned<T>::type;
    static_assert(static_cast<T>(static_cast<U>(static_cast<T>(-1))) ==
                  static_cast<T>(-1), "Casting a signed T value to unsigned and"
                               " back again must result in the original value");
    HPBC_PRECONDITION2(modulus > 0);
    HPBC_PRECONDITION2(0 <= a && a < modulus);
    HPBC_PRECONDITION2(0 <= b && b < modulus);

#if defined(HURCHALLA_AVOID_CSELECT)
    static_assert((static_cast<T>(-1) >> 1) == static_cast<T>(-1),
                          "Arithmetic right shift is required but unavailable");
    T tmp = static_cast<T>(a - b);
    // if tmp is negative, use a bit mask of all 1s.  Otherwise use all 0s.
    U mask = static_cast<U>(tmp >> ut_numeric_limits<T>::digits);
    U masked_modulus = static_cast<U>(mask & static_cast<U>(modulus));
    U result = static_cast<U>(static_cast<U>(tmp) + masked_modulus);
    HPBC_ASSERT2(result == (impl_modular_subtraction_unsigned<U,PTAG>::call(
               static_cast<U>(a), static_cast<U>(b), static_cast<U>(modulus))));
#else
    U result= impl_modular_subtraction_unsigned<U,PTAG>::call(static_cast<U>(a),
                                    static_cast<U>(b), static_cast<U>(modulus));
#endif

    return static_cast<T>(result);
  }
};


}}  // end namespace

#endif
