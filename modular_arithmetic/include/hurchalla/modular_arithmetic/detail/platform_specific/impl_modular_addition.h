// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace detail {


// -----------------------------------------------------------------------------
// Template function versions.  By C++ rules, these have second priority for
// function argument matching - the non-template versions below it have first
// priority.
// -----------------------------------------------------------------------------


#if 0
// --- Version #1 ---
// This is an easy to understand version of the templated function.
//
// However, on x86, Version #1's use of (modulus - b) will very often require
// two uops because the compiler will see 'modulus' stays constant over a long
// period of time, and hence it will choose not to overwrite its register -
// necessitating an extra uop for a copy of the register.  This fact may be
// inconsequential, for example if the compiler is able to hoist both the copy 
// and the subtraction out of a loop, but if 'b' does not stay constant within
// a loop, then hoisting is not possible and Version #1 would require 2 uops
// for the subtraction in a perhaps critical loop.
// In contrast, Version #2 uses (b - modulus) which requires only one uop when
// 'b' does not stay constant within a loop.  This fact is why we slightly
// prefer Version #2, and why we have disabled Version #1 and enabled Version #2
// via #if.  Note that if 'b' stays constant for a period of time, or if we are
// compiling for ARM, we can expect total uops to be the same between the two
// function versions.  In all situations, we can expect the latency of the two
// versions to be the same.

template <typename T> HURCHALLA_FORCE_INLINE
T impl_modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // i.e. the input must be prereduced
    HPBC_PRECONDITION2(b<modulus);  // i.e. the input must be prereduced

    // We want essentially-  result = (a+b < modulus) ? a+b : a+b-modulus
    //   But due to the potential for overflow on a+b, we need to instead test
    //   the alternative predicate (a < modulus-b), which gives us our desired
    //   result without any problem of overflow.  So we can and should use:
    //   result = (a < modulus-b) ? a+b : a+b-modulus
    T tmp = static_cast<T>(modulus - b);
    T sum = static_cast<T>(a + b);
    T tmp2 = static_cast<T>(a - tmp);
    T result = (a < tmp) ? sum : tmp2;

    HPBC_POSTCONDITION2(static_cast<T>(0) <= result && result < modulus);
    return result;
}
#else
// --- Version #2 ---
// This is a more difficult to understand version of the template function.  The
// proof of this function's correctness is given by the theorem in the comments
// at the end of this file.  See the notes at the end of those comments to
// understand the implementation details.
template <typename T> HURCHALLA_FORCE_INLINE
T impl_modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // the input must be prereduced
    HPBC_PRECONDITION2(b<modulus);  // the input must be prereduced

    using U = typename extensible_make_unsigned<T>::type;
    U sum = static_cast<U>(static_cast<U>(a) + static_cast<U>(b));
    U tmp = static_cast<U>(static_cast<U>(b) - static_cast<U>(modulus));
    U result = static_cast<U>(static_cast<U>(a) + tmp);
# if 0
    result = (result >= static_cast<U>(a)) ? sum : result;
# else
    HURCHALLA_CMOV(result >= static_cast<U>(a), result, sum);
# endif

    HPBC_POSTCONDITION2(static_cast<U>(0) <= result &&
                                              result < static_cast<U>(modulus));
    return static_cast<T>(result);
}
#endif



// -----------------------------------------------------------------------------
// These non-template functions get first priority for argument matching.
// -----------------------------------------------------------------------------

// These inline asm functions implement version #2 of the template function.
// MSVC doesn't support inline asm, so we skip it.

#if (defined(HURCHALLA_ALLOW_INLINE_ASM_ALL) || \
     defined(HURCHALLA_ALLOW_INLINE_ASM_MODADD)) && \
    defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER)
HURCHALLA_FORCE_INLINE std::uint32_t impl_modular_addition_prereduced_inputs(
                        std::uint32_t a, std::uint32_t b, std::uint32_t modulus)
{
    using std::uint32_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    // By calculating tmp outside of the __asm__, we allow the compiler to
    // potentially loop hoist tmp, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint32_t sum = a + b;
    uint32_t tmp = b - modulus;
    uint32_t tmp2 = a;  // we prefer not to overwrite an input (a)
    __asm__ ("addl %[tmp], %[tmp2] \n\t"       /* tmp2 = a + tmp */
             "cmovbl %[tmp2], %[sum] \n\t"     /* sum = (tmp2<a) ? tmp2 : sum */
             : [tmp2]"+&r"(tmp2), [sum]"+r"(sum)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [tmp]"r"(tmp)
# else
             : [tmp]"rm"(tmp)
# endif
             : "cc");
    uint32_t result = sum;

    HPBC_POSTCONDITION2(result < modulus);  // uint32_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
              impl_modular_addition_prereduced_inputs<uint32_t>(a, b, modulus));
    return result;
}

HURCHALLA_FORCE_INLINE std::uint64_t impl_modular_addition_prereduced_inputs(
                        std::uint64_t a, std::uint64_t b, std::uint64_t modulus)
{
    using std::uint64_t;
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    // By calculating tmp outside of the __asm__, we allow the compiler to
    // potentially loop hoist tmp, if this function is inlined into a loop.
    // https://en.wikipedia.org/wiki/Loop-invariant_code_motion
    uint64_t sum = a + b;
    uint64_t tmp = b - modulus;
    uint64_t tmp2 = a;  // we prefer not to overwrite an input (a)
    __asm__ ("addq %[tmp], %[tmp2] \n\t"       /* tmp2 = a + tmp */
             "cmovbq %[tmp2], %[sum] \n\t"     /* sum = (tmp2<a) ? tmp2 : sum */
             : [tmp2]"+&r"(tmp2), [sum]"+r"(sum)
# if defined(__clang__)        /* https://bugs.llvm.org/show_bug.cgi?id=20197 */
             : [tmp]"r"(tmp)
# else
             : [tmp]"rm"(tmp)
# endif
             : "cc");
    uint64_t result = sum;

    HPBC_POSTCONDITION2(result < modulus);  // uint64_t guarantees result>=0.
    HPBC_POSTCONDITION2(result ==
              impl_modular_addition_prereduced_inputs<uint64_t>(a, b, modulus));
    return result;
}
#endif





// ---------Theorem and proof for Version #2--------
// The constant "R" used below represents the value
// R = 2^(ut_numeric_limits<T>::digits).  For example, if T is uint64_t,
// then R = 2^64.
// We'll use a psuedo-cast notation of (Z)x to indicate when we are treating x
// as an infinite precision signed integer - i.e. a member of the set Z of
// mathematical integers.
// The notation "%%" used below should be interpreted as a conceptual modulo
// operator that will always produce a non-negative remainder for the result.
// This is slightly different from the actual C/C++ modulo operator "%", which
// produces a negative remainder if the dividend is negative.

// Theorem:  Require 0 <= a < m, 0 <= b < m, and 0 < m < R; and let
// tmp = ((Z)b - (Z)m) %% R.
//    If ((Z)a + (Z)tmp >= R), then
//        ((Z)a + (Z)b) %% m == a+b-m.
//    else
//        ((Z)a + (Z)b) %% m == a+b.
//
// Proof:
// As a precondition, we know  0 <= b < m.
// [1] Therefore  -(Z)m <= (Z)b - (Z)m < 0
// As a precondition, we know  m < R,  thus
// [2]  -R < -(Z)m
// Combining [1] and [2],  -R < (Z)b - (Z)m < 0.
// Adding R to all parts,  0 < (Z)b - (Z)m + R < R.
// We can see this expression is bound between 0 and R, so
// (Z)b - (Z)m + R == ((Z)b - (Z)m + R) %% R.
// [3] Thus  ((Z)b - (Z)m) %% R == (Z)b - (Z)m + R.
//
// C/C++ (and assembly) performs unsigned addition and subtraction modulo R, and
// of course produces an unsigned non-negative result.
// Therefore, tmp = b - m == ((Z)b - (Z)m) %% R.  Thus 0 <= tmp < R, and by [3],
// [4]  tmp = b - m == (Z)b - (Z)m + R.
//
// We would like to test whether  (Z)a + (Z)b >= (Z)m,  but this test can't be
// directly evaluated in C/C++/assembly due to potential overflow on  a+b.
// However, we can re-express this test as  (Z)b - (Z)m + R >= R - (Z)a,  and
// combining this with [4], this test becomes  (Z)tmp >= R - (Z)a.  We can then
// rearrange this last test into the test  (Z)a + (Z)tmp >= R.
// These tests are all equivalent, and so
// [5]  (Z)a + (Z)tmp >= R  implies  (Z)a + (Z)b >= (Z)m.  And
// [6]  (Z)a + (Z)tmp < R   implies  (Z)a + (Z)b < (Z)m.
// Note that we can easily evaluate the test  (Z)a + (Z)tmp >= R  using
// C/C++/assembly, by performing the addition  a + tmp,  and detecting whether
// or not the add overflows.
//
// As preconditions, we know  0 <= a < m,  0 <= b < m,  and  m < R.
// [7]  Therefore  0 <= (Z)a + (Z)b < (Z)m + (Z)m.
// [8]  Assume  (Z)a + (Z)b >= (Z)m:
//    Then by [7]          (Z)m <= (Z)a + (Z)b < (Z)m + (Z)m,  and
//                            0 <= (Z)a + (Z)b - (Z)m < (Z)m.  Thus,
//    ((Z)a + (Z)b - (Z)m) %% m == (Z)a + (Z)b - (Z)m,  and thus
//           ((Z)a + (Z)b) %% m == (Z)a + (Z)b - (Z)m.
//    Since (Z)m < R,         0 <= (Z)a + (Z)b - (Z)m < R.  Thus,
//    ((Z)a + (Z)b - (Z)m) %% R == (Z)a + (Z)b - (Z)m,  and therefore in C/C++,
//                    a + b - m == (Z)a + (Z)b - (Z)m.  This gives us
//           ((Z)a + (Z)b) %% m == a + b - m.
// [9]  Assume  (Z)a + (Z)b < (Z)m:
//    Then by [7]             0 <= (Z)a + (Z)b < (Z)m.  Thus,
//           ((Z)a + (Z)b) %% m == (Z)a + (Z)b
//    Since (Z)m < R,         0 <= (Z)a + (Z)b < R.  Thus,
//           ((Z)a + (Z)b) %% R == (Z)a + (Z)b,  and therefore in C/C++,
//                        a + b == (Z)a + (Z)b.  This gives us
//           ((Z)a + (Z)b) %% m == a + b.
//
// [10]  Combining [5] with [8],  if ((Z)a + (Z)tmp >= R)  then
//                               ((Z)a + (Z)b) %% m == a+b-m.
// [11]  Combining [6] with [9],  if ((Z)a + (Z)tmp < R)  then
//                               ((Z)a + (Z)b) %% m == a+b.
//
// Implementation notes:
// As stated above, in C/C++/assembly, we can test if  (Z)a + (Z)tmp >= R  by
//    detecting if the addition a + tmp  overflows.
// To detect overflow on the addition (a+tmp):
//    In assembly, we can perform the add and then look at the carry flag to see
//    if it overflowed.  In C, we can test if (a+tmp >= a).  If false, then the
//    add overflowed.  Note: with luck the compiler will recognize this C idiom
//    and produce assembly that is the same as what we would write.  Clang
//    seems to do fairly well at this, gcc less well, and icc least well.
//
// Putting together [4], [10], and [11], in C++ we could write
// modular_addition_prereduced(unsigned T a, T b, T m) {
//    static_assert( !(ut_numeric_limits<T>::is_signed), "" );
//    T tmp = b-m;
//    return (a+tmp >= a) ? a+b : a+tmp;
// }


}}  // end namespace

#endif
