// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_ADDITION_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_IMPL_MODULAR_ADDITION_H_INCLUDED


#include "hurchalla/modular_arithmetic/traits/extensible_make_unsigned.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


// See version #2 of the template function below for corresponding C++ code.
#if defined(HURCHALLA_TARGET_ISA_X86_64) && !defined(_MSC_VER) && \
    !defined(__clang__) // msvc doesn't support inline asm, and clang compiles
                        // template #2 to optimal asm (with more flexibility)
inline uint32_t impl_modular_addition_prereduced_inputs(uint32_t a, uint32_t b,
                                                               uint32_t modulus)
{
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint32_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint32_t guarantees b>=0.

    // [input operands]:  anon_reg0 = a, anon_reg1 = b, anon_reg2 = modulus
    // [output operands]: result = anon_reg0
    // [clobber list]:    both add and sub clobber the FLAGS register ["cc"]
    uint32_t result, sum;
    __asm__ ("leal (%q2,%q3), %1 \n\t"
             "subl %4, %0 \n\t"
             "addl %3, %0 \n\t"
             "cmovael %1, %0 \n\t"
             : "=&r"(result), "=&r"(sum)
             : "%0"(a), "r"(b), "r"(modulus)
             : "cc");
    return result;
}
// See version #2 of the template function below for corresponding C++ code.
inline uint64_t impl_modular_addition_prereduced_inputs(uint64_t a, uint64_t b,
                                                               uint64_t modulus)
{
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(a<modulus);  // uint64_t guarantees a>=0.
    HPBC_PRECONDITION2(b<modulus);  // uint64_t guarantees b>=0.

    // [input operands]:  anon_reg0 = a, anon_reg1 = b, anon_reg2 = modulus
    // [output operands]: result = anon_reg0
    // [clobber list]:    both add and sub clobber the FLAGS register ["cc"]
    uint64_t result, sum;
    __asm__ ("leaq (%2,%3), %1 \n\t"
             "subq %4, %0 \n\t"
             "addq %3, %0 \n\t"
             "cmovaeq %1, %0 \n\t"
             : "=&r"(result), "=&r"(sum)
             : "%0"(a), "r"(b), "r"(modulus)
             : "cc");
    return result;
}
#endif


// I expect the following ARM64 asm code is correct, but I'm unable to test it
// for now.  Therefore I've disabled the arm asm until it's tested.  In any
// event, clang produces optimal asm from template #2 (but gcc does not).
#if 0
// See version #2 of the template function below for corresponding C++ code.
#if !defined(_MSC_VER) && (defined(HURCHALLA_TARGET_ISA_ARM_64)
inline uint64_t impl_modular_addition_prereduced_inputs(uint64_t a, uint64_t b,
                                                               uint64_t modulus)
{
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(0<=a && a<modulus);  // i.e. the input must be prereduced
    HPBC_PRECONDITION2(0<=b && b<modulus);  // i.e. the input must be prereduced

    // [input operands]:  anon_reg0 = a, anon_reg1 = b, anon_reg2 = modulus
    // [output operands]: result = anon_reg0
    // [clobber list]:    both add and sub clobber the FLAGS register ["cc"]
    uint64_t result, sum;
    __asm__ ("add %2, %3, %1 \n\t"
             "sub %4, %2, %0 \n\t"
             "adds %3, %0, %0 \n\t"
             "csel %1, %0, %0, hs \n\t"
             : "=&r"(result), "=&r"(sum)
             : "%0"(a), "r"(b), "r"(modulus)
             : "cc");
    return result;
}
#endif
#endif  // #if 0




// -----------------------------------------------------------------------------
// Template function versions.  By C++ rules, these have second priority for
// function argument matching - the non-template versions above have first
// priority.
// -----------------------------------------------------------------------------


// --- Version #1 ---
// This is an easy to understand version of the templated function, which would
// in actuality work correctly for any integer type.  However, we restrict it
// with enable_if, so that the fastest template version will be used for the
// appropriate type T.
// The enable_if condition of !(std::is_integral<T>::value) ensures that this
// function is enabled only for the integral types that std::type_traits doesn't
// recognize - i.e. non-native integral types, such as for example
// boost::multiprecision::int256_t.
template <typename T>
#ifdef _MSC_VER  // MSVC(as of 06/20) produces quite sub-optimal asm for Version
    T            // #2.  Therefore we want MSVC to always use this Version #1.
#else
    typename std::enable_if< !(std::is_integral<T>::value) &&
                         (ma_numeric_limits<T>::is_signed), T >::type
#endif
impl_modular_addition_prereduced_inputs(T a, T b, T modulus)
{
    static_assert(ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION2(modulus>0);
    HPBC_PRECONDITION2(0<=a && a<modulus);  // i.e. the input must be prereduced
    HPBC_PRECONDITION2(0<=b && b<modulus);  // i.e. the input must be prereduced

    // We want essentially-  result = (a+b < modulus) ? a+b : a+b-modulus
    //   But due to potential overflow on a+b we need to write it as follows
    T tmp = static_cast<T>(modulus - b);
    T result = (a < tmp) ? static_cast<T>(a+b) : static_cast<T>(a - tmp);
    return result;
}


// --- Version #2 ---
// This is a somewhat difficult to understand version of the template function.
// It can in general compile to assembly that is a cycle or two faster than the
// easy-to-understand template Version #1.  However, we don't always enable this
// template (via std::enable_if), because it uses casts between signed and
// unsigned; these casts potentially might not be free for a non-native signed
// type.  The proof of this function's correctness is given by the theorem in
// the comments just below.  See the notes at the end of those comments to
// understand how this function implements the theorem.
#ifndef _MSC_VER   // MSVC (as of 06/20) produces sub-optimal asm from this
    template <typename T>
    typename std::enable_if< (std::is_integral<T>::value) ||
                             !(ma_numeric_limits<T>::is_signed), T >::type
    impl_modular_addition_prereduced_inputs(T a, T b, T modulus)
    {
        static_assert(ma_numeric_limits<T>::is_integer, "");
        HPBC_PRECONDITION2(modulus>0);
        HPBC_PRECONDITION2(0<=a && a<modulus);  // the input must be prereduced
        HPBC_PRECONDITION2(0<=b && b<modulus);  // the input must be prereduced
    
        using U = typename extensible_make_unsigned<T>::type;
        U sum = static_cast<U>(static_cast<U>(a) + static_cast<U>(b));
        U tmp = static_cast<U>(static_cast<U>(a) - static_cast<U>(modulus));
        U tmp2 = static_cast<U>(static_cast<U>(b) + tmp);
        U result = (tmp2 < static_cast<U>(b)) ? tmp2 : sum;
        return static_cast<T>(result);
    }
#endif

// ---------Theorem and proof--------
// The constant "R" used below represents the value
// R = 2^(std::numeric_limits<T>::digits).  For example, if T is uint64_t,
// then R = 2^64.
// We'll use a psuedo-cast notation of (integer)x to indicate when we are
// treating x as an infinite precision signed integer - i.e. a member of the set
// Z of mathematical integers.
// The notation "%%" used below should be interpreted as a conceptual modulo
// operator that will always produce a non-negative remainder for the result.
// This is slightly different from the actual C/C++ modulo operator "%", which
// produces a negative remainder if the dividend is negative.

// Theorem:  Require 0 <= a < m, 0 <= b < m, and 0 < m < R; and let
// tmp = ((integer)a - (integer)m) %% R.
//    If ((integer)b + (integer)tmp >= R), then
//        ((integer)a + (integer)b) %% m == b+a-m.
//    else
//        ((integer)a + (integer)b) %% m == a+b.
//
// Proof:
// As a precondition, we know  0 <= a < m.
// [1] Therefore  -(integer)m <= (integer)a - (integer)m < 0
// As a precondition, we know  m < R,  thus
// [2]  -R < -(integer)m
// Combining [1] and [2],  -R < (integer)a - (integer)m < 0.
// Adding R to all parts,  0 < (integer)a - (integer)m + R < R.
// We can see this expression is bound between 0 and R, so
// (integer)a - (integer)m + R == ((integer)a - (integer)m + R) %% R
// [3] Thus  ((integer)a - (integer)m) %% R == (integer)a - (integer)m + R.
//
// C/C++ (and assembly) performs unsigned addition and subtraction modulo R, and
// of course produces an unsigned non-negative result.
// Therefore,  tmp = a - m == ((integer)a - (integer)m) %% R.  Thus
// 0 <= tmp < R.   And by [3],
// [4]  tmp = a - m == (integer)a - (integer)m + R.
//
// We want to test whether  (integer)a + (integer)b >= (integer)m,  because
// if this test is true, then  ((integer)a + (integer)b) %% m == a+b-m.  If it
// is false, then  ((integer)a + (integer)b) %% m == a+b.
// We can re-express this test as  (integer)a - (integer)m + R >= R - (integer)b
// Combining this with [4], the test becomes  (integer)tmp >= R - (integer)b,
// which we can rearrange into the test
// [5]  (integer)b + (integer)tmp >= R.
// The original form of the test can't be directly evaluated in C/C++ or in
// assembly, but [5] can be easily evaluted by performing the addition  b + tmp,
// and detecting whether or not the addition overflows.
//
// [6] If b+tmp overflows, then  (integer)b + (integer)a >= (integer)m,  and
//    thus  ((integer)a + (integer)b) %% m == a+b-m.
//    We know this because we showed above that testing for overflow is
//    equivalent to testing whether  (integer)b + (integer)tmp >= R,  which is
//    equivalent to testing whether  (integer)a + (integer)b >= (integer)m.
// [7] If b+tmp does not overflow, then  (integer)b + (integer)a < (integer)m,
//    and thus  ((integer)a + (integer)b) %% m == a+b.

// Notes:
// To detect overflow on the add (b+tmp):
//    In assembly, we can perform the add and then look at the carry flag to see
//    if it overflowed.  In C, we can test if (b+tmp < b).  If true, then the
//    add overflowed.  Note: with luck the compiler will recognize this C idiom
//    and produce assembly that is the same as what we would write.  Clang
//    typically does fairly well at this, gcc less well, and icc least well.
//
// Putting together [6] and [7], in C++ we could write
// modular_addition_prereduced(unsigned T a, T b, T m) {
//    static_assert( !(std::numeric_limits<T>::is_signed), "" );
//    T tmp = a-m;
//    return (b+tmp < b) ? b+tmp : a+b;
// }


}}  // end namespace

#endif
