// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontgomeryDefault.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/montgomery_pow.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/util/traits/is_equality_comparable.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>
#include <array>
#include <cstddef>

namespace hurchalla {


// T must be a signed or unsigned integral type.
// Usually you would not specify MontyType and instead accept the default.
template<class T, class MontyType = typename detail::MontgomeryDefault<T>::type>
class MontgomeryForm final {
    const MontyType impl;
    using U = typename MontyType::uint_type;
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!(ut_numeric_limits<U>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<U>::digits >=
                  ut_numeric_limits<T>::digits, "");
public:
    using IntegerType = T;

    // If you need to compare MontgomeryValues for equality or inequality, call
    // getCanonicalValue() and compare the resulting CanonicalValues.
    using MontgomeryValue = typename MontyType::montvalue_type;
    static_assert(std::is_default_constructible<MontgomeryValue>::value, "");
    static_assert(std::is_copy_constructible<MontgomeryValue>::value, "");

    // A CanonicalValue is the unique representation of any value in Montgomery
    // form; a value of type MontgomeryValue may or may not be a unique
    // representation.  Since it represents a unique value, you can compare
    // CanonicalValues for equality and inequality.  Use getCanonicalValue() to
    // get a CanonicalValue from a MontgomeryValue.
    // CanonicalValue is implicitly convertible to MontgomeryValue; thus you can
    // use a CanonicalValue anywhere that requires a MontgomeryValue.
    using CanonicalValue = typename MontyType::canonvalue_type;
    static_assert(std::is_default_constructible<CanonicalValue>::value, "");
    static_assert(std::is_copy_constructible<CanonicalValue>::value, "");
    static_assert(is_equality_comparable<CanonicalValue>::value, "");
    static_assert(
               std::is_convertible<CanonicalValue, MontgomeryValue>::value, "");

    // A FusingValue is used for the addend or subtrahend for fused-multiply-
    // add (fmadd) or fused-multiply-subtract (fmsub).
    using FusingValue = typename MontyType::fusingvalue_type;
    static_assert(std::is_default_constructible<FusingValue>::value, "");
    static_assert(std::is_copy_constructible<FusingValue>::value, "");
    static_assert(std::is_convertible<FusingValue, MontgomeryValue>::value, "");


    explicit MontgomeryForm(T modulus) : impl(static_cast<U>(modulus))
    {
        // Usually odd modulus is required, but since MontyWrappedStandardMath
        // is an exception to the rule, we do not require odd modulus here.
        // We leave it to the constructor of the MontyType member (impl) to
        // enforce the requirement of an odd modulus, if applicable.
        //HPBC_PRECONDITION(modulus % 2 == 1);
        HPBC_PRECONDITION(modulus > 1);
    }

    // Returns the largest valid modulus allowed for the constructor.
    static constexpr T max_modulus()
    {
        return (MontyType::max_modulus() >
                                    static_cast<U>(ut_numeric_limits<T>::max()))
                ? ((ut_numeric_limits<T>::max() % 2 == 0)
                    ? static_cast<T>(ut_numeric_limits<T>::max() - 1)
                    : ut_numeric_limits<T>::max())
                : static_cast<T>(MontyType::max_modulus());
    }

    // Returns the modulus given to the constructor
    HURCHALLA_FORCE_INLINE
    T getModulus() const { return static_cast<T>(impl.getModulus()); }

    // Returns the converted value of the standard number 'a' into monty form.
    // Requires a >= 0.  (Note there is no restriction on how large 'a' can be.)
    HURCHALLA_FORCE_INLINE
    MontgomeryValue convertIn(T a) const
    {
        HPBC_PRECONDITION(a >= 0);
        return impl.convertIn(static_cast<U>(a));
    }

    // Converts (montgomery value) x into a "normal" number; returns the result.
    // Guarantees 0 <= result < modulus.
    HURCHALLA_FORCE_INLINE
    T convertOut(MontgomeryValue x) const
    {
        T a = static_cast<T>(impl.convertOut(x));
        HPBC_POSTCONDITION(0 <= a && a < getModulus());
        return a;
    }

    // Returns a unique (canonical) value representing the equivalence class of
    // x modulo the modulus.  You can not directly compare MontgomeryValues, but
    // you can call getCanonicalValue(), and then use standard equality or
    // inequality operators to compare the resulting CanonicalValues.
    HURCHALLA_FORCE_INLINE
    CanonicalValue getCanonicalValue(MontgomeryValue x) const
    {
        return impl.getCanonicalValue(x);
    }
    // Returns a value representing the equivalence class of x modulo the
    // modulus, which can be used as an argument for the fused-multiply-add/sub
    // functions "fmadd" and "fmsub".
    HURCHALLA_FORCE_INLINE
    FusingValue getFusingValue(MontgomeryValue x) const
    {
        return impl.getFusingValue(x);
    }

    // Returns the canonical monty value that represents the type T value 1.
    // The call is equivalent to getCanonicalValue(convertIn(static_cast<T>(1)))
    // but it's more efficient (essentially zero cost) and more convenient.
    HURCHALLA_FORCE_INLINE
    CanonicalValue getUnityValue() const
    {
        return impl.getUnityValue();
    }
    // Returns the canonical monty value that represents the type T value 0.
    // The call is equivalent to getCanonicalValue(convertIn(static_cast<T>(0)))
    // but it's more efficient (essentially zero cost) and more convenient.
    HURCHALLA_FORCE_INLINE
    CanonicalValue getZeroValue() const
    {
        return impl.getZeroValue();
    }
    // Returns the canonical monty value that represents the type T value
    // modulus-1 (which equals -1 (mod modulus)).  The call is equivalent to
    // getCanonicalValue(convertIn(static_cast<T>(modulus - 1))), but it's more
    // efficient (essentially zero cost) and more convenient.
    HURCHALLA_FORCE_INLINE
    CanonicalValue getNegativeOneValue() const
    {
        return impl.getNegativeOneValue();
    }


    // Returns the modular sum of (the montgomery values) x and y.  Performance
    // note: add may have lower latency than subtract (it should never have
    // higher latency), and subtract may use fewer uops than add (it should
    // never use more uops).
    HURCHALLA_FORCE_INLINE
    MontgomeryValue add(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.add(x, y);
    }
    // Returns the modular sum of the montgomery value x and the canonical
    // value y.  Performance note: this function is sometimes more efficient
    // than the above add() function and should never be less efficient, so it
    // can be useful to call getCanonicalValue outside of a loop to get 'y' as
    // a CanonicalValue, when you are calling add() inside a loop.
    HURCHALLA_FORCE_INLINE
    MontgomeryValue add(MontgomeryValue x, CanonicalValue y) const
    {
        MontgomeryValue ret = impl.add(x, y);
        HPBC_ASSERT(getCanonicalValue(ret) ==
                    getCanonicalValue(add(x, static_cast<MontgomeryValue>(y))));
        return ret;
    }
    HURCHALLA_FORCE_INLINE
    MontgomeryValue add(CanonicalValue x, MontgomeryValue y) const
    {
        return add(y, x);
    }
    // Adding CanonicalValues can be more efficient than the above functions
    HURCHALLA_FORCE_INLINE
    CanonicalValue add(CanonicalValue x, CanonicalValue y) const
    {
        CanonicalValue ret = impl.add(x, y);
        HPBC_ASSERT(ret == getCanonicalValue(add(
            static_cast<MontgomeryValue>(x), static_cast<MontgomeryValue>(y))));
        return ret;
    }


    // Returns the modular difference of (the montgomery values) x and y.  More
    // precisely, x minus y.
    HURCHALLA_FORCE_INLINE
    MontgomeryValue subtract(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.subtract(x, y);
    }
    // Returns the modular difference of the montgomery value x minus the
    // canonical value y.  Performance note: this function is sometimes more
    // efficient than the above subtract() function and should never be less
    // efficient, so it can be useful to call getCanonicalValue outside of a
    // loop to get 'y' as a CanonicalValue, when you are calling subtract()
    // inside a loop.
    HURCHALLA_FORCE_INLINE
    MontgomeryValue subtract(MontgomeryValue x, CanonicalValue y) const
    {
        MontgomeryValue ret = impl.subtract(x, y);
        HPBC_ASSERT(getCanonicalValue(ret) ==
               getCanonicalValue(subtract(x, static_cast<MontgomeryValue>(y))));
        return ret;
    }
    HURCHALLA_FORCE_INLINE
    MontgomeryValue subtract(CanonicalValue x, MontgomeryValue y) const
    {
        MontgomeryValue ret = impl.subtract(x, y);
        HPBC_ASSERT(getCanonicalValue(ret) ==
               getCanonicalValue(subtract(static_cast<MontgomeryValue>(x), y)));
        return ret;
    }
    // Subtracting CanonicalValues can be more efficient than the above funcs
    HURCHALLA_FORCE_INLINE
    CanonicalValue subtract(CanonicalValue x, CanonicalValue y) const
    {
        CanonicalValue ret = impl.subtract(x, y);
        HPBC_ASSERT(ret == getCanonicalValue(subtract(
            static_cast<MontgomeryValue>(x), static_cast<MontgomeryValue>(y))));
        return ret;
    }

    // Returns either the modular subtraction x-y, or y-x.  It is unspecified
    // which of the two subtractions will be returned.  For situations where you
    // don't care which subtraction gets performed, unorderedSubtract() will
    // usually perform slightly better than this->subtract() [in terms of number
    // of instructions, number of registers used, and possibly total latency].
    HURCHALLA_FORCE_INLINE
    MontgomeryValue unorderedSubtract(MontgomeryValue x,
                                       MontgomeryValue y) const
    {
        return impl.unordered_subtract(x, y);
    }
    HURCHALLA_FORCE_INLINE
    MontgomeryValue unorderedSubtract(MontgomeryValue x,
                                       CanonicalValue y) const
    {
        return impl.unordered_subtract(x, y);
    }
    HURCHALLA_FORCE_INLINE
    MontgomeryValue unorderedSubtract(CanonicalValue x,
                                       MontgomeryValue y) const
    {
        return impl.unordered_subtract(x, y);
    }

    // Returns the modular negation of the montgomery value x.
    HURCHALLA_FORCE_INLINE
    MontgomeryValue negate(MontgomeryValue x) const
    {
        return impl.negate(x);
    }
    // Returns the modular negation of the canonical value x.
    HURCHALLA_FORCE_INLINE
    CanonicalValue negate(CanonicalValue x) const
    {
        // we dont expect an impl would be able to do better than this. The
        // obvious optimization (modulus-x) generally is no good for x==0.
        return subtract(getZeroValue(), x);
    }


    // Returns the modular product of (the montgomery values) x and y.
    // Performance note: when calling this function and deciding which variable
    // to supply for its first argument and which to supply for its second
    // argument, you should prefer to use your variable which has changed less
    // recently for the second argument.  If possible, it is best to use a non-
    // loop carried dependency for the second argument.  This can improve
    // preformance for some Monty types.
    // Usually you don't want to specify PTAG (just accept the default).  For
    // advanced use: PTAG can be LowlatencyTag or LowuopsTag
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y) const
    {
        // note: the compiler should remove isZero calculations during dead code
        // elimination; isZero is unused here and impl.multiply is forced inline
        bool isZero;
        return impl.multiply(x, y, isZero, PTAG());
    }

    // This overload is an optimized equivalent of doing
    //    MontgomeryValue result = multiply(x, y);
    //    bool resultIsZero = (getCanonicalValue(result) == getZeroValue());
    // Performance note: same as given for multiply() above.
    // Note on the optimization: for some Monty Types (e.g. MontyQuarterRange),
    // setting resultIsZero via delegation to a lower level (as we do here) lets
    // us avoid an extra conditional move for getCanonicalValue().
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y,
                                                       bool& resultIsZero) const
    {
        MontgomeryValue ret = impl.multiply(x, y, resultIsZero, PTAG());
        HPBC_POSTCONDITION(resultIsZero ==
                                      (getCanonicalValue(ret)==getZeroValue()));
        return ret;
    }

    // "Fused multiply-subtract" operation:  Returns the modular evaluation of
    // (x * y) - z.  Note that z must be a canonical value.
    // Performance note 1: see the note for multiply(), which applies here too.
    // Performance note 2: This function usually has the benefits of both lower
    // latency and fewer uops than subtract(multiply(x, y), z).  You can expect
    // that this function will never be less efficient, and thus you should
    // almost always prefer to use it over the combination of subtract/multiply.
    // An exception to this rule about efficiency might occur if you are forced
    // to call getCanonicalValue for every time you call this function (i.e. you
    // are unable to reuse z), but even in that situation this function may
    // still have lower latency.
    // Usually you don't want to specify PTAG (just accept the default).  For
    // advanced use: PTAG can be LowlatencyTag or LowuopsTag
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y,
                                                         CanonicalValue z) const
    {
        MontgomeryValue ret = impl.fmsub(x, y, z, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                                getCanonicalValue(subtract(multiply(x, y), z)));
        return ret;
    }
    // The same comments as above apply, but this function uses a FusingValue
    // parameter z, which can be more efficient for some MontyTypes.  It is
    // never less efficient.
    // Note: when using a MontyType that receives a performance benefit from
    // using a FusingValue, we can note that getFusingValue() is slower than
    // getCanonicalValue() would have been.  This is rarely a reason to prefer
    // the overload of fmsub that uses a CanonicalValue though.  Normally
    // getFusingValue and getCanonicalValue would be called outside a loop, and
    // fmsub would be called inside the loop- which is the place performance
    // usually matters.
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y,
                                                            FusingValue z) const
    {
        MontgomeryValue ret = impl.fmsub(x, y, z, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                                getCanonicalValue(subtract(multiply(x, y), z)));
        return ret;
    }

    // "Fused multiply-add" operation:  Returns the modular evaluation of
    // (x * y) + z.  Note that z must be a canonical value.
    // Performance note: this function is generally more efficient than
    // add(multiply(x, y), z);  see the performance note for fmsub() - those
    // comments apply here too.  If you have the choice between fmsub or fmadd,
    // you should prefer to use fmsub since it uses either one less assembly
    // instruction or one less register, while having the same latency.  You can
    // change an fmadd into an fmsub if you negate z (subtraction by a negative
    // is equivalent to addition).  This is a beneficial change to make if you
    // have an fmadd in a loop which reuses the same value of z - you can
    // perform a single negate of z outside the loop and change the fmadd to
    // fmsub inside the loop.  You should not make this change if you are unable
    // to reuse z (i.e. don't negate z alongside every call to fmsub).
    // Usually you don't want to specify PTAG (just accept the default).  For
    // advanced use: PTAG can be LowlatencyTag or LowuopsTag
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
                                                         CanonicalValue z) const
    {
        MontgomeryValue ret = impl.fmadd(x, y, z, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                                     getCanonicalValue(add(multiply(x, y), z)));
        return ret;
    }
    // The same comments as above apply, but this function uses a FusingValue
    // parameter z, which can be more efficient for some MontyTypes.  It is
    // never less efficient.
    // Note: the comments above fmsub() regarding initializing a FusingValue
    // apply here as well.
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
                                                            FusingValue z) const
    {
        MontgomeryValue ret = impl.fmadd(x, y, z, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                                     getCanonicalValue(add(multiply(x, y), z)));
        return ret;
    }


    // Returns the modular evaluation of  x * x.
    // Performance notes:
    //   For some MontyTypes, this function can offer better performance than
    // multiply(x, x), so long as type T is the same size or smaller than the
    // CPU's native integer register size.  If type T is larger than the CPU's
    // integer register size, this function is likely to perform worse or no
    // better than multiply(x, x).
    template <class PTAG = LowlatencyTag>
    HURCHALLA_FORCE_INLINE MontgomeryValue square(MontgomeryValue x) const
    {
        MontgomeryValue ret = impl.square(x, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) == getCanonicalValue(
                                                               multiply(x, x)));
        return ret;
    }
    // "Fused square with subtract" operation.  Returns the modular evaluation
    // of (x * x) - cv.  Note that this function is usually but not always a
    // good performance replacement for fmsub(x, x, cv); see the comments above
    // square() for details.
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue fusedSquareSub(MontgomeryValue x, CanonicalValue cv) const
    {
        MontgomeryValue ret = impl.fusedSquareSub(x, cv, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                               getCanonicalValue(subtract(multiply(x, x), cv)));
        return ret;
    }
    // "Fused square with add" operation.  Returns the modular evaluation
    // of (x * x) + cv.  Note that this function is usually but not always a
    // good performance replacement for fmsub(x, x, cv); see the comments above
    // square() for details.
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue fusedSquareAdd(MontgomeryValue x, CanonicalValue cv) const
    {
        MontgomeryValue ret = impl.fusedSquareAdd(x, cv, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                                    getCanonicalValue(add(multiply(x, x), cv)));
        return ret;
    }


    // Calculates and returns the modular exponentiation of the montgomery value
    // 'base' to the power of (the type T variable) 'exponent'.
    MontgomeryValue pow(MontgomeryValue base, T exponent) const
    {
        HPBC_PRECONDITION(exponent >= 0);
        std::array<MontgomeryValue, 1> bases = {{ base }};
        std::array<MontgomeryValue, 1> result =
                detail::montgomery_array_pow<typename MontyType::MontyTag,
                                   MontgomeryForm>::pow(*this, bases, exponent);
        return result[0];
    }

    // This is a specially optimized version of the pow() function above.
    // It computes the results of multiple bases raised to the same power, and
    // takes advantage of CPU instruction level parallelism for efficiency.  You
    // can expect that calling this function multiple times with a small/optimal
    // value for NUM_BASES will be more efficient than calling this function a
    // single time using a large/suboptimal value for NUM_BASES.  Typically you
    // might expect an optimal NUM_BASES to be somewhere in the range of 3 to 6,
    // but you need to benchmark to find the best efficiency on your CPU.  FYI,
    // you should probably not expect to use a value of NUM_BASES significantly
    // greater than the number of (non-SIMD) integer multiply instructions that
    // can be simultaneously in-flight on your CPU on a single thread (via
    // pipelined and/or superscalar hardware multiply).  For example on the
    // Intel Skylake CPU, 3 IMUL instructions can be in-flight at once per core,
    // and so NUM_BASES == 3 is likely to be fairly close to optimal, whereas
    // NUM_BASES == 10 would likely be suboptimal; you need to benchmark to find
    // an optimal value though, and that value may differ from expectations.
    // Also note that as NUM_BASES increases, this function requires more and
    // more space in the instruction cache, and this can negatively affect
    // performance for the rest of your program.  Therefore, when measured
    // performance with two different values of NUM_BASES is similar, you should
    // almost always prefer the smaller NUM_BASES value.
    template <std::size_t NUM_BASES>
    std::array<MontgomeryValue, NUM_BASES>
    pow(std::array<MontgomeryValue, NUM_BASES>& bases, T exponent) const
    {
        HPBC_PRECONDITION(exponent >= 0);
        return detail::montgomery_array_pow<typename MontyType::MontyTag,
                                   MontgomeryForm>::pow(*this, bases, exponent);
    }


    // Returns the "greatest common divisor" of the standard representations
    // (non-montgomery) of both x and the modulus, using the supplied functor.
    // The functor must take two integral arguments of the same type and return
    // the gcd of its two arguments.  [Usually you would make the functor's
    // operator() a templated function, where the template parameter represents
    // the integral argument type.  Or more simply, you can just use a lambda,
    // with 'auto' type for its function parameters.]
    // Calling  gcd_with_modulus(x)  is more efficient than computing the
    // equivalent value  gcd_functor(convertOut(x), modulus).
    template <class F> HURCHALLA_FORCE_INLINE
    T gcd_with_modulus(MontgomeryValue x, const F& gcd_functor) const
    {
        return static_cast<T>(impl.gcd_with_modulus(x, gcd_functor));
    }


    // Returns  a % modulus.  A convenience function for better performance.
    // If you have already instantiated this MontgomeryForm, then calling
    // remainder() should be faster than directly computing  a % modulus,
    // even if your CPU has extremely fast division (like many new CPUs).
    T remainder(T a) const
    {
        HPBC_PRECONDITION(a >= 0);
        return static_cast<T>(impl.remainder(static_cast<U>(a)));
    }
};


} // end namespace

#endif
