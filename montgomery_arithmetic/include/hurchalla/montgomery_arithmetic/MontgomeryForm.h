// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontgomeryDefault.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/montgomery_pow.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <array>
#include <cstddef>

namespace hurchalla {


// When using the default MontyType, T must be signed or unsigned integral type.
// A custom MontyType may have different requirements for type T (e.g. that T is
// an unsigned integral type).
template<class T, class MontyType = typename detail::MontgomeryDefault<T>::type>
class MontgomeryForm {
    const MontyType impl;
    using U = typename MontyType::uint_type;
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!(ut_numeric_limits<U>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<U>::digits >=
                  ut_numeric_limits<T>::digits, "");
public:
    using MontgomeryValue = typename MontyType::montvalue_type;
    using T_type = T;

    class CanonicalValue : public MontgomeryValue {
        friend MontgomeryForm;
        HURCHALLA_FORCE_INLINE
        explicit CanonicalValue(MontgomeryValue val) : MontgomeryValue(val) {}
    public:
        HURCHALLA_FORCE_INLINE CanonicalValue() : MontgomeryValue() {}
        HURCHALLA_FORCE_INLINE
        friend bool operator==(const CanonicalValue& x, const CanonicalValue& y)
        {
            return x.value == y.value;
        }
        HURCHALLA_FORCE_INLINE
        friend bool operator!=(const CanonicalValue& x, const CanonicalValue& y)
        {
            return !(x == y);
        }
    };

    explicit MontgomeryForm(T modulus) : impl(static_cast<U>(modulus))
    {
        // usually odd modulus is required, but not for MontyWrappedStandardMath
        //HPBC_PRECONDITION(modulus % 2 == 1);
        HPBC_PRECONDITION(modulus > 1);
    }
    MontgomeryForm(const MontgomeryForm&) = delete;
    MontgomeryForm& operator=(const MontgomeryForm&) = delete;

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
        MontgomeryValue ret = impl.getCanonicalValue(x);
        return CanonicalValue(ret);
    }
    // Returns the canonical monty value that represents the type T value 1.
    // The call is equivalent to getCanonicalValue(convertIn(static_cast<T>(1)))
    // but it's more efficient (essentially zero cost) and more convenient.
    HURCHALLA_FORCE_INLINE
    CanonicalValue getUnityValue() const
    {
        MontgomeryValue ret = impl.getUnityValue();
        HPBC_ASSERT(impl.isCanonical(ret));
        return CanonicalValue(ret);
    }
    // Returns the canonical monty value that represents the type T value 0.
    // The call is equivalent to getCanonicalValue(convertIn(static_cast<T>(0)))
    // but it's more efficient (essentially zero cost) and more convenient.
    HURCHALLA_FORCE_INLINE
    CanonicalValue getZeroValue() const
    {
        MontgomeryValue ret = impl.getZeroValue();
        HPBC_ASSERT(impl.isCanonical(ret));
        return CanonicalValue(ret);
    }
    // Returns the canonical monty value that represents the type T value
    // modulus-1 (which equals -1 (mod modulus)).  The call is equivalent to
    // getCanonicalValue(convertIn(static_cast<T>(modulus - 1))), but it's more
    // efficient (essentially zero cost) and more convenient.
    HURCHALLA_FORCE_INLINE
    CanonicalValue getNegativeOneValue() const
    {
        MontgomeryValue ret = impl.getNegativeOneValue();
        HPBC_ASSERT(impl.isCanonical(ret));
        return CanonicalValue(ret);
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
    // Returns the modular sum of the montgomery value x and the *canonical*
    // value y.  Performance note: this function is sometimes more efficient
    // than the above add() function and should never be less efficient, so it
    // can be useful to call getCanonicalValue outside of a loop to get 'y' as
    // a CanonicalValue, when you are calling add() inside a loop.
    HURCHALLA_FORCE_INLINE
    MontgomeryValue add(MontgomeryValue x, CanonicalValue y) const
    {
        MontgomeryValue ret = impl.add_canonical_value(x, y);
        HPBC_ASSERT(getCanonicalValue(ret) ==
                    getCanonicalValue(add(x, static_cast<MontgomeryValue>(y))));
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
    // *canonical* value y.  Performance note: this function is sometimes more
    // efficient than the above subtract() function and should never be less
    // efficient, so it can be useful to call getCanonicalValue outside of a
    // loop to get 'y' as a CanonicalValue, when you are calling subtract()
    // inside a loop.
    HURCHALLA_FORCE_INLINE
    MontgomeryValue subtract(MontgomeryValue x, CanonicalValue y) const
    {
        MontgomeryValue ret = impl.subtract_canonical_value(x, y);
        HPBC_ASSERT(getCanonicalValue(ret) ==
               getCanonicalValue(subtract(x, static_cast<MontgomeryValue>(y))));
        // all implementations of subtract_canonical_value() guarantee:
        HPBC_ASSERT(impl.isCanonical(x) ? impl.isCanonical(ret) : true);
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

    // Returns the modular negation of the montgomery value x.
    HURCHALLA_FORCE_INLINE
    MontgomeryValue negate(MontgomeryValue x) const
    {
        return subtract(getZeroValue(), x);
    }
    // Returns the modular negation of the canonical value x.
    HURCHALLA_FORCE_INLINE
    CanonicalValue negate(CanonicalValue x) const
    {
        MontgomeryValue ret = subtract(getZeroValue(), x);
        // Since both getZeroValue() and x are canonical, subtract() guarantees:
        HPBC_ASSERT(impl.isCanonical(ret));
        return CanonicalValue(ret);
    }

    // Returns the modular product of (the montgomery values) x and y.
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

    // Calculates and returns the modular exponentiation of the montgomery value
    // 'base' to the power of (the type T variable) 'exponent'.
    MontgomeryValue pow(MontgomeryValue base, T exponent) const
    {
        HPBC_PRECONDITION(exponent >= 0);
        return detail::montgomery_pow(*this, base, exponent);
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
    // more space in the instruction cache, and this can negatively effect
    // performance for the rest of your program.  Therefore, when measured
    // performance with two different values of NUM_BASES is similar, you should
    // almost always prefer the smaller NUM_BASES value.
    template <std::size_t NUM_BASES>
    std::array<MontgomeryValue, NUM_BASES>
    pow(std::array<MontgomeryValue, NUM_BASES>& bases, T exponent) const
    {
        HPBC_PRECONDITION(exponent >= 0);
        return detail::montgomery_pow(*this, bases, exponent);
    }

    // Returns the greatest common denominator of the standard representations
    // (non-montgomery) of both x and the modulus, using the supplied functor.
    // The functor must take two integral arguments of the same type and return
    // the gcd of its two arguments.
    // Calling  gcd_with_modulus(x)  is more efficient than computing the
    // equivalent value  gcd_functor(convertOut(x), modulus).
    template <class F> HURCHALLA_FORCE_INLINE
    T gcd_with_modulus(MontgomeryValue x, const F& gcd_functor) const
    {
        return static_cast<T>(impl.gcd_with_modulus(x, gcd_functor));
    }

    // "Fused multiply-subtract" operation:  Returns the modular evaluation of
    // (x * y) - z.  Note that z must be a canonical value.
    // Performance note: This function usually has the benefits of both lower
    // latency and fewer uops than subtract(multiply(x, y), z).  You can expect
    // that this function will never be less efficient, and thus you should
    // almost always prefer to use it over the combination of subtract/multiply.
    // An exception to this rule about efficiency might occur if you are forced
    // to call getCanonicalValue for every time you call this function (i.e. you
    // are unable to reuse z), but even in that situation this function will
    // likely have lower latency.
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

    // "Multiply, with test for result congruent to zero":
    // The operation's results are identical to
    //    MontgomeryValue result = multiply(x, y);
    //    bool isZero = (getCanonicalValue(result) == getZeroValue());
    // Note on the optimization: for certain Monty Types (e.g.
    // MontyQuarterRange), setting isZero via delegation to a lower level (as we
    // do here) lets us avoid an extra conditional move for getCanonicalValue().
    template <class PTAG = LowlatencyTag> HURCHALLA_FORCE_INLINE
    MontgomeryValue multiplyIsZero(MontgomeryValue x, MontgomeryValue y,
                                                             bool& isZero) const
    {
        MontgomeryValue ret = impl.multiply(x, y, isZero, PTAG());
        HPBC_POSTCONDITION(isZero == (getCanonicalValue(ret)==getZeroValue()));
        return ret;
    }
};


} // end namespace

#endif
