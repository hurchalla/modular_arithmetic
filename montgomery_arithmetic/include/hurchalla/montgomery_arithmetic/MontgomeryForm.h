// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/MontgomeryDefault.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


// When using the default MontyType, T must be signed or unsigned integral type.
// A custom MontyType may have different requirements for type T (e.g. that T is
// an unsigned integral type).
template<typename T, class MontyType = typename MontgomeryDefault<T>::type>
class MontgomeryForm {
    const MontyType impl;
    using U = typename MontyType::template_param_type;
    static_assert(modular_arithmetic::ma_numeric_limits<U>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<U>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<U>::digits >=
                  modular_arithmetic::ma_numeric_limits<T>::digits, "");
public:
    using MontgomeryValue = typename MontyType::montvalue_type;
    using T_type = T;

    class CanonicalValue : public MontgomeryValue {
        friend MontgomeryForm;
        explicit CanonicalValue(MontgomeryValue val) : MontgomeryValue(val) {}
    public:
        //CanonicalValue() : MontgomeryValue() {}
        friend bool operator==(const CanonicalValue& x, const CanonicalValue& y)
        {
            return x.value == y.value;
        }
        friend bool operator!=(const CanonicalValue& x, const CanonicalValue& y)
        {
            return !(x == y);
        }
    };

    explicit MontgomeryForm(T modulus) : impl(static_cast<U>(modulus))
    {
        HPBC_PRECONDITION(modulus % 2 == 1);  // modulus must be odd
        HPBC_PRECONDITION(modulus > 1);
    }
    MontgomeryForm(const MontgomeryForm&) = delete;
    MontgomeryForm& operator=(const MontgomeryForm&) = delete;

    // Returns the largest valid modulus allowed for the constructor.
    static constexpr T max_modulus()
    {
        return (MontyType::max_modulus() >
            static_cast<U>(modular_arithmetic::ma_numeric_limits<T>::max()))
            ? ((modular_arithmetic::ma_numeric_limits<T>::max() == 0)
                ? static_cast<T>(modular_arithmetic::ma_numeric_limits<T>::max()
                                 - 1)
                : modular_arithmetic::ma_numeric_limits<T>::max())
            : static_cast<T>(MontyType::max_modulus());
    }

    // Returns the modulus given to the constructor
    T getModulus() const { return static_cast<T>(impl.getModulus()); }

    // Returns the converted value of the standard number 'a' into monty form.
    // Requires a >= 0.
    MontgomeryValue convertIn(T a) const
    {
        HPBC_PRECONDITION(a >= 0);
        return impl.convertIn(static_cast<U>(a));
    }

    // Converts (montgomery value) x into a "normal" number; returns the result.
    // Guarantees 0 <= result < modulus.
    T convertOut(MontgomeryValue x) const
    {
        T a = static_cast<T>(impl.convertOut(x));
        HPBC_POSTCONDITION(0 <= a && a < static_cast<T>(impl.getModulus()));
        return a;
    }

    // Returns a unique (canonical) value representing the equivalence class of
    // x modulo the modulus.  You can not directly compare MontgomeryValues, but
    // you can call getCanonicalValue(), and then use standard equality or
    // inequality operators to compare the resulting CanonicalValues.
    CanonicalValue getCanonicalValue(MontgomeryValue x) const
    {
        MontgomeryValue ret = impl.getCanonicalValue(x);
        return CanonicalValue(ret);
    }
    // Returns the canonical monty value that represents the type T value 1.
    // The call is equivalent to getCanonicalValue(convertIn(static_cast<T>(1)))
    // but it's more efficient (essentially zero cost) and more convenient.
    CanonicalValue getUnityValue() const
    {
        MontgomeryValue ret = impl.getUnityValue();
        HPBC_ASSERT(impl.isCanonical(ret));
        return CanonicalValue(ret);
    }
    // Returns the canonical monty value that represents the type T value 0.
    // The call is equivalent to getCanonicalValue(convertIn(static_cast<T>(0)))
    // but it's more efficient (essentially zero cost) and more convenient.
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
    MontgomeryValue add(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.add(x, y);
    }
    // Returns the modular sum of the montgomery value x and the *canonical*
    // value y.  Performance note: this function is sometimes more efficient
    // than the above add() function and should never be less efficient, so it
    // can be useful to call getCanonicalValue outside of a loop to get 'y' as
    // a CanonicalValue, when you are calling add() inside a loop.
    MontgomeryValue add(MontgomeryValue x, CanonicalValue y) const
    {
        MontgomeryValue ret = impl.add_canonical_value(x, y);
        HPBC_ASSERT(getCanonicalValue(ret) ==
                    getCanonicalValue(add(x, static_cast<MontgomeryValue>(y))));
        return ret;
    }

    // Returns the modular difference of (the montgomery values) x and y.  More
    // precisely, x minus y.
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
    // don't care which subtraction gets performed, unordered_subtract() will
    // usually perform slightly better than this->subtract() [in terms of number
    // of instructions, number of registers used, and possibly total latency].
    MontgomeryValue unordered_subtract(MontgomeryValue x,
                                       MontgomeryValue y) const
    {
        return impl.unordered_subtract(x, y);
    }

    // Returns the modular negation of the montgomery value x.
    MontgomeryValue negate(MontgomeryValue x) const
    {
        return subtract(getZeroValue(), x);
    }
    // Returns the modular negation of the canonical value x.
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
    template <class PTAG = LowlatencyTag>
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
        // This is an optimized version of Algorithm 14.76, from
        // Applied Handbook of Cryptography- http://cacr.uwaterloo.ca/hac/
        // See also: hurchalla/modular_arithmetic/internal/impl_modular_pow.h
        MontgomeryValue result = (exponent & static_cast<T>(1)) ?
                                                    base : impl.getUnityValue();
        while (exponent > static_cast<T>(1))
        {
            exponent = static_cast<T>(exponent >> static_cast<T>(1));
            base = multiply<LowuopsTag>(base, base);
            // The multiply above is a loop carried dependency.  Thus, a second
            // loop carried dependency with the same length can be essentially
            // free due to instruction level parallelism, so long as it does not
            // introduce any branch mispredictions.
            // So we will always compute the second multiply, instead of
            // conditionally computing it, and we will encourage the compiler to
            // use a (branchless) conditional move instruction.
            // We use lowlatencyTag below because the "result" loop carried
            // dependency depends upon both multiply and a conditional move,
            // whereas "base" above depends only on multiply and thus is tagged
            // for lowuops since it is less likely to be a latency bottleneck.
            MontgomeryValue tmp = multiply<LowlatencyTag>(result, base);
            result = (exponent & static_cast<T>(1)) ? tmp : result;
        }
        return result;
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
    template <class PTAG = LowlatencyTag>
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
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
                                                         CanonicalValue z) const
    {
        MontgomeryValue ret = impl.fmadd(x, y, z, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                                     getCanonicalValue(add(multiply(x, y), z)));
        return ret;
    }

    // "Fused add-multiply" operation:  Returns the modular evaluation of
    // (x + y) * z.  Note that y must be a canonical value.
    // Performance note: for certain MontyTypes (e.g. MontyHalfRange and
    // MontySixthRange) this function has the twin benefits of lower latency and
    // fewer uops than multiply(add(x, y), z).  For all other MontyTypes, you
    // can expect this function's efficiency will be at least as good as
    // multiply(add(x, y), z).  Since there are potential benefits without any
    // expected downside, you should prefer to use this function.
    // FSMUL note: If what you really want/need is an fsmul (fused
    // subtract-multiply), you can negate y prior to this call to get its
    // equivalent; see the comments above fmadd for why you should perform
    // negate outside of a loop.
    // Usually you don't want to specify PTAG (just accept the default).  For
    // advanced use: PTAG can be LowlatencyTag or LowuopsTag
    template <class PTAG = LowlatencyTag>
    MontgomeryValue famul(MontgomeryValue x, CanonicalValue y,
                                                        MontgomeryValue z) const
    {
        // note: the compiler should remove isZero calculations during dead code
        // elimination; isZero is unused here and impl.famul() is forced inline.
        bool isZero;
        MontgomeryValue ret = impl.famul(x, y, z, isZero, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                                     getCanonicalValue(multiply(add(x, y), z)));
        return ret;
    }

    // "Fused add-multiply, with test for result congruent to zero":
    // This peculiar function is an optimization that is useful for Pollard-rho
    // factorization, though perhaps it might be useful elsewhere too.
    // The operation's results are identical to
    //    MontgomeryValue result = famul(x, y, z);
    //    bool isZero = (getCanonicalValue(result) == getZeroValue());
    // Note on the optimization: for certain Monty Types (MontyQuarterRange and
    // MontySixthRange), setting isZero via delegation to a lower level (as we
    // do here) lets us avoid an extra conditional move for getCanonicalValue().
    template <class PTAG = LowlatencyTag>
    MontgomeryValue famulIsZero(MontgomeryValue x, CanonicalValue y,
                                          MontgomeryValue z, bool& isZero) const
    {
        MontgomeryValue ret = impl.famul(x, y, z, isZero, PTAG());
        HPBC_POSTCONDITION(getCanonicalValue(ret) ==
                                     getCanonicalValue(multiply(add(x, y), z)));
        HPBC_POSTCONDITION(isZero == (getCanonicalValue(ret)==getZeroValue()));
        return ret;
    }
    // "Multiply, with test for result congruent to zero".  See comments on
    // famulIsZero() for rationale.
    template <class PTAG = LowlatencyTag>
    MontgomeryValue multiplyIsZero(MontgomeryValue x, MontgomeryValue y,
                                                             bool& isZero) const
    {
        MontgomeryValue ret = impl.multiply(x, y, isZero, PTAG());
        HPBC_POSTCONDITION(isZero == (getCanonicalValue(ret)==getZeroValue()));
        return ret;
    }
};


}} // end namespace

#endif
