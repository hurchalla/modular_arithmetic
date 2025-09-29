// Copyright (c) 2024-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// This file exists to provide faster compilation of unit tests.
// Its contents probably shouldn't be used as an example for anything, except
// perhaps how to make unit tests compile faster when run time performance
// doesn't matter, via using polymorphism to de-templatize.


#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_ABSTRACT_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_ABSTRACT_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/is_equality_comparable.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/cselect_on_bit.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>
#include <array>
#include <vector>
#include <cstddef>
#include <cstdint>

namespace hurchalla {



// value types for AbstractMontgomeryForm
template <typename U>
struct AMFValueTypes {
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!(ut_numeric_limits<U>::is_signed), "");
    // regular montgomery value type
    class V {
        U value;
     protected:
        explicit V(U a) : value(a) {}
        U get() const { return value; }
        template<bool> friend class AbstractMontgomeryForm;
     public:
        // This next constructor purposely does not initialize 'value' - the
        // contents are undefined until the object is assigned to.
        V() = default;
        void cmov(bool cond, V v)
        {
            using PerfTag = CSelectDefaultTag;
              // value = cond ? v.value : value
            value = ::hurchalla::conditional_select<U, PerfTag>(cond,v.value,value);
        }
        template <int BITNUM>
        static V cselect_on_bit_ne0(uint64_t num, V v1, V v2)
        {
            U sel = ::hurchalla::cselect_on_bit<BITNUM>::ne_0(num, v1.get(), v2.get());
            return V(sel);
        }
        template <int BITNUM>
        static V cselect_on_bit_eq0(uint64_t num, V v1, V v2)
        {
            U sel = ::hurchalla::cselect_on_bit<BITNUM>::eq_0(num, v1.get(), v2.get());
            return V(sel);
        }
    };
    // canonical montgomery value type
    struct C : public V {
        C() = default;
        friend bool operator==(const C& x, const C& y)
            { return x.get() == y.get(); }
        friend bool operator!=(const C& x, const C& y)
            { return !(x == y); }

        template <int BITNUM>
        static C cselect_on_bit_ne0(uint64_t num, C c1, C c2)
        {
            U sel = ::hurchalla::cselect_on_bit<BITNUM>::ne_0(num, c1.get(), c2.get());
            return C(sel);
        }
        template <int BITNUM>
        static C cselect_on_bit_eq0(uint64_t num, C c1, C c2)
        {
            U sel = ::hurchalla::cselect_on_bit<BITNUM>::eq_0(num, c1.get(), c2.get());
            return C(sel);
        }
     protected:
        template<bool> friend class AbstractMontgomeryForm;
        explicit C(U a) : V(a) {}
    };
    // fusing montgomery value (addend/subtrahend for fmadd/fmsub)
    struct FV : public V {
        FV() = default;
     protected:
        template<bool> friend class AbstractMontgomeryForm;
        explicit FV(U a) : V(a) {}
    };
};




// This polymorphic base class is intended as an aid for unit testing; it
// is useful as a much faster compiling alternative to the MontgomeryForm
// template class.  Note that this class significantly sacrifices run-time
// performance for the sake of compilation speed.
//
template <bool useSignedT>
class AbstractMontgomeryForm {

protected:
#if (HURCHALLA_COMPILER_HAS_UINT128_T())
    using T = typename std::conditional<useSignedT, __int128_t, __uint128_t>::type;
#else
    using T = typename std::conditional<useSignedT, std::int64_t, std::uint64_t>::type;
#endif
    static_assert(ut_numeric_limits<T>::is_integer, "");

    using U = typename extensible_make_unsigned<T>::type;

public:
    using IntegerType = T;
//    using RU = typename extensible_make_unsigned<T>::type;

    using MontgomeryValue = typename AMFValueTypes<U>::V;
    static_assert(std::is_default_constructible<MontgomeryValue>::value, "");
    static_assert(std::is_copy_constructible<MontgomeryValue>::value, "");

    using CanonicalValue = typename AMFValueTypes<U>::C;
    static_assert(std::is_default_constructible<CanonicalValue>::value, "");
    static_assert(std::is_copy_constructible<CanonicalValue>::value, "");
    static_assert(is_equality_comparable<CanonicalValue>::value, "");
    static_assert(
               std::is_convertible<CanonicalValue, MontgomeryValue>::value, "");

    using FusingValue = typename AMFValueTypes<U>::FV;
    static_assert(std::is_default_constructible<FusingValue>::value, "");
    static_assert(std::is_copy_constructible<FusingValue>::value, "");
    static_assert(std::is_convertible<FusingValue, MontgomeryValue>::value, "");


    AbstractMontgomeryForm() = default;
    AbstractMontgomeryForm(const AbstractMontgomeryForm&) = default;
    AbstractMontgomeryForm& operator=(const AbstractMontgomeryForm&) = default;
    AbstractMontgomeryForm(AbstractMontgomeryForm&&) = default;
    AbstractMontgomeryForm& operator=(AbstractMontgomeryForm&&) = default;

    virtual ~AbstractMontgomeryForm() = default;


    virtual IntegerType max_modulus() const = 0;

    virtual IntegerType getModulus() const = 0;

    virtual CanonicalValue getCanonicalValue(MontgomeryValue x) const = 0;

    virtual FusingValue getFusingValue(MontgomeryValue x) const = 0;

    virtual CanonicalValue getUnityValue() const = 0;

    virtual CanonicalValue getZeroValue() const = 0;

    virtual CanonicalValue getNegativeOneValue() const = 0;

    virtual MontgomeryValue add(MontgomeryValue x, MontgomeryValue y) const = 0;

    virtual MontgomeryValue add(MontgomeryValue x, CanonicalValue y) const = 0;

    virtual MontgomeryValue add(CanonicalValue x, MontgomeryValue y) const = 0;

    virtual CanonicalValue add(CanonicalValue x, CanonicalValue y)
        const = 0;

    virtual MontgomeryValue unorderedSubtract(MontgomeryValue x,
        MontgomeryValue y) const = 0;
    
    virtual MontgomeryValue unorderedSubtract(MontgomeryValue x,
        CanonicalValue y) const = 0;

    virtual MontgomeryValue unorderedSubtract(CanonicalValue x,
        MontgomeryValue y) const = 0;

    virtual MontgomeryValue negate(MontgomeryValue x) const = 0;

    virtual CanonicalValue negate(CanonicalValue x) const = 0;

    virtual MontgomeryValue two_times(MontgomeryValue x) const = 0;

    virtual CanonicalValue two_times(CanonicalValue x) const = 0;

    virtual MontgomeryValue halve(MontgomeryValue x) const = 0;

    virtual CanonicalValue halve(CanonicalValue x) const = 0;

    virtual MontgomeryValue pow(MontgomeryValue base, IntegerType exponent)
        const = 0;

    virtual MontgomeryValue two_pow(IntegerType exponent) const = 0;

private:
    virtual MontgomeryValue convertIn(IntegerType a,
        bool useLowlatencyTag) const = 0;

    virtual IntegerType convertOut(MontgomeryValue x,
        bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue subtract(MontgomeryValue x, MontgomeryValue y,
        bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue subtract(MontgomeryValue x, CanonicalValue y,
        bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue subtract(CanonicalValue x, MontgomeryValue y,
        bool useLowlatencyTag) const = 0;

    virtual CanonicalValue subtract(CanonicalValue x, CanonicalValue y,
        bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue multiply2(MontgomeryValue x, MontgomeryValue y,
        bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue multiply2(MontgomeryValue x, MontgomeryValue y,
        bool& resultIsZero, bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y,
        CanonicalValue z, bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y,
        FusingValue z, bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
        CanonicalValue z, bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
        FusingValue z, bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue square(MontgomeryValue x,
        bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue fusedSquareSub(MontgomeryValue x, CanonicalValue cv,
        bool useLowlatencyTag) const = 0;

    virtual MontgomeryValue fusedSquareAdd(MontgomeryValue x, CanonicalValue cv,
        bool useLowlatencyTag) const = 0;

    virtual IntegerType remainder(IntegerType a,
        bool useLowlatencyTag) const = 0;

    virtual CanonicalValue inverse(MontgomeryValue x,
        bool useLowlatencyTag) const = 0;

    virtual std::vector<MontgomeryValue> vectorPow(
        const std::vector<MontgomeryValue>& bases, IntegerType exponent) const = 0;

    virtual IntegerType euclidean_gcd_with_modulus(MontgomeryValue x) const = 0;


public:
// adapters for functions that have template params; we wrap the virtual funcs
// since virtual funcs can't be templated.

    template <class PTAG = LowuopsTag>
    MontgomeryValue convertIn(IntegerType a) const
    {
        return convertIn(a, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <class PTAG = LowuopsTag>
    IntegerType convertOut(MontgomeryValue x) const
    {
        return convertOut(x, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <class PTAG = LowuopsTag>
    MontgomeryValue subtract(MontgomeryValue x, MontgomeryValue y) const
    {
        return subtract(x, y, std::is_same<PTAG, LowlatencyTag>::value);
    }
    template <class PTAG = LowuopsTag>
    MontgomeryValue subtract(MontgomeryValue x, CanonicalValue y) const
    {
        return subtract(x, y, std::is_same<PTAG, LowlatencyTag>::value);
    }
    template <class PTAG = LowuopsTag>
    MontgomeryValue subtract(CanonicalValue x, MontgomeryValue y) const
    {
        return subtract(x, y, std::is_same<PTAG, LowlatencyTag>::value);
    }
    template <class PTAG = LowuopsTag>
    CanonicalValue subtract(CanonicalValue x, CanonicalValue y) const
    {
        return subtract(x, y, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y) const
    {
        return multiply2(x, y, std::is_same<PTAG, LowlatencyTag>::value);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y,
        bool& resultIsZero) const
    {
        return multiply2(x, y, resultIsZero, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y,
        CanonicalValue z) const
    {
        return fmsub(x, y, z, std::is_same<PTAG, LowlatencyTag>::value);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y,
        FusingValue z) const
    {
        return fmsub(x, y, z, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
        CanonicalValue z) const
    {
        return fmadd(x, y, z, std::is_same<PTAG, LowlatencyTag>::value);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
        FusingValue z) const
    {
        return fmadd(x, y, z, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue square(MontgomeryValue x) const
    {
         return square(x, std::is_same<PTAG, LowlatencyTag>::value);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fusedSquareSub(MontgomeryValue x, CanonicalValue cv) const
    {
        return fusedSquareSub(x, cv, std::is_same<PTAG, LowlatencyTag>::value);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fusedSquareAdd(MontgomeryValue x, CanonicalValue cv) const
    {
        return fusedSquareAdd(x, cv, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <class PTAG = LowlatencyTag>
    CanonicalValue inverse(MontgomeryValue x) const
    {
        return inverse(x, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <class PTAG = LowlatencyTag>
    IntegerType remainder(IntegerType a) const
    {
        return remainder(a, std::is_same<PTAG, LowlatencyTag>::value);
    }

    template <std::size_t NUM_BASES>
    std::array<MontgomeryValue, NUM_BASES>
    pow(const std::array<MontgomeryValue, NUM_BASES>& bases, IntegerType exponent) const
    {
        std::vector<MontgomeryValue> bases_vec(bases.begin(), bases.end());
        auto answer_vec = vectorPow(bases_vec, exponent);
        std::array<MontgomeryValue, NUM_BASES> answer;
        for (std::size_t i=0; i<NUM_BASES; ++i)
            answer[i] = answer_vec[i];
        return answer;
    }

    template <class F>
    T gcd_with_modulus(MontgomeryValue x, const F& /* gcd_functor, ignored */) const
    {
        return euclidean_gcd_with_modulus(x);
    }

};


} // end namespace

#endif
