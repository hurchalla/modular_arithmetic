// Copyright (c) 2024 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_NOFORCEINLINE_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_NOFORCEINLINE_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"

namespace hurchalla {


template<class MF>
class NoForceInlineMontgomeryForm final {
    MF impl;
    using T = typename MF::IntegerType;
public:
    using IntegerType = T;
    using MontgomeryValue = typename MF::MontgomeryValue;
    using CanonicalValue = typename MF::CanonicalValue;
    using FusingValue = typename MF::FusingValue;

    explicit NoForceInlineMontgomeryForm(T modulus) : impl(modulus) {}

    static constexpr T max_modulus()
    {
        return MF::max_modulus();
    }
    T getModulus() const
    {
        return impl.getModulus();
    }
    MontgomeryValue convertIn(T a) const
    {
        return impl.convertIn(a);
    }
    T convertOut(MontgomeryValue x) const
    {
        return impl.convertOut(x);
    }
    CanonicalValue getCanonicalValue(MontgomeryValue x) const
    {
        return impl.getCanonicalValue(x);
    }
    FusingValue getFusingValue(MontgomeryValue x) const
    {
        return impl.getFusingValue(x);
    }
    CanonicalValue getUnityValue() const
    {
        return impl.getUnityValue();
    }
    CanonicalValue getZeroValue() const
    {
        return impl.getZeroValue();
    }
    CanonicalValue getNegativeOneValue() const
    {
        return impl.getNegativeOneValue();
    }

    MontgomeryValue add(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.add(x, y);
    }
    MontgomeryValue add(MontgomeryValue x, CanonicalValue y) const
    {
        return impl.add(x, y);
    }
    MontgomeryValue add(CanonicalValue x, MontgomeryValue y) const
    {
        return impl.add(y, x);
    }
    CanonicalValue add(CanonicalValue x, CanonicalValue y) const
    {
        return impl.add(x, y);
    }

    MontgomeryValue subtract(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.subtract(x, y);
    }
    MontgomeryValue subtract(MontgomeryValue x, CanonicalValue y) const
    {
        return impl.subtract(x, y);
    }
    MontgomeryValue subtract(CanonicalValue x, MontgomeryValue y) const
    {
        return impl.subtract(x, y);
    }
    CanonicalValue subtract(CanonicalValue x, CanonicalValue y) const
    {
        return impl.subtract(x, y);
    }
/*
template <typename T> struct always_false {
    static constexpr bool value = false;
};
template <typename T>
void error_compile_type() const
{
    static_assert(always_false<T>::value, "");
}
*/
//hurchalla::detail::MontyQRValueTypes<unsigned long long>::V
    MontgomeryValue unorderedSubtract(MontgomeryValue x,
                                       MontgomeryValue y) const
    {
//        error_compile_type<MF>();
        return impl.unorderedSubtract(x, y);
    }
    MontgomeryValue unorderedSubtract(MontgomeryValue x,
                                       CanonicalValue y) const
    {
        return impl.unorderedSubtract(x, y);
    }
    MontgomeryValue unorderedSubtract(CanonicalValue x,
                                       MontgomeryValue y) const
    {
        return impl.unorderedSubtract(x, y);
    }

    MontgomeryValue negate(MontgomeryValue x) const
    {
        return impl.negate(x);
    }
    CanonicalValue negate(CanonicalValue x) const
    {
        return impl.negate(x);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y) const
    {
        return impl.template multiply<PTAG>(x, y);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y,
                                                       bool& resultIsZero) const
    {
        return impl.template multiply<PTAG>(x, y, resultIsZero);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y,
                                                         CanonicalValue z) const
    {
        return impl.template fmsub<PTAG>(x, y, z);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y,
                                                            FusingValue z) const
    {
        return impl.template fmsub<PTAG>(x, y, z);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
                                                         CanonicalValue z) const
    {
        return impl.template fmadd<PTAG>(x, y, z);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y,
                                                            FusingValue z) const
    {
        return impl.template fmadd<PTAG>(x, y, z);
    }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue square(MontgomeryValue x) const
    {
        return impl.template square<PTAG>(x);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fusedSquareSub(MontgomeryValue x, CanonicalValue cv) const
    {
        return impl.template fusedSquareSub<PTAG>(x, cv);
    }
    template <class PTAG = LowlatencyTag>
    MontgomeryValue fusedSquareAdd(MontgomeryValue x, CanonicalValue cv) const
    {
        return impl.template fusedSquareAdd<PTAG>(x, cv);
    }

    MontgomeryValue pow(MontgomeryValue base, T exponent) const
    {
        return impl.pow(base, exponent);
    }

    template <std::size_t NUM_BASES>
    std::array<MontgomeryValue, NUM_BASES>
    pow(std::array<MontgomeryValue, NUM_BASES>& bases, T exponent) const
    {
        return impl.pow(bases, exponent);
    }

    template <class F>
    T gcd_with_modulus(MontgomeryValue x, const F& gcd_functor) const
    {
        return impl.gcd_with_modulus(x, gcd_functor);
    }

    T remainder(T a) const
    {
        return impl.remainder(a);
    }
};


} // end namespace

#endif
