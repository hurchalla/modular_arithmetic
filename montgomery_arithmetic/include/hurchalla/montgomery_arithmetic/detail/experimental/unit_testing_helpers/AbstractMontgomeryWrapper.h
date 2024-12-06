// Copyright (c) 2024 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_ABSTRACT_MONTGOMERY_WRAPPER_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_ABSTRACT_MONTGOMERY_WRAPPER_H_INCLUDED


//#include "AbstractMontgomeryForm.h"
#include <memory>

namespace hurchalla {

// AMF should be AbstractMontgomeryForm<true> or AbstractMontgomeryForm<false>
template <class AMF>
class AbstractMontgomeryWrapper final {
    std::unique_ptr<const AMF> pimpl;
public:
    using IntegerType = typename AMF::IntegerType;
    using MontgomeryValue = typename AMF::MontgomeryValue;
    using CanonicalValue = typename AMF::CanonicalValue;
    using FusingValue = typename AMF::FusingValue;

    explicit AbstractMontgomeryWrapper(std::unique_ptr<const AMF> pimpl_)
        : pimpl(std::move(pimpl_)) {}

    IntegerType max_modulus() const { return pimpl->max_modulus(); }
    IntegerType getModulus() const { return pimpl->getModulus(); }
    MontgomeryValue convertIn(IntegerType a) const
        { return pimpl->convertIn(a); }
    IntegerType convertOut(MontgomeryValue x) const
        { return pimpl->convertOut(x); }
    CanonicalValue getCanonicalValue(MontgomeryValue x) const
        { return pimpl->getCanonicalValue(x); }
    FusingValue getFusingValue(MontgomeryValue x) const
        { return pimpl->getFusingValue(x); }
    CanonicalValue getUnityValue() const
        { return pimpl->getUnityValue(); }
    CanonicalValue getZeroValue() const
        { return pimpl->getZeroValue(); }
    CanonicalValue getNegativeOneValue() const
        { return pimpl->getNegativeOneValue(); }
    MontgomeryValue add(MontgomeryValue x, MontgomeryValue y) const
        { return pimpl->add(x, y); }
    MontgomeryValue add(MontgomeryValue x, CanonicalValue y) const
        { return pimpl->add(x, y); }
    MontgomeryValue add(CanonicalValue x, MontgomeryValue y) const
        { return pimpl->add(x, y); }
    CanonicalValue add(CanonicalValue x, CanonicalValue y) const
        { return pimpl->add(x, y); }
    MontgomeryValue subtract(MontgomeryValue x, MontgomeryValue y) const
        { return pimpl->subtract(x, y); }
    MontgomeryValue subtract(MontgomeryValue x, CanonicalValue y) const
        { return pimpl->subtract(x, y); }
    MontgomeryValue subtract(CanonicalValue x, MontgomeryValue y) const
        { return pimpl->subtract(x, y); }
    CanonicalValue subtract(CanonicalValue x, CanonicalValue y) const
        { return pimpl->subtract(x, y); }
    MontgomeryValue unorderedSubtract(MontgomeryValue x, MontgomeryValue y) const
        { return pimpl->unorderedSubtract(x, y); }
    MontgomeryValue unorderedSubtract(MontgomeryValue x, CanonicalValue y) const
        { return pimpl->unorderedSubtract(x, y); }
    MontgomeryValue unorderedSubtract(CanonicalValue x, MontgomeryValue y) const
        { return pimpl->unorderedSubtract(x, y); }
    MontgomeryValue negate(MontgomeryValue x) const
        { return pimpl->negate(x); }
    CanonicalValue negate(CanonicalValue x) const
        { return pimpl->negate(x); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y) const
        { return pimpl->template multiply<PTAG>(x, y); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue multiply(MontgomeryValue x, MontgomeryValue y, bool& resultIsZero) const
        { return pimpl->template multiply<PTAG>(x, y, resultIsZero); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y, CanonicalValue z) const
        { return pimpl->template fmsub<PTAG>(x, y, z); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmsub(MontgomeryValue x, MontgomeryValue y, FusingValue z) const
        { return pimpl->template fmsub<PTAG>(x, y, z); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y, CanonicalValue z) const
        { return pimpl->template fmadd<PTAG>(x, y, z); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fmadd(MontgomeryValue x, MontgomeryValue y, FusingValue z) const
        { return pimpl->template fmadd<PTAG>(x, y, z); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue square(MontgomeryValue x) const
        { return pimpl->template square<PTAG>(x); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fusedSquareSub(MontgomeryValue x, CanonicalValue cv) const
        { return pimpl->template fusedSquareSub<PTAG>(x, cv); }

    template <class PTAG = LowlatencyTag>
    MontgomeryValue fusedSquareAdd(MontgomeryValue x, CanonicalValue cv) const
        { return pimpl->template fusedSquareAdd<PTAG>(x, cv); }

    template <class F>
    IntegerType gcd_with_modulus(MontgomeryValue x, const F& gcd_functor) const
        { return pimpl->gcd_with_modulus(x, gcd_functor); }

    template <std::size_t NUM_BASES>
    std::array<MontgomeryValue, NUM_BASES>
    pow(std::array<MontgomeryValue, NUM_BASES>& bases, IntegerType exponent) const
        { return pimpl->pow(bases, exponent); }

    MontgomeryValue pow(MontgomeryValue base, IntegerType exponent) const
        { return pimpl->pow(base, exponent); }

    IntegerType remainder(IntegerType a) const
        { return pimpl->remainder(a); }
};


} // end namespace

#endif
