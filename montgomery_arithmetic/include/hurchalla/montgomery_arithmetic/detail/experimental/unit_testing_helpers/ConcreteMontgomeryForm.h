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


#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_CONCRETE_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_CONCRETE_MONTGOMERY_FORM_H_INCLUDED


#include <utility>


#include "AbstractMontgomeryForm.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/is_equality_comparable.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>
#include <array>
#include <vector>
#include <cstddef>
#include <utility>

namespace hurchalla {


template <class MF, std::size_t... POW_ARRAY_SIZES>
class ConcreteMontgomeryForm final : public AbstractMontgomeryForm<ut_numeric_limits<typename MF::IntegerType>::is_signed> {
    const MF mf;

public:
    using Parent = AbstractMontgomeryForm<ut_numeric_limits<typename MF::IntegerType>::is_signed>;
    using T = typename Parent::IntegerType;
    using V = typename Parent::MontgomeryValue;
    using C = typename Parent::CanonicalValue;
    using FV = typename Parent::FusingValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");

private:
    using U = typename extensible_make_unsigned<T>::type;

    struct OpenV : public V {
#ifndef _MSC_VER
        auto get() const -> decltype(V::get()) { return V::get(); }
        // for explanation of OT declaration, see
        // https://stackoverflow.com/questions/26435084/how-to-get-the-return-type-of-a-member-function-from-within-a-class
        using OT = decltype((std::declval<OpenV>().*std::declval<decltype(&OpenV::get)>())());
#else
        using OT = decltype((std::declval<V>().*std::declval<decltype(&V::get)>())());
        OT get() const { return V::get(); }
#endif
        OpenV() = default;
        explicit OpenV(OT a) : V(a) {}
        explicit OpenV(V x) : V(x) {}
    };
    struct OpenC : public C {
#ifndef _MSC_VER
        auto get() const -> decltype(C::get()) { return C::get(); }
        using OT = decltype((std::declval<OpenC>().*std::declval<decltype(&OpenC::get)>())());
#else
        using OT = decltype((std::declval<C>().*std::declval<decltype(&C::get)>())());
        OT get() const { return C::get(); }
#endif
        OpenC() = default;
        explicit OpenC(OT a) : C(a) {}
        explicit OpenC(C x) : C(x) {}
    };
    struct OpenFV : public FV {
#ifndef _MSC_VER
        auto get() const -> decltype(FV::get()) { return FV::get(); }
        using OT = decltype((std::declval<OpenFV>().*std::declval<decltype(&OpenFV::get)>())());
#else
        using OT = decltype((std::declval<FV>().*std::declval<decltype(&FV::get)>())());
        OT get() const { return FV::get(); }
#endif
        OpenFV() = default;
        explicit OpenFV(OT a) : FV(a) {}
        explicit OpenFV(FV x) : FV(x) {}
    };

    using MFT = typename MF::IntegerType;
    using MFV = typename MF::MontgomeryValue;
    using MFC = typename MF::CanonicalValue;
    using MFFV = typename MF::FusingValue;
    static_assert(ut_numeric_limits<MFT>::is_integer, "");
    static_assert(ut_numeric_limits<T>::max() >= ut_numeric_limits<MFT>::max() , "");
    static_assert(ut_numeric_limits<T>::min() <= ut_numeric_limits<MFT>::min() , "");
    static_assert((ut_numeric_limits<T>::is_signed && ut_numeric_limits<MFT>::is_signed)
          || (!ut_numeric_limits<T>::is_signed && !ut_numeric_limits<MFT>::is_signed) , "");


    struct OpenMFV : public MFV {
#ifndef _MSC_VER
        auto get() const -> decltype(MFV::get()) { return MFV::get(); }
        using OT = decltype((std::declval<OpenMFV>().*std::declval<decltype(&OpenMFV::get)>())());
#else
        using OT = decltype((std::declval<MFV>().*std::declval<decltype(&MFV::get)>())());
        OT get() const { return MFV::get(); }
#endif
        OpenMFV() = default;
        explicit OpenMFV(OT a) : MFV(a) {}
        explicit OpenMFV(MFV x) : MFV(x) {}
        explicit OpenMFV(OpenV x) : MFV(static_cast<OT>(x.get()))
        {
            static_assert(!ut_numeric_limits<typename OpenV::OT>::is_signed, "");
            if HURCHALLA_CPP17_CONSTEXPR (ut_numeric_limits<OT>::is_signed) {
                using S = typename extensible_make_signed<typename OpenV::OT>::type;
                S s = static_cast<S>(x.get());
                HPBC_CLOCKWORK_ASSERT2(ut_numeric_limits<OT>::min() <= s &&
                              s <= ut_numeric_limits<OT>::max());
            } else {
                HPBC_CLOCKWORK_ASSERT2(0 <= x.get() && x.get() <= ut_numeric_limits<OT>::max());
            }
            static_assert(static_cast<OT>(-1) ==
                static_cast<OT>(static_cast<typename OpenV::OT>(static_cast<OT>(-1))), "");
        }
    };
    struct OpenMFC : public MFC {
#ifndef _MSC_VER
        auto get() const -> decltype(MFC::get()) { return MFC::get(); }
        using OT = decltype((std::declval<OpenMFC>().*std::declval<decltype(&OpenMFC::get)>())());
#else
        using OT = decltype((std::declval<MFC>().*std::declval<decltype(&MFC::get)>())());
        OT get() const { return MFC::get(); }
#endif
        OpenMFC() = default;
        explicit OpenMFC(OT a) : MFC(a) {}
        explicit OpenMFC(MFC x) : MFC(x) {}
        explicit OpenMFC(OpenC x) : MFC(static_cast<OT>(x.get()))
        {
            static_assert(!ut_numeric_limits<typename OpenC::OT>::is_signed, "");
            if HURCHALLA_CPP17_CONSTEXPR (ut_numeric_limits<OT>::is_signed) {
                using S = typename extensible_make_signed<typename OpenC::OT>::type;
                S s = static_cast<S>(x.get());
                HPBC_CLOCKWORK_ASSERT2(ut_numeric_limits<OT>::min() <= s &&
                              s <= ut_numeric_limits<OT>::max());
            } else {
                HPBC_CLOCKWORK_ASSERT2(0 <= x.get() && x.get() <= ut_numeric_limits<OT>::max());
            }
            static_assert(static_cast<OT>(-1) ==
                static_cast<OT>(static_cast<typename OpenC::OT>(static_cast<OT>(-1))), "");
        }
    };
    struct OpenMFFV : public MFFV {
#ifndef _MSC_VER
        auto get() const -> decltype(MFFV::get()) { return MFFV::get(); }
        using OT = decltype((std::declval<OpenMFFV>().*std::declval<decltype(&OpenMFFV::get)>())());
#else
        using OT = decltype((std::declval<MFFV>().*std::declval<decltype(&MFFV::get)>())());
        OT get() const { return MFFV::get(); }
#endif
        OpenMFFV() = default;
        explicit OpenMFFV(OT a) : MFFV(a) {}
        explicit OpenMFFV(MFFV x) : MFFV(x) {}
        explicit OpenMFFV(OpenFV x) : MFFV(static_cast<OT>(x.get()))
        {
            static_assert(!ut_numeric_limits<typename OpenFV::OT>::is_signed, "");
            if HURCHALLA_CPP17_CONSTEXPR (ut_numeric_limits<OT>::is_signed) {
                using S = typename extensible_make_signed<typename OpenFV::OT>::type;
                S s = static_cast<S>(x.get());
                HPBC_CLOCKWORK_ASSERT2(ut_numeric_limits<OT>::min() <= s &&
                              s <= ut_numeric_limits<OT>::max());
            } else {
                HPBC_CLOCKWORK_ASSERT2(0 <= x.get() && x.get() <= ut_numeric_limits<OT>::max());
            }
            static_assert(static_cast<OT>(-1) ==
                static_cast<OT>(static_cast<typename OpenFV::OT>(static_cast<OT>(-1))), "");
        }
    };

public:

    explicit ConcreteMontgomeryForm(T modulus) : mf(static_cast<MFT>(modulus))
    {
        HPBC_CLOCKWORK_PRECONDITION2(ut_numeric_limits<MFT>::min() <= modulus &&
                           modulus <= ut_numeric_limits<MFT>::max());
    }

    ConcreteMontgomeryForm(const ConcreteMontgomeryForm&) = delete;
    ConcreteMontgomeryForm& operator=(const ConcreteMontgomeryForm&) = delete;
    ConcreteMontgomeryForm(ConcreteMontgomeryForm&&) = delete;
    ConcreteMontgomeryForm& operator=(ConcreteMontgomeryForm&&) = delete;

    virtual ~ConcreteMontgomeryForm() override = default;


    virtual T max_modulus() const override
    {
        return MF::max_modulus();
    }

    virtual T getModulus() const override
    {
        return mf.getModulus();
    }

    virtual C getCanonicalValue(V x) const override
    {
        OpenMFC mfc(mf.getCanonicalValue(OpenMFV(OpenV(x))));
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual FV getFusingValue(V x) const override
    {
        OpenMFFV mffv(mf.getFusingValue(OpenMFV(OpenV(x))));
        // note: mffv.get() might be signed or unsigned; OpenFV::OT is unsigned
        return OpenFV(static_cast<typename OpenFV::OT>(mffv.get()));
    }

    virtual C getUnityValue() const override
    {
        OpenMFC mfc(mf.getUnityValue());
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual C getZeroValue() const override
    {
        OpenMFC mfc(mf.getZeroValue());
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual C getNegativeOneValue() const override
    {
        OpenMFC mfc(mf.getNegativeOneValue());
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual V add(V x, V y) const override
    {
        OpenMFV mfv(mf.add(OpenMFV(OpenV(x)), OpenMFV(OpenV(y))));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V add(V x, C y) const override
    {
        OpenMFV mfv(mf.add(OpenMFV(OpenV(x)), OpenMFC(OpenC(y))));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V add(C x, V y) const override
    {
        OpenMFV mfv(mf.add(OpenMFC(OpenC(x)), OpenMFV(OpenV(y))));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual C add(C x, C y) const override
    {
        OpenMFC mfc(mf.add(OpenMFC(OpenC(x)), OpenMFC(OpenC(y))));
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual V unorderedSubtract(V x, V y) const override
    {
        OpenMFV mfv(mf.unorderedSubtract(OpenMFV(OpenV(x)), OpenMFV(OpenV(y))));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V unorderedSubtract(V x, C y) const override
    {
        OpenMFV mfv(mf.unorderedSubtract(OpenMFV(OpenV(x)), OpenMFC(OpenC(y))));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V unorderedSubtract(C x, V y) const override
    {
        OpenMFV mfv(mf.unorderedSubtract(OpenMFC(OpenC(x)), OpenMFV(OpenV(y))));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V negate(V x) const override
    {
        OpenMFV mfv(mf.negate(OpenMFV(OpenV(x))));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual C negate(C x) const override
    {
        OpenMFC mfc(mf.negate(OpenMFC(OpenC(x))));
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual V two_times(V x) const override
    {
        OpenMFV mfv(mf.two_times(OpenMFV(OpenV(x))));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual C two_times(C x) const override
    {
        OpenMFC mfc(mf.two_times(OpenMFC(OpenC(x))));
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual V two_pow(T exponent) const override
    {
        HPBC_CLOCKWORK_PRECONDITION2(0 <= exponent);
        if (ut_numeric_limits<T>::max() > ut_numeric_limits<MFT>::max()) {
            // kind of an unavoidable hack, so that AbstractMontgomeryForm has
            // the same contract for two_pow() as MongomeryForm, which allows
            // exponent to have any value of T >= 0.
            static constexpr MFT mft_max = ut_numeric_limits<MFT>::max();
            if (exponent > mft_max) {
                MFV maxpow = mf.two_pow(mft_max);
                MFV accum = mf.getUnityValue();
                do {
                    accum = mf.multiply(accum, maxpow);
                    exponent -= mft_max;
                } while (exponent > static_cast<T>(mft_max));
                accum = mf.multiply(accum, mf.two_pow(static_cast<MFT>(exponent)));
                OpenMFV result(accum);
                // note: result.get() might be signed or unsigned; OpenV::OT is unsigned
                return OpenV(static_cast<typename OpenV::OT>(result.get()));
            }
        }
        OpenMFV mfv(mf.two_pow(static_cast<MFT>(exponent)));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V pow(V base, T exponent) const override
    {
        HPBC_CLOCKWORK_PRECONDITION2(0 <= exponent);
        if (ut_numeric_limits<T>::max() > ut_numeric_limits<MFT>::max()) {
            // kind of an unavoidable hack, so that AbstractMontgomeryForm has
            // the same contract for pow() as MongomeryForm, which allows
            // exponent to have any value of T >= 0.
            static constexpr MFT mft_max = ut_numeric_limits<MFT>::max();
            if (exponent > mft_max) {
                MFV mfv_base = OpenMFV(OpenV(base));
                MFV maxpow = mf.pow(mfv_base, mft_max);
                MFV accum = mf.getUnityValue();
                do {
                    accum = mf.multiply(accum, maxpow);
                    exponent -= mft_max;
                } while (exponent > static_cast<T>(mft_max));
                accum = mf.multiply(accum, mf.pow(mfv_base, static_cast<MFT>(exponent)));
                OpenMFV result(accum);
                // note: result.get() might be signed or unsigned; OpenV::OT is unsigned
                return OpenV(static_cast<typename OpenV::OT>(result.get()));
            }
        }
        OpenMFV mfv(mf.pow(OpenMFV(OpenV(base)), static_cast<MFT>(exponent)));
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

private:
    virtual V convertIn(T a, bool useLowlatencyTag) const override
    {
        HPBC_CLOCKWORK_PRECONDITION2(0 <= a);
        if (ut_numeric_limits<T>::max() > ut_numeric_limits<MFT>::max()) {
            // kind of an unavoidable hack, so that AbstractMontgomeryForm has
            // the same contract for convertIn() as MongomeryForm, which allows
            // 'a' to have any value of T >= 0.
            if (a > ut_numeric_limits<MFT>::max())
                a = a % getModulus();
        }
        OpenMFV mfv;
        if (useLowlatencyTag)
            mfv = OpenMFV(mf.template convertIn<LowlatencyTag>(static_cast<MFT>(a)));
        else
            mfv = OpenMFV(mf.template convertIn<LowuopsTag>(static_cast<MFT>(a)));
        // note that mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual T convertOut(V x, bool useLowlatencyTag) const override
    {
        T a;
        if (useLowlatencyTag)
            a = mf.template convertOut<LowlatencyTag>(OpenMFV(OpenV(x)));
        else
            a = mf.template convertOut<LowuopsTag>(OpenMFV(OpenV(x)));
        return a;
    }

    virtual V subtract(V x, V y, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template subtract<LowlatencyTag>(OpenMFV(OpenV(x)), OpenMFV(OpenV(y))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template subtract<LowuopsTag>(OpenMFV(OpenV(x)), OpenMFV(OpenV(y))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V subtract(V x, C y, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template subtract<LowlatencyTag>(OpenMFV(OpenV(x)), OpenMFC(OpenC(y))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template subtract<LowuopsTag>(OpenMFV(OpenV(x)), OpenMFC(OpenC(y))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V subtract(C x, V y, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template subtract<LowlatencyTag>(OpenMFC(OpenC(x)), OpenMFV(OpenV(y))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template subtract<LowuopsTag>(OpenMFC(OpenC(x)), OpenMFV(OpenV(y))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual C subtract(C x, C y, bool useLowlatencyTag) const override
    {
        OpenMFC mfc;
        if (useLowlatencyTag) {
            OpenMFC mfc2(mf.template subtract<LowlatencyTag>(OpenMFC(OpenC(x)), OpenMFC(OpenC(y))));
            mfc = mfc2;
        } else {
            OpenMFC mfc2(mf.template subtract<LowuopsTag>(OpenMFC(OpenC(x)), OpenMFC(OpenC(y))));
            mfc = mfc2;
        }
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual V multiply2(V x, V y, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template multiply<LowlatencyTag>(OpenMFV(OpenV(x)), OpenMFV(OpenV(y))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template multiply<LowuopsTag>(OpenMFV(OpenV(x)), OpenMFV(OpenV(y))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V multiply2(V x, V y, bool& resultIsZero, bool useLowlatencyTag)
        const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template multiply<LowlatencyTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), resultIsZero));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template multiply<LowuopsTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), resultIsZero));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V fmsub(V x, V y, C z, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template fmsub<LowlatencyTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), OpenMFC(OpenC(z))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template fmsub<LowuopsTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), OpenMFC(OpenC(z))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V fmsub(V x, V y, FV z, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template fmsub<LowlatencyTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), OpenMFFV(OpenFV(z))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template fmsub<LowuopsTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), OpenMFFV(OpenFV(z))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V fmadd(V x, V y, C z, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template fmadd<LowlatencyTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), OpenMFC(OpenC(z))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template fmadd<LowuopsTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), OpenMFC(OpenC(z))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V fmadd(V x, V y, FV z, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template fmadd<LowlatencyTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), OpenMFFV(OpenFV(z))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template fmadd<LowuopsTag>(OpenMFV(OpenV(x)),
                OpenMFV(OpenV(y)), OpenMFFV(OpenFV(z))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V square(V x, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template square<LowlatencyTag>(OpenMFV(OpenV(x))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template square<LowuopsTag>(OpenMFV(OpenV(x))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V fusedSquareSub(V x, C cv, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template fusedSquareSub<LowlatencyTag>(OpenMFV(OpenV(x)),
                OpenMFC(OpenC(cv))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template fusedSquareSub<LowuopsTag>(OpenMFV(OpenV(x)),
                OpenMFC(OpenC(cv))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual V fusedSquareAdd(V x, C cv, bool useLowlatencyTag) const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template fusedSquareAdd<LowlatencyTag>(OpenMFV(OpenV(x)),
                OpenMFC(OpenC(cv))));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template fusedSquareAdd<LowuopsTag>(OpenMFV(OpenV(x)),
                OpenMFC(OpenC(cv))));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    virtual C inverse(V x, bool useLowlatencyTag) const override
    {
        OpenMFC mfc;
        if (useLowlatencyTag) {
            OpenMFC mfc2(mf.template inverse<LowlatencyTag>(OpenMFV(OpenV(x))));
            mfc = mfc2;
        } else {
            OpenMFC mfc2(mf.template inverse<LowuopsTag>(OpenMFV(OpenV(x))));
            mfc = mfc2;
        }
        // note: mfc.get() might be signed or unsigned; OpenC::OT is unsigned
        return OpenC(static_cast<typename OpenC::OT>(mfc.get()));
    }

    virtual T remainder(T a, bool useLowlatencyTag) const override
    {
        HPBC_CLOCKWORK_PRECONDITION2(0 <= a);
        if (ut_numeric_limits<T>::max() > ut_numeric_limits<MFT>::max()) {
            // kind of an unavoidable hack, so that AbstractMontgomeryForm has
            // the same contract for remainder() as MongomeryForm, which allows
            // 'a' to have any value of T >= 0.
            if (a > ut_numeric_limits<MFT>::max())
                return a % getModulus();
        }
        T result;
        if (useLowlatencyTag)
            result = mf.template remainder<LowlatencyTag>(static_cast<MFT>(a));
        else
            result = mf.template remainder<LowuopsTag>(static_cast<MFT>(a));
        return result;
    }

    virtual V divideBySmallPowerOf2(C cx, int exponent, bool useLowlatencyTag)
        const override
    {
        OpenMFV mfv;
        if (useLowlatencyTag) {
            OpenMFV mfv2(mf.template divideBySmallPowerOf2<LowlatencyTag>(
                                                 OpenMFC(OpenC(cx)), exponent));
            mfv = mfv2;
        } else {
            OpenMFV mfv2(mf.template divideBySmallPowerOf2<LowuopsTag>(
                                                 OpenMFC(OpenC(cx)), exponent));
            mfv = mfv2;
        }
        // note: mfv.get() might be signed or unsigned; OpenV::OT is unsigned
        return OpenV(static_cast<typename OpenV::OT>(mfv.get()));
    }

    // This class (ConcreteMontgomeryForm) only supports calling vectorPow()
    // with a std::vector that has size equal to one of the sizes given by the
    // POW_ARRAY_SIZES variadic template argument for this class. Using a size
    // that was not specified is a programmer error. Since the only intended
    // use case of ConcreteMontgomeryForm is for unit testing, we can expect
    // the client unit test code will ensure the vector size of any call to
    // vectorPow() matches a (variadic argument) size that the client provided
    // when constructing this class.
    //
    // Note: the purpose of vectorPow() is to work around MontgomeryForm's use
    // of a templated array size in its pow() function - given that a virtual
    // function like vectorPow can't be templated - by having vectorPow take a
    // std::vector instead of a std::array. VectorPowHelper redirects the call
    // to MongomeryForm's pow funtion, assuming (and asserting) that the vector
    // size matches one of ConcreteMontgomeryForm's variadic argument sizes.
    // It's a clunky work-around, but it's probably usable (though ugly) in the
    // controlled environment of unit testing. It seems unlikely that any good
    // solution can exist to the problem of calling MongomeryForm's templated
    // pow function (given that pow should not be changed) via a virtual
    // function. Ultimately we're trying to unit test MongomeryForm, so it's
    // not an option to emulate the templated pow function - it has to be
    // called to be tested.
    template <std::size_t...> struct VectorPowHelper;

    // an adapter we use to call array pow, when we know at compile time
    // the exact size that the vector 'bases' will have.
    template <std::size_t A>
    static void fixed_size_vector_pow(const MF& mf,
        const std::vector<V>& bases, T exponent, std::vector<V>& answers)
    {
        HPBC_CLOCKWORK_ASSERT(bases.size() == A);   // if this fails, it is because
        // ConcreteMontgomeryForm was constructed with template variadic size_t
        // argument POW_ARRAY_SIZES that did not include the size of bases (as
        // used in this run-time assertion).
        // Most likely some code called AbstractMontgomeryWrapper or
        // AbstractMontgomeryForm's templated pow function, using for the pow
        // template argument (implicit or explicit) an array size that was not
        // included in this class's POW_ARRAY_SIZES.  The pow() call causes
        // this function to run, and a missing size in the values used for
        // POW_ARRAY_SIZES would cause the assertion above to fail.
        // If this is failing in unit testing, (at the time of this writing)
        // likely culprits for having a mismatch are test_MontgomeryForm.h
        // and test_MontgomeryForm_extra.cpp.

        std::array<MFV, A> arr;
        for (std::size_t i=0; i<A; ++i)
            arr[i] = OpenMFV(OpenV(bases[i]));

        std::array<MFV, A> result;
        if (ut_numeric_limits<T>::max() > ut_numeric_limits<MFT>::max()) {
            // kind of an unavoidable hack, so that AbstractMontgomeryForm has
            // the same contract for pow() as MongomeryForm, which allows
            // exponent to have any value of T >= 0.
            static constexpr MFT mft_max = ut_numeric_limits<MFT>::max();
            if (exponent > mft_max) {
                std::array<MFV, A> maxpow = mf.pow(arr, mft_max);
                result.fill(mf.getUnityValue());
                do {
                    for (std::size_t i=0; i<A; ++i)
                        result[i] = mf.multiply(result[i], maxpow[i]);
                    exponent -= mft_max;
                } while (exponent > static_cast<T>(mft_max));
                std::array<MFV, A> tmp = mf.pow(arr, static_cast<MFT>(exponent));
                for (std::size_t i=0; i<A; ++i)
                    result[i] = mf.multiply(tmp[i], result[i]);
            } else
                result = mf.pow(arr, static_cast<MFT>(exponent));
        } else
            result = mf.pow(arr, static_cast<MFT>(exponent));

        answers.resize(A);
        for (std::size_t i=0; i<A; ++i) {
            OpenMFV omfv(result[i]);
            // note: omfv.get() might be signed or unsigned; OpenV::OT is unsigned
            answers[i] = OpenV(static_cast<typename OpenV::OT>(omfv.get()));
        }
    }

    template <std::size_t A, std::size_t... B> struct VectorPowHelper<A, B...> {
        static void call(const MF& mf, const std::vector<V>& bases, T exponent, std::vector<V>& answers)
        {
            if (bases.size() == A)
                fixed_size_vector_pow<A>(mf, bases, exponent, answers);
            else
                VectorPowHelper<B...>::call(mf, bases, exponent, answers);
        }
    };
    template <std::size_t A> struct VectorPowHelper<A> {
        static void call(const MF& mf, const std::vector<V>& bases, T exponent, std::vector<V>& answers)
        {
            fixed_size_vector_pow<A>(mf, bases, exponent, answers);
        }
    };

    virtual std::vector<V> vectorPow(const std::vector<V>& bases, T exponent)
        const override
    {
        std::vector<V> answers;
        VectorPowHelper<POW_ARRAY_SIZES...>::call(mf, bases, exponent, answers);
        return answers;
    }


    struct EuclideanGcdFunctor {
        template <typename T1>
        T1 operator()(T1 a, T1 b) const
        {
            static_assert(ut_numeric_limits<T1>::is_integer, "");
            static_assert(!ut_numeric_limits<T1>::is_signed, "");
            HPBC_CLOCKWORK_PRECONDITION2(a > 0 || b > 0);
            while (a != 0) {
                T1 tmp = a;
                a = static_cast<T1>(b % a);
                b = tmp;
            }
            HPBC_CLOCKWORK_POSTCONDITION2(b > 0);
            return b;
        }
    };
    virtual T euclidean_gcd_with_modulus(V x) const override
    {
        return mf.gcd_with_modulus(OpenMFV(OpenV(x)), EuclideanGcdFunctor());
    }

};


} // end namespace

#endif
