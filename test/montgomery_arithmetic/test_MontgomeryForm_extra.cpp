// Copyright (c) 2020-2024 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "test_MontgomeryForm.h"
#include "NoForceInlineMontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/MontyFullRangeMasked.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/AbstractMontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/ConcreteMontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/AbstractMontgomeryWrapper.h"
#include "gtest/gtest.h"


namespace {


#if 0
 template <class T, class Monty> using MF =
    hurchalla::MontgomeryForm<T, Monty>;
#else
 template<class T, class Monty> using MF =
    hurchalla::NoForceInlineMontgomeryForm<hurchalla::MontgomeryForm<T, Monty>>;
#endif

// test the 'unusual' Montgomery types, which are MontyWrappedStandardMath and
// the experimental class MontyFullRangeMasked.

TEST(MontgomeryArithmetic, MontyWrappedStandardMath) {
    test_custom_monty<MF, hurchalla::detail::MontyWrappedStandardMath>();
}


#ifdef HURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT
// MontyFullRangeMasked is experimental, so we skip it when we're not doing
// extensive (heavyweight) testing.
TEST(MontgomeryArithmetic, MontyFullRangeMasked) {
    test_custom_monty<MF, hurchalla::detail::MontyFullRangeMasked>();
}

// The group of classes: ConcreteMontgomeryForm, AbstractMontgomeryForm, and
// AbstractMontgomeryWrapper, are experimental, so we skip testing them when
// we're not doing extensive (heavyweight) testing.
TEST(MontgomeryArithmetic, MontyVirtual) {
    {
        using ConcreteMF = hurchalla::ConcreteMontgomeryForm<hurchalla::MontgomeryForm<__uint128_t>,
                                                             TESTABLE_ARRAY_POW_SIZES()>;
        using Wrapper = hurchalla::AbstractMontgomeryWrapper<ConcreteMF::Parent>;
        test_MontgomeryForm<Wrapper, ConcreteMF>();
    }
    {
        using ConcreteMF = hurchalla::ConcreteMontgomeryForm<hurchalla::MontgomeryForm<uint32_t>,
                                                             TESTABLE_ARRAY_POW_SIZES()>;
        using Wrapper = hurchalla::AbstractMontgomeryWrapper<ConcreteMF::Parent>;
        test_MontgomeryForm<Wrapper, ConcreteMF>();
    }
    {
        using ConcreteMF = hurchalla::ConcreteMontgomeryForm<hurchalla::MontgomeryForm<int32_t>,
                                                             TESTABLE_ARRAY_POW_SIZES()>;
        using Wrapper = hurchalla::AbstractMontgomeryWrapper<ConcreteMF::Parent>;
        test_MontgomeryForm<Wrapper, ConcreteMF>();
    }
}
#endif


} // end anonymous namespace
