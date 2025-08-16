// Copyright (c) 2020-2024 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "test_MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/MontyFullRangeMasked.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/AbstractMontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/ConcreteMontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/unit_testing_helpers/AbstractMontgomeryWrapper.h"
#include "hurchalla/util/compiler_macros.h"
#include "gtest/gtest.h"


namespace {


// For unit testing, we want fast compile times, so it helps to use the version
// of MontgomeryForm that generally doesn't do force inlining.
#if 1
constexpr bool forceInlineAllFunctions = false;
#else
// note: even the default template arg for MontgomeryForm wouldn't have force
// inlined everything for uint128_t or int128_t (we would expect the functions
// to have too many instructions for it to be a good idea).  So for some T we
// get more inlining than the default, when this #else is enabled.
constexpr bool forceInlineAllFunctions = true;
#endif

template <class T, class Monty> using MF =
    hurchalla::MontgomeryForm<T, forceInlineAllFunctions, Monty>;



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
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    {
        using ConcreteMF = hurchalla::ConcreteMontgomeryForm<hurchalla::MontgomeryForm<__uint128_t>,
                                                             TESTABLE_ARRAY_POW_SIZES()>;
        using Wrapper = hurchalla::AbstractMontgomeryWrapper<ConcreteMF::Parent>;
        test_MontgomeryForm<Wrapper, ConcreteMF>();
    }
#endif
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
