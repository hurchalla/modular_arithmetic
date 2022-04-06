// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "test_MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"
#include "hurchalla/montgomery_arithmetic/detail/experimental/MontyFullRangeMasked.h"
#include "gtest/gtest.h"


namespace {


// test the 'unusual' Montgomery types, which are MontyWrappedStandardMath and
// the experimental class MontyFullRangeMasked.

TEST(MontgomeryArithmetic, MontyWrappedStandardMath) {
    test_custom_monty<hurchalla::detail::MontyWrappedStandardMath>();
}

TEST(MontgomeryArithmetic, MontyFullRangeMasked) {
    test_custom_monty<hurchalla::detail::MontyFullRangeMasked>();
}


} // end anonymous namespace
