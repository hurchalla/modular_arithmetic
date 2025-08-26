// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_EXTENSIONS_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_EXTENSIONS_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <cstddef>

namespace hurchalla { namespace detail { 


// Implementation helper functions that shouldn't be exposed in the
// MontgomeryForm API.


template <class MF, class PTAG>
struct MontgomeryFormExtensions final {

    using RU = typename MF::RU;
    // conceptually, R = 1 << (ut_numeric_limits<RU>::digits)
    static_assert(ut_numeric_limits<RU>::is_integer, "");
    static_assert(!(ut_numeric_limits<RU>::is_signed), "");

    using MontgomeryValue = typename MF::MontgomeryValue;

    HURCHALLA_FORCE_INLINE
    static MontgomeryValue convertInExtended(const MF& mf, RU a)
    {
        return mf.impl.template convertInExtended<PTAG>(a);
    }

    // note: magicValue is R cubed mod N  (in normal integer form)
    HURCHALLA_FORCE_INLINE
    static RU getMagicValue(const MF& mf)
    {
        return mf.impl.getMagicValue();
    }

    HURCHALLA_FORCE_INLINE
    static MontgomeryValue
    convertInExtended_aTimesR(const MF& mf, RU a, RU magicValue)
    {
        HPBC_CLOCKWORK_PRECONDITION1(magicValue == getMagicValue(mf));
        return mf.impl.template convertInExtended_aTimesR<PTAG>(a, magicValue);
    }

    // this shifts RsquaredModN by exponent (rather than multiplying by
    // (1<<exponent)) before calling REDC as usual.
    // The amount RsquaredModN can be shifted is limited by the bit width of
    // RsquaredModN's type - shifting more would be undefined behavior.
    // Thus the (exponent) shift is limited to 0 <= shift < digitsR.
    HURCHALLA_FORCE_INLINE
    static MontgomeryValue twoPowLimited(const MF& mf, size_t exponent)
    {
        HPBC_CLOCKWORK_PRECONDITION1(exponent < ut_numeric_limits<RU>::digits);
        return mf.impl.template twoPowLimited<PTAG>(exponent);
    }

    // this shifts RcubedModN by exponent (rather than multiplying by
    // (1<<exponent)) before calling REDC as usual.
    // Similarly to twoPowLimited, the exponent shift must be limited
    // to 0 <= shift < digitsR.
    HURCHALLA_FORCE_INLINE
    static MontgomeryValue
    RTimesTwoPowLimited(const MF& mf, size_t exponent, RU magicValue)
    {
        HPBC_CLOCKWORK_PRECONDITION1(exponent < ut_numeric_limits<RU>::digits);
        return mf.impl.template RTimesTwoPowLimited<PTAG>(exponent, magicValue);
    }
};


}} // end namespace

#endif
