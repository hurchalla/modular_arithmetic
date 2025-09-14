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

    using RU = typename MF::MontType::uint_type;
    // conceptually, R = 1 << (ut_numeric_limits<RU>::digits)
    static_assert(ut_numeric_limits<RU>::is_integer, "");
    static_assert(!(ut_numeric_limits<RU>::is_signed), "");

    using CanonicalValue = typename MF::CanonicalValue;
    using MontgomeryValue = typename MF::MontgomeryValue;
    using SquaringValue = typename MF::MontType::squaringvalue_type;

    HURCHALLA_FORCE_INLINE
    static MontgomeryValue convertInExtended(const MF& mf, RU a)
    {
        return mf.impl.template convertInExtended<PTAG>(a);
    }

    // note: montvalueR is the Montgomery representation of R.
    //       In normal integer form it is literally R squared mod N.
    HURCHALLA_FORCE_INLINE
    static CanonicalValue getMontvalueR(const MF& mf)
    {
        return mf.impl.getMontvalueR();
    }

    // this first shifts x by exponent, which is equivalent to
    // multiplying x by 2^exponent, and then it completes the
    // mont mul as usual by calling REDC.
    // -- IMPORTANT NOTE -- because (2^exponent) is an integer domain
    // value rather than a montgomery domain value, the returned
    // result viewed as an integer value is
    // REDC((x_int * R) * (2^exponent)) == (x_int * (2^exponent) * R) * R^(-1)
    // To counteract the inverse R factor, so that you get what most likely
    // you wanted, being just plain   (x_int * (2^exponent) * R),
    // you need to ensure that x has an extra factor of R built into it it,
    // rather than just the normal single factor of x_int * R.  To build an
    // extra factor of R into x, you first get  montR = getMontvalueR(mf),
    // and then you do a normal montgomery multiply of x and montR.
    HURCHALLA_FORCE_INLINE
    static MontgomeryValue twoPowLimited_times_x(const MF& mf, size_t exponent, CanonicalValue x)
    {
        HPBC_CLOCKWORK_PRECONDITION(exponent < ut_numeric_limits<RU>::digits);
        return mf.impl.template twoPowLimited_times_x<PTAG>(exponent, x);
    }
    HURCHALLA_FORCE_INLINE
    static MontgomeryValue twoPowLimited_times_x_v2(const MF& mf, size_t exponent, CanonicalValue x)
    {
        HPBC_CLOCKWORK_PRECONDITION(0 < exponent && exponent <= ut_numeric_limits<RU>::digits);
        return mf.impl.template twoPowLimited_times_x_v2<PTAG>(exponent, x);
    }

    // note: magicValue is R cubed mod N  (in normal integer form)
    HURCHALLA_FORCE_INLINE
    static RU getMagicValue(const MF& mf)
    {
        return mf.impl.template getMagicValue<PTAG>();
    }

    HURCHALLA_FORCE_INLINE
    static MontgomeryValue
    convertInExtended_aTimesR(const MF& mf, RU a, RU magicValue)
    {
        HPBC_CLOCKWORK_PRECONDITION(magicValue == getMagicValue(mf));
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
        HPBC_CLOCKWORK_PRECONDITION(exponent < ut_numeric_limits<RU>::digits);
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
        HPBC_CLOCKWORK_PRECONDITION(exponent < ut_numeric_limits<RU>::digits);
        return mf.impl.template RTimesTwoPowLimited<PTAG>(exponent, magicValue);
    }


    HURCHALLA_FORCE_INLINE
    static SquaringValue getSquaringValue(const MF& mf, MontgomeryValue x)
    {
        return mf.impl.getSquaringValue(x);
    }

    HURCHALLA_FORCE_INLINE
    static SquaringValue squareSV(const MF& mf, SquaringValue sv)
    {
        return mf.impl.template squareSV<PTAG>(sv);
    }

    HURCHALLA_FORCE_INLINE
    static MontgomeryValue
    squareToMontgomeryValue(const MF& mf, SquaringValue sv)
    {
        return mf.impl.template squareToMontgomeryValue<PTAG>(sv);
    }

    HURCHALLA_FORCE_INLINE
    static MontgomeryValue getMontgomeryValue(const MF& mf, SquaringValue sv)
    {
        return mf.impl.getMontgomeryValue(sv);
    }
};


}} // end namespace

#endif
