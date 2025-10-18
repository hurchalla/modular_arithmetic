// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_EXPERIMENTAL_API_MONTGOMERY_TWO_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_EXPERIMENTAL_API_MONTGOMERY_TWO_POW_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/montgomery_two_pow.h"
#include <cstddef>
#include <array>

namespace experimental_api {


// Calculates the integer pow(2, n), modulo the modulus of mf, and returns the
// result in MongomeryForm representation.
//
// MF can be any MontgomeryForm type (see MontgomeryForm.h), and U can be any
// integer type.  ('n' is the exponent to use)
//
template <class MF, typename U>
typename MF::MontgomeryValue montgomery_two_pow(const MF& mf, U n)
{
    // Rather than calling this function, you could just directly call
    // mf.two_pow(n), as done in the next line.
    return mf.two_pow(n);

    // Implementation note: the above function call internally just delegates to
    // return hurchalla::detail::montgomery_two_pow::call(mf, n);
    // It uses novel optimizations of the k-ary exponentiation algorithm
    // ( https://en.wikipedia.org/wiki/Exponentiation_by_squaring )
    // that rely on a hard-coded base 2.
}


// An array version of the above function - you can expect it to always have
// significantly higher throughput than the above.  (In benchmarks I have
// observed it to have a performance advantage of anywhere from 1.4x to 3x
// higher throughput, depending on the CPU type and whether 64 or 128 bit
// integer types are calculated)
//
// For each array index 'i' from 0 to ARRAY_SIZE-1, this function calculates
// the integer result[i] = pow(2, n[i])  modulo the modulus of mf[i], and
// returns this result array; the result array is in MontgomeryForm
// representation.
//
// MF can be any MontgomeryForm type (see MontgomeryForm.h), and U can be any
// integer type.
//
template <class MF, typename U, size_t ARRAY_SIZE>
std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
array_montgomery_two_pow(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
{
    // Implementation note: at the moment this API function is the only easy way
    // to get the array version of Montgomery two pow (MontgomeryForm.h does not
    // have an *array* two_pow member function).
    // At some point in the next 8 months I expect to create a SIMD version of
    // MontgomeryForm, and at that time the SIMD MontgomeryForm will become the
    // preferred API to use to access the (high throughput) array version of
    // Montgomery two_pow.

    return hurchalla::detail::montgomery_two_pow::call(mf, n);
}


} // end namespace

#endif
