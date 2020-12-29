// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_ALIASES_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_FORM_ALIASES_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySixthRange.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


// These are convenient aliases for some high performance MontgomeryForm types.
// For usage, see USAGE below.
// In general though, we should note that rather than these aliases, the most
// convenient name to use is simply MontgomeryForm<T>.  And we should further
// note that these aliases do not always improve upon the performance of
// MontgomeryForm<T> - the performance gain (or lack of gain) you can expect
// from these aliases depends on the particular alias:
// (1) The aliases Monty64, Monty32, and Monty128 map to the same classes as
// MontgomeryForm<uint64_t>, MontgomeryForm<uint32_t>, and
// MontgomeryForm<__uint128_t>.  Thus these aliases offer no improvement on
// performance over MontgomeryForm<T>.
// (2) The aliases Monty62/61, Monty30/29, and Monty126/125 do offer significant
// performance improvements over, respectively, MontgomeryForm<uint64_t>,
// MontgomeryForm<uint32_t>, and MontgomeryForm<__uint128_t>.  If you know
// either at compile-time or via run-time checking that your modulus will be
// small enough to allow you to use one of these aliases (see below), then you
// might expect perhaps a 5-20% performance improvement.  However run-time
// checking adds complexity and extra code size (since you need an alternate
// code path for when the modulus is not small enough).
// (3) The alias Monty61 improves upon Monty62 only for the famul() member
// function; it's otherwise equivalent to Monty62.  Likewise for Monty29 vs.
// Monty30, and likewise for Monty125 vs. Monty126.
// (4) The alias Monty63 improves upon Monty64 (and thus the equivalent
// MontgomeryForm<uint64_t>) only for the famul() member function; it's
// otherwise equivalent to Monty64 and MontgomeryForm<uint64_t>.  Likewise for
// Monty31 vs. Monty32/MontgomeryForm<uint32_t>, and likewise for Monty127 vs.
// Monty128/MontygomeryForm<__uint128_t>.

// USAGE:
// The number in the alias name denotes the size limit (by power of 2) of the
// modulus you can use to construct an object of the alias type.  For example,
// Monty63 allows any modulus < 2^63, and it is an error to use any larger size
// modulus with this type.  An alias type with a smaller modulus limit will
// often perform better than an alias with a larger limit; the smaller alias
// will never perform worse than the larger alias, when given the same modulus.


// -- When your type is uint64_t --
using Monty64 = MontgomeryForm<
                       std::uint64_t, detail::MontyFullRange<std::uint64_t>>;
using Monty63 = MontgomeryForm<
                       std::uint64_t, detail::MontyHalfRange<std::uint64_t>>;
using Monty62 = MontgomeryForm<
                       std::uint64_t, detail::MontyQuarterRange<std::uint64_t>>;
// technically, Monty61 works for any modulus < 2^(61.415)
using Monty61 = MontgomeryForm<
                       std::uint64_t, detail::MontySixthRange<std::uint64_t>>;


// -- When your type is uint32_t --
#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  error "HURCHALLA_TARGET_BIT_WIDTH must be defined (normally compiler_macros.h automatically detects it)"
#endif
#if (HURCHALLA_TARGET_BIT_WIDTH < 64)
  using Monty32 = MontgomeryForm<
                       std::uint32_t, detail::MontyFullRange<std::uint32_t>>;
  using Monty31 = MontgomeryForm<
                       std::uint32_t, detail::MontyHalfRange<std::uint32_t>>;
  using Monty30 = MontgomeryForm<
                       std::uint32_t, detail::MontyQuarterRange<std::uint32_t>>;
  // technically, Monty29 works for any modulus < 2^(29.415)
  using Monty29 = MontgomeryForm<
                       std::uint32_t, detail::MontySixthRange<std::uint32_t>>;
#else
  using Monty32 = MontgomeryForm<
                       std::uint32_t, detail::MontySixthRange<std::uint64_t>>;
  using Monty31 = MontgomeryForm<
                       std::uint32_t, detail::MontySixthRange<std::uint64_t>>;
  using Monty30 = MontgomeryForm<
                       std::uint32_t, detail::MontySixthRange<std::uint64_t>>;
  using Monty29 = MontgomeryForm<
                       std::uint32_t, detail::MontySixthRange<std::uint64_t>>;
#endif


// -- When your type is __uint128_t --
#if (HURCHALLA_COMPILER_HAS_UINT128_T())
  using Monty128 = MontgomeryForm<
                 std::__uint128_t, detail::MontyFullRange<std::__uint128_t>>;
  using Monty127 = MontgomeryForm<
                 std::__uint128_t, detail::MontyHalfRange<std::__uint128_t>>;
  using Monty126 = MontgomeryForm<
                 std::__uint128_t, detail::MontyQuarterRange<std::__uint128_t>>;
  // technically, Monty125 works for any modulus < 2^(125.415)
  using Monty125 = MontgomeryForm<
                 std::__uint128_t, detail::MontySixthRange<std::__uint128_t>>;
#endif


}} // end namespace

#endif
