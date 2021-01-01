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


// These are aliases for high performance MontgomeryForm types that are tailored
// to the size of the modulus you plan to use.  But unless you wish to squeeze
// out every possible performance advantage, you will generally find it more
// convenient to simply use MontgomeryForm<T>, instead of using these aliases.
//
// USAGE:
// The number in the Monty alias template denotes the size limit (by power of 2)
// of the modulus you can use to construct an object of the alias type.  For
// example, Monty<uint64_t, 63> allows any modulus < 2^63, and it would be an
// error to use any modulus >= 2^63 with this type.  Relatively speaking, an
// alias type with a smaller modulus limit will often perform better than an
// alias with a larger limit; you can expect that the smaller alias will never
// perform worse than the larger alias, when given the same modulus.
//
// PERFORMANCE DETAILS:
// Note that these aliases do not always improve upon the performance of
// MontgomeryForm<T> - the performance gain (or lack of gain) you can expect
// from these aliases depends on the particular alias:
// (1) The aliases Monty<uint32_t, 32>, Monty<uint64_t, 64>, and
// Monty<__uint128_t, 128> map to respectively the same classes as
// MontgomeryForm<uint32_t>, MontgomeryForm<uint64_t>, and
// MontgomeryForm<__uint128_t>.  Thus these aliases offer no improvement on
// performance over MontgomeryForm<T>.
// (2) The aliases Monty<uint32_t, 30>, Monty<uint32_t, 29>, Monty<uint64_t, 62>,
// Monty<uint64_t, 61>, Monty<__uint128_t, 126>, and Monty<__uint128_t, 125> 
// all offer significant performance improvements over their corresponding
// classes of MontgomeryForm<uint32_t>, MontgomeryForm<uint64_t>, and
// MontgomeryForm<__uint128_t>.  If you know either at compile-time or via run-
// time checking that your modulus will be small enough to allow you to use one
// of these aliases (see below), then you might expect these aliases to provide
// perhaps a 5-20% performance gain over MontgomeryForm<T>.  However if you use
// run-time checking, your checking is of course extra complexity and extra code
// size (since you need an alternate code path for when the modulus is not small
// enough).
// (3) The alias Monty<uint64_t, 61> improves upon Monty<uint64_t, 62> only for
// the famul() member function; otherwise it's the same as Monty<uint64_t, 62>.
// Likewise for Monty<uint32_t, 29> vs. Monty<uint32_t, 30>, and likewise for
// Monty<__uint128_t, 125> vs. Monty<__uint128_t, 126>.
// (4) The alias Monty<uint64_t, 63> improves upon Monty<uint64_t, 64> (and thus
// the equivalent MontgomeryForm<uint64_t>) only for the famul() member
// function; it's otherwise equivalent to both Monty<uint64_t, 64> and
// MontgomeryForm<uint64_t>.  Likewise for Monty<uint32_t, 31> vs.
// Monty<uint32_t, 32> (and MontgomeryForm<uint32_t>), and likewise for
// Monty<__uint128_t, 127> vs. Monty<__uint128_t, 128> (and
// MontygomeryForm<__uint128_t>).



// Don't directly use anything from the detail namespace here.  Instead, use the
// Monty template alias at the bottom of this file.
namespace detail {

  // primary template
  template<typename T, int BITS> struct MontyAlias {};

  // specializations:

  template<> struct MontyAlias<std::uint64_t, 64> { using type = 
      MontgomeryForm<std::uint64_t, detail::MontyFullRange<std::uint64_t>>;
  };
  template<> struct MontyAlias<std::uint64_t, 63> { using type = 
      MontgomeryForm<std::uint64_t, detail::MontyHalfRange<std::uint64_t>>;
  };
  template<> struct MontyAlias<std::uint64_t, 62> { using type = 
      MontgomeryForm<std::uint64_t, detail::MontyQuarterRange<std::uint64_t>>;
  };
  // technically, MontyAlias<61> works for any modulus < 2^(61.415)
  template<> struct MontyAlias<std::uint64_t, 61> { using type = 
      MontgomeryForm<std::uint64_t, detail::MontySixthRange<std::uint64_t>>;
  };

#ifndef HURCHALLA_TARGET_BIT_WIDTH
#   error "HURCHALLA_TARGET_BIT_WIDTH must be defined (normally compiler_macros.h automatically detects it)"
#endif
#if (HURCHALLA_TARGET_BIT_WIDTH < 64)
  template<> struct MontyAlias<std::uint32_t, 32> { using type = 
      MontgomeryForm<std::uint32_t, detail::MontyFullRange<std::uint32_t>>;
  };
  template<> struct MontyAlias<std::uint32_t, 31> { using type = 
      MontgomeryForm<std::uint32_t, detail::MontyHalfRange<std::uint32_t>>;
  };
  template<> struct MontyAlias<std::uint32_t, 30> { using type = 
      MontgomeryForm<std::uint32_t, detail::MontyQuarterRange<std::uint32_t>>;
  };
  // technically, MontyAlias<29> works for any modulus < 2^(29.415)
  template<> struct MontyAlias<std::uint32_t, 29> { using type = 
      MontgomeryForm<std::uint32_t, detail::MontySixthRange<std::uint32_t>>;
  };
#else
  template<> struct MontyAlias<std::uint32_t, 32> { using type = 
      MontgomeryForm<std::uint32_t, detail::MontySixthRange<std::uint64_t>>;
  };
  template<> struct MontyAlias<std::uint32_t, 31> { using type = 
      MontgomeryForm<std::uint32_t, detail::MontySixthRange<std::uint64_t>>;
  };
  template<> struct MontyAlias<std::uint32_t, 30> { using type = 
      MontgomeryForm<std::uint32_t, detail::MontySixthRange<std::uint64_t>>;
  };
  // technically, MontyAlias<29> works for any modulus < 2^(29.415)
  template<> struct MontyAlias<std::uint32_t, 29> { using type = 
      MontgomeryForm<std::uint32_t, detail::MontySixthRange<std::uint64_t>>;
  };
#endif

#if (HURCHALLA_COMPILER_HAS_UINT128_T())
  template<> struct MontyAlias<__uint128_t, 128> { using type = 
      MontgomeryForm<__uint128_t, detail::MontyFullRange<__uint128_t>>;
  };
  template<> struct MontyAlias<__uint128_t, 127> { using type = 
      MontgomeryForm<__uint128_t, detail::MontyHalfRange<__uint128_t>>;
  };
  template<> struct MontyAlias<__uint128_t, 126> { using type = 
      MontgomeryForm<__uint128_t, detail::MontyQuarterRange<__uint128_t>>;
  };
  // technically, MontyAlias<125> works for any modulus < 2^(125.415)
  template<> struct MontyAlias<__uint128_t, 125> { using type = 
      MontgomeryForm<__uint128_t, detail::MontySixthRange<__uint128_t>>;
  };
#endif

}  // end namespace detail


template <typename T, int BITS>
using Monty = typename detail::MontyAlias<T, BITS>::type;


}} // end namespace

#endif
