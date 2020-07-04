// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_SIZED_UINT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_SIZED_UINT_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


template <int BITS> struct sized_uint { using type = void; };
template<> struct sized_uint<8>   { using type = uint8_t; };
template<> struct sized_uint<16>  { using type = uint16_t; };
template<> struct sized_uint<32>  { using type = uint32_t; };
template<> struct sized_uint<64>  { using type = uint64_t; };

#if (HURCHALLA_COMPILER_HAS_UINT128_T())
  // Some compilers support __uint128_t, so we'll specialize with it if possible
  template<> struct sized_uint<128> { using type = __uint128_t; };
#elif (HURCHALLA_TARGET_BIT_WIDTH > = 128)
  // presumably if a CPU target is >= 128bit, type uint128_t will exist.
  template<> struct sized_uint<128> { using type = uint128_t; };
#endif


}} // end namespace

#endif
