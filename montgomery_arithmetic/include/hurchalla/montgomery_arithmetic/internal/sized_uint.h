
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_SIZED_UINT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_SIZED_UINT_H_INCLUDED


#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


template <int BITS> struct sized_uint { using type = void; };
template<> struct sized_uint<8>  { using type = uint8_t; };
template<> struct sized_uint<16> { using type = uint16_t; };
template<> struct sized_uint<32> { using type = uint32_t; };
template<> struct sized_uint<64> { using type = uint64_t; };
// uint128_t would be nice, but it's unsupported on most compilers
//template<> struct sized_uint<128> { using type = uint128_t; 
//                                   static constexpr bool is_valid = true; };


}} // end namespace

#endif
