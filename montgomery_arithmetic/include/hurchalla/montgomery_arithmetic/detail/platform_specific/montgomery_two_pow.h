// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_TWO_POW_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_TWO_POW_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/impl_montgomery_two_pow.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyTags.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstddef>
#include <array>
#include <type_traits>

namespace hurchalla { namespace detail {


class Tag_montgomery_two_pow_gcc {};
class Tag_montgomery_two_pow_clang {};

// Essentially Tag_montgomery_two_pow_big is for __uint128_t (or bigger?)
class Tag_montgomery_two_pow_big {};
// Essentially Tag_montgomery_two_pow_small is for uint64_t or smaller
class Tag_montgomery_two_pow_small {};



// Primary template for tagging different perf tunings
template <class MontyTag, class CompilerTag, class PowSizeTag>
struct tagged_montgomery_two_pow {};




// These specializations call impl_montgomery_two_pow, using the template
// arguments that were found to perform best in benchmarks for a corresponding
// configuration of compiler/uint_type/Mont_type


// Fyi, the meaning of the impl_montgomery_two_pow::call template parameters are
// scalar call
//  template <class MF, typename U,
//            bool USE_SLIDING_WINDOW_OPTIMIZATION,
//            size_t TABLE_BITS, size_t CODE_SECTION,
//            bool USE_SQUARING_VALUE_OPTIMIZATION>
//
// array call
//  template <class MF, typename U,
//            size_t ARRAY_SIZE, size_t TABLE_BITS, size_t CODE_SECTION,
//            bool USE_SQUARING_VALUE_OPTIMIZATION>


// -- the following best performance tunings were measured with mac M2 --


// Partial Specialization: clang and big uint pow.
// Intended for MontgomeryFull, but catches all non-specialized monty types
template <class MontyTag> struct tagged_montgomery_two_pow
   <MontyTag, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 34, true>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 30, false>(mf, n);
  }
};
// Full Specialization: clang and big uint pow and MontgomeryHalf.
template <> struct tagged_montgomery_two_pow
   <TagMontyHalfrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 22, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: clang and big uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
   <TagMontyQuarterrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 22, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 30, false>(mf, n);
  }
};


// Partial specialization: gcc and big uint pow.
// Intended for MontgomeryFull, but catches all non-specialized monty types
template <class MontyTag> struct tagged_montgomery_two_pow
   <MontyTag, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 33, true>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 30, false>(mf, n);
  }
};
// Full Specialization: gcc and big uint pow and MontgomeryHalf.
template <> struct tagged_montgomery_two_pow
   <TagMontyHalfrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 33, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: gcc and big uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
   <TagMontyQuarterrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 22, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};



// Partial Specialization: clang and small uint pow.
// Intended for MontgomeryFull, but catches all non-specialized monty types
template <class MontyTag> struct tagged_montgomery_two_pow
   <MontyTag, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 33, true>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: clang and small uint pow and MontgomeryHalf.
template <> struct tagged_montgomery_two_pow
   <TagMontyHalfrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 24, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 29, false>(mf, n);
  }
};
// Full Specialization: clang and small uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
   <TagMontyQuarterrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 24, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};


// Partial Specialization: gcc and small uint pow and MontgomeryFull.
// Intended for MontgomeryFull, but catches all non-specialized monty types
template <class MontyTag> struct tagged_montgomery_two_pow
   <MontyTag, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 23, true>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 31, false>(mf, n);
  }
};
// Full Specialization: gcc and small uint pow and MontgomeryHalf.
template <> struct tagged_montgomery_two_pow
   <TagMontyHalfrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 24, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 28, false>(mf, n);
  }
};
// Full Specialization: gcc and small uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
   <TagMontyQuarterrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<MF, U, false, 0, 24, false>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<MF, U, ARRAY_SIZE, 0, 29, false>(mf, n);
  }
};




struct montgomery_two_pow {

  // Calculate pow(2, n), modulo the modulus of mf, and return the result in
  // montgomeryform representation.
  template <class MF, typename T>
  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, T nt)
  {
    HPBC_CLOCKWORK_PRECONDITION(nt >= 0);
    using U = typename extensible_make_unsigned<T>::type;
    U n = static_cast<U>(nt);

    using MontyTag = typename MF::MontType::MontyTag;
    using RU = typename MontgomeryFormExtensions<MF, LowlatencyTag>::RU;
    constexpr bool isBigPow = ut_numeric_limits<RU>::digits >
                              HURCHALLA_TARGET_BIT_WIDTH;
    using SizeTag = typename std::conditional<isBigPow,
                Tag_montgomery_two_pow_big, Tag_montgomery_two_pow_small>::type;
#if (!defined(__GNUC__) || defined(__clang__))
    using Compiler = Tag_montgomery_two_pow_clang;   // clang tuning
#else
    using Compiler = Tag_montgomery_two_pow_gcc;     // gcc tuning
#endif
    
    return tagged_montgomery_two_pow<MontyTag, Compiler, SizeTag>::call(mf, n);
  }


  // Helper function - delegated Array version of montgomery two pow
  template <class MF, typename U, std::size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  helper(const std::array<MF,ARRAY_SIZE>& mf, const std::array<U,ARRAY_SIZE>& n)
  {
    static_assert(hurchalla::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hurchalla::ut_numeric_limits<U>::is_signed, "");

    using MontyTag = typename MF::MontType::MontyTag;
    using RU = typename MontgomeryFormExtensions<MF, LowlatencyTag>::RU;
    constexpr bool isBigPow = ut_numeric_limits<RU>::digits >
                              HURCHALLA_TARGET_BIT_WIDTH;
    using SizeTag = typename std::conditional<isBigPow,
                Tag_montgomery_two_pow_big, Tag_montgomery_two_pow_small>::type;
#if (!defined(__GNUC__) || defined(__clang__))
    using Compiler = Tag_montgomery_two_pow_clang;   // clang tuning
#else
    using Compiler = Tag_montgomery_two_pow_gcc;     // gcc tuning
#endif

    return tagged_montgomery_two_pow<MontyTag, Compiler, SizeTag>::call(mf, n);
  }


  // Array version of montgomery two pow, for unsigned T
  template <class MF, typename T, std::size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE static
  typename std::enable_if<!(hurchalla::ut_numeric_limits<T>::is_signed),
                          std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
                          >::type
  call(const std::array<MF,ARRAY_SIZE>& mf, const std::array<T,ARRAY_SIZE>& nt)
  {
      return helper(mf, nt);
  }

  // Array version of montgomery two pow, for signed T
  template <class MF, typename T, std::size_t ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE static
  typename std::enable_if<(hurchalla::ut_numeric_limits<T>::is_signed),
                          std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
                          >::type
  call(const std::array<MF,ARRAY_SIZE>& mf, const std::array<T,ARRAY_SIZE>& nt)
  {
      using U = typename extensible_make_unsigned<T>::type;
      std::array<U, ARRAY_SIZE> n;
      HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<ARRAY_SIZE; ++i) {
        HPBC_CLOCKWORK_PRECONDITION(nt[i] >= 0);
        n[i] = static_cast<U>(nt[i]);
      }
      return helper(mf, n);
  }
};


}} // end namespace

#endif
