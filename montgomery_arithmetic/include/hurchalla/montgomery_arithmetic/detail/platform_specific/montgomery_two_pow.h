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


// -- Generic tunings, measured with mac M2 with no asm enabled --

// Full Specialization: clang and big uint pow and MontgomeryFull.
template <> struct tagged_montgomery_two_pow
   <TagMontyFullrange, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<true, 0, 3, MF, U>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<0, 0, MF, U, ARRAY_SIZE>(mf, n);
  }
};
// Partial specialization: clang and big uint pow.
template <class MontyTag> struct tagged_montgomery_two_pow
       <MontyTag, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<true, 0, 3, MF, U>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<0, 2, MF, U, ARRAY_SIZE>(mf, n);
  }
};
// Partial specialization: clang and small uint pow.
template <class MontyTag> struct tagged_montgomery_two_pow
       <MontyTag, Tag_montgomery_two_pow_clang, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<true, 0, 1, MF, U>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<0, 0, MF, U, ARRAY_SIZE>(mf, n);
  }
};


// Full specialization: gcc and big uint pow and MontgomeryQuarter.
template <> struct tagged_montgomery_two_pow
  <TagMontyQuarterrange, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<true, 0, 1, MF, U>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<0, 0, MF, U, ARRAY_SIZE>(mf, n);
  }
};
// Partial specialization: gcc and big uint pow.
template <class MontyTag> struct tagged_montgomery_two_pow
       <MontyTag, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_big>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<false, 0, 2, MF, U>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<0, 0, MF, U, ARRAY_SIZE>(mf, n);
  }
};
// Partial specialization: gcc and small uint pow.
template <class MontyTag> struct tagged_montgomery_two_pow
       <MontyTag, Tag_montgomery_two_pow_gcc, Tag_montgomery_two_pow_small>
{
  template <class MF, typename U>  HURCHALLA_FORCE_INLINE
  static typename MF::MontgomeryValue call(const MF& mf, U n)
  {
    return impl_montgomery_two_pow::call<true, 0, 3, MF, U>(mf, n);
  }
  template <class MF, typename U, std::size_t ARRAY_SIZE> HURCHALLA_FORCE_INLINE
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  call(const std::array<MF, ARRAY_SIZE>& mf, const std::array<U, ARRAY_SIZE>& n)
  {
    return impl_montgomery_two_pow::call<0, 0, MF, U, ARRAY_SIZE>(mf, n);
  }
};



struct montgomery_two_pow {

  // Calculate pow(2, n), modulo the modulus of mf, and return the result in
  // montgomeryform representation.
  template <class MF, typename T>
  static typename MF::MontgomeryValue call(const MF& mf, T nt)
  {
    HPBC_PRECONDITION(nt >= 0);
    using U = typename extensible_make_unsigned<T>::type;
    U n = static_cast<U>(nt);

    using MontyTag = typename MF::MontyTag;
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
  static std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
  HURCHALLA_FORCE_INLINE
  helper(const std::array<MF,ARRAY_SIZE>& mf, const std::array<U,ARRAY_SIZE>& n)
  {
    static_assert(hurchalla::ut_numeric_limits<U>::is_integer, "");
    static_assert(!hurchalla::ut_numeric_limits<U>::is_signed, "");

    using MontyTag = typename MF::MontyTag;
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
  static
  typename std::enable_if<!(hurchalla::ut_numeric_limits<T>::is_signed),
                          std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
                          >::type
  call(const std::array<MF,ARRAY_SIZE>& mf, const std::array<T,ARRAY_SIZE>& nt)
  {
      return helper<MF, T, ARRAY_SIZE>(mf, nt);
  }

  // Array version of montgomery two pow, for signed T
  template <class MF, typename T, std::size_t ARRAY_SIZE>
  static
  typename std::enable_if<(hurchalla::ut_numeric_limits<T>::is_signed),
                          std::array<typename MF::MontgomeryValue, ARRAY_SIZE>
                          >::type
  call(const std::array<MF,ARRAY_SIZE>& mf, const std::array<T,ARRAY_SIZE>& nt)
  {
      using U = typename extensible_make_unsigned<T>::type;
      std::array<U, ARRAY_SIZE> n;
      HURCHALLA_REQUEST_UNROLL_LOOP for (int i=0; i<ARRAY_SIZE; ++i) {
        HPBC_PRECONDITION(nt[i] >= 0);
        n[i] = static_cast<U>(nt[i]);
      }
      return helper<MF, U, ARRAY_SIZE>(mf, n);
  }
};


}} // end namespace

#endif
