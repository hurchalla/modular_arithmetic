// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/impl_array_get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/unsigned_square_to_hilo_product.h"
#include "hurchalla/util/count_leading_zeros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"

#include "experimental_montgomery_pow_2kary.h"

#include <iostream>
#include <stdexcept>
#include <chrono>
#include <utility>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <cstdlib>


#if defined(HURCHALLA_CLOCKWORK_ENABLE_ASSERTS) || defined(HURCHALLA_UTIL_ENABLE_ASSERTS)
#  warning "asserts are enabled and will slow performance"
#endif



// this utility function vector_to_stdarray() is adapted from
// https://stackoverflow.com/a/18497366
template<std::size_t ... N>
struct seq
{
   using type = seq<N...>;
   template <std::size_t I> struct push_back : seq<N..., I> {};
};
template<std::size_t N>
struct genseq : genseq<N-1>::type::template push_back<N-1> {};
template<>
struct genseq<0> : seq<> {};

template<class T, std::size_t... N>
std::array<T, sizeof...(N)> vector_to_stdarray_impl(const std::vector<T>& vec, seq<N...>)
{
   HPBC_CLOCKWORK_PRECONDITION2(vec.size() >= sizeof...(N));
   return { vec[N]... };
}
template<std::size_t SIZE, class T>
std::array<T, SIZE> vector_to_stdarray(const std::vector<T>& vec)
{
   HPBC_CLOCKWORK_PRECONDITION2(vec.size() >= SIZE);
   return vector_to_stdarray_impl(vec, typename genseq<SIZE>::type{} );
}

template<class T, std::size_t... N>
std::array<T, sizeof...(N)> array_to_stdarray_impl(const T* arr, seq<N...>)
{
   return { arr[N]... };
}
template<class T, std::size_t SIZE>
std::array<T, SIZE> array_to_stdarray(const T* arr)
{
   return array_to_stdarray_impl(arr, typename genseq<SIZE>::type{} );
}




// Note: uint_to_string() and string_to_uint() provide an easy way to
// do I/O with 128 bit (or larger) integer types.

template <typename U>
std::string uint_to_string(U number)
{
   namespace hc = hurchalla;
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   if (number == 0)
      return std::string("0");
   std::string str;
   while (number > 0) {
      char digit = static_cast<char>((number % 10) + '0');
      str.push_back(digit);
      number = number / 10;
   }
   return std::string(str.rbegin(), str.rend());
}

template <typename U>
std::string uint_to_octal_string(U number)
{
   namespace hc = hurchalla;
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   if (number == 0)
      return std::string("0");
   std::string str;
   while (number > 0) {
      char digit = static_cast<char>((number % 8) + '0');
      str.push_back(digit);
      number = number / 8;
   }
   str.push_back('0');  // octal numbers are prefixed with '0' (we will be reversing str on next line)
   return std::string(str.rbegin(), str.rend());
}


struct STUException : public std::runtime_error
{
   STUException(std::string const& msg) : std::runtime_error(msg) {}
};

template <typename U>
U string_to_uint(const std::string& str)
{
   namespace hc = hurchalla;
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   constexpr U maxU = hc::ut_numeric_limits<U>::max();
   U number = 0;
   for (const auto& c : str) {
      char ch = static_cast<char>(c);
      if (ch < '0' || ch > '9')
         throw STUException("string_to_uint() called with invalid argument:"
               " non-digit character found in 'str'");
      U digit = static_cast<U>(ch - '0');
      if (number > (maxU - digit) / 10)
         throw STUException("string_to_uint() called with invalid argument:"
               " the contents of 'str' would convert to a value that is too"
               " large to fit in type 'U'");
      number = 10 * number + digit;
   }
   return number;
}




template <typename U>
U generate_random_value(std::mt19937_64& gen,
                        std::uniform_int_distribution<uint64_t>& distrib64)
{
   static_assert(hurchalla::ut_numeric_limits<U>::is_integer, "");
   static_assert(!hurchalla::ut_numeric_limits<U>::is_signed, "");
   static_assert(hurchalla::ut_numeric_limits<U>::digits <= 128, "");
   if HURCHALLA_CPP17_CONSTEXPR (hurchalla::ut_numeric_limits<U>::digits > 64) {
      uint64_t u1 = distrib64(gen);
      uint64_t u2 = distrib64(gen);
      using P = typename hurchalla::safely_promote_unsigned<U>::type;
      U val = static_cast<U>((static_cast<P>(u2) << 64u) | u1);
      return val;
   } else {
      return static_cast<U>(distrib64(gen));
   }
}





template <size_t TABLE_BITS, bool USE_SLIDING_WINDOW_OPTIMIZATION,
          size_t CODE_SECTION, class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION,
          typename ST>
int test_correctness_pow(ST seed)
{
   namespace hc = hurchalla;
   using U = typename MontType::IntegerType;
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");
   using V = typename MontType::MontgomeryValue;

   constexpr U maxU = hc::ut_numeric_limits<U>::max();
   U range = static_cast<U>(100);

   constexpr U maxMF = MontType::max_modulus();
   auto mod_range = static_cast<U>(range);
   if (mod_range >= maxMF)
      mod_range = maxMF - 1;

   std::mt19937_64 gen(seed);
   std::uniform_int_distribution<uint64_t> distrib64;

   {
      U mod = 123;
      MontType mf(mod);

      U exponent = 0;
      U base = 0;
      V mont_base = mf.convertIn(base);
      auto mont_result = hc::experimental::experimental_montgomery_pow_2kary::call<
                   MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                   USE_SQUARING_VALUE_OPTIMIZATION>(mf, mont_base, exponent);
      if (mf.getCanonicalValue(mont_result) != mf.getUnityValue()) {
         std::cout << "bug1 in montgomery_pow_2kary found: got wrong result for ";
         std::cout << uint_to_string(base) << "^" << uint_to_string(exponent) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }

      exponent = generate_random_value<U>(gen, distrib64);
      if (exponent == 0)
         ++exponent;
      if (exponent < 128)
         exponent += 128;
      mont_result = hc::experimental::experimental_montgomery_pow_2kary::call<
                   MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                   USE_SQUARING_VALUE_OPTIMIZATION>(mf, mont_base, exponent);
      if (mf.getCanonicalValue(mont_result) != mf.getZeroValue()) {
         std::cout << "bug2 in montgomery_pow_2kary found: got wrong result for ";
         std::cout << uint_to_string(base) << "^" << uint_to_string(exponent) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }

      base = 1;
      mont_base = mf.convertIn(base);
      mont_result = hc::experimental::experimental_montgomery_pow_2kary::call<
                   MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                   USE_SQUARING_VALUE_OPTIMIZATION>(mf, mont_base, exponent);
      if (mf.getCanonicalValue(mont_result) != mf.getUnityValue()) {
         std::cout << "bug3 in montgomery_pow_2kary found: got wrong result for ";
         std::cout << uint_to_string(base) << "^" << uint_to_string(exponent) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }

      base = mod - 1;
      mont_base = mf.convertIn(base);
      mont_result = hc::experimental::experimental_montgomery_pow_2kary::call<
                   MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                   USE_SQUARING_VALUE_OPTIMIZATION>(mf, mont_base, exponent);
      bool success = (exponent % 2 == 0) ? mf.getCanonicalValue(mont_result) == mf.getUnityValue()
                                         : mf.getCanonicalValue(mont_result) == mf.getNegativeOneValue();
      if (!success) {
         std::cout << "bug4 in montgomery_pow_2kary found: got wrong result for ";
         std::cout << uint_to_string(base) << "^" << uint_to_string(exponent) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }
   }


   for (U i = 0; i < mod_range - 2; ++i) {
      U mod = (i % 2 == 0) ? maxMF - i : i + 2;
      MontType mf(mod);
      U exponent = (i % 3 == 0) ? i : maxU - i;
      U base = generate_random_value<U>(gen, distrib64);
      while (base >= mod)
         base /= 2;
      V mont_base = mf.convertIn(base);
      auto mont_result = hc::experimental::experimental_montgomery_pow_2kary::call<
                MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION>(mf, mont_base, exponent);
      U result = mf.convertOut(mont_result);
      U standard_result = mf.convertOut(mf.pow(mont_base, exponent));
      if (result != standard_result) {
         std::cout << "bug5 in montgomery_pow_2kary found: got wrong result for ";
         std::cout << uint_to_string(base) << "^" << uint_to_string(exponent) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }
   }
   for (U i = 0; i < mod_range/2; ++i) {
      U mod = generate_random_value<U>(gen, distrib64);
      while (mod > maxMF)
         mod /= 2;
      if (mod % 2 == 0)
         ++mod;
      if (mod < 3)
         mod = 3;
      MontType mf(mod);
      U exponent = generate_random_value<U>(gen, distrib64);
      U base = generate_random_value<U>(gen, distrib64);
      V mont_base = mf.convertIn(base);
      auto mont_result = hc::experimental::experimental_montgomery_pow_2kary::call<
                MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION>(mf, mont_base, exponent);
      U result = mf.convertOut(mont_result);
      U standard_result = mf.convertOut(mf.pow(mont_base, exponent));
      if (result != standard_result) {
         std::cout << "bug6 in montgomery_pow_2kary found: got wrong result for ";
         std::cout << uint_to_string(base) << "^" << uint_to_string(exponent) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }
   }
   return 0;
}


template <size_t TABLE_BITS, size_t CODE_SECTION, size_t ARRAY_SIZE, class MontType,
          bool USE_SQUARING_VALUE_OPTIMIZATION, typename ST>
int test_correctness_array_pow(ST seed)
{
   namespace hc = hurchalla;
   using U = typename MontType::IntegerType;
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");
   using V = typename MontType::MontgomeryValue;

   constexpr U range = static_cast<U>(100);

   constexpr U maxMF = MontType::max_modulus();
   auto mod_range = static_cast<U>(range);
   if (mod_range >= maxMF)
      mod_range = maxMF - 1;

   std::mt19937_64 gen(seed);
   std::uniform_int_distribution<uint64_t> distrib64;

   mod_range -= 16;
   for (U i = 0; i < mod_range - 2; ++i) {
      U mod = (i % 2 == 0) ? maxMF - i : i + 2;
      // We use std::vector to indirectly make a MontType array, since
      // MontType has no default constructor.
      std::vector<MontType> mf_vec;
      std::array<U, ARRAY_SIZE> exponent_arr;
      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         if (mod >= 3 + 2*j)
            mf_vec.emplace_back(mod - 2*j);
         else
            mf_vec.emplace_back(3);
         exponent_arr[j] = mod + j * 100000;   // overflow is ok here
      }
      std::array<MontType, ARRAY_SIZE> mf_arr = vector_to_stdarray<ARRAY_SIZE>(mf_vec);

      std::array<U, ARRAY_SIZE> base_arr;
      std::array<V, ARRAY_SIZE> mont_base_arr;
      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         base_arr[j] = generate_random_value<U>(gen, distrib64);
         while (base_arr[j] >= mf_arr[j].getModulus())
            base_arr[j] /= 2;
         mont_base_arr[j] = mf_arr[j].convertIn(base_arr[j]);
      }

      auto mont_result_arr = hc::experimental::experimental_montgomery_pow_2kary::call
                  <MontType, U, ARRAY_SIZE, TABLE_BITS, CODE_SECTION, USE_SQUARING_VALUE_OPTIMIZATION>
                  (mf_arr, mont_base_arr, exponent_arr);

      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         U result = mf_arr[j].convertOut(mont_result_arr[j]);
         U standard_result = mf_arr[j].convertOut(mf_arr[j].pow(mont_base_arr[j], exponent_arr[j]));
         if (result != standard_result) {
            std::cout << "bug7 in array montgomery_pow_2kary found: got wrong result for ";
            std::cout << uint_to_string(base_arr[j]) << "^" << uint_to_string(exponent_arr[j]) << " (mod " <<
                  uint_to_string(mf_arr[j].getModulus()) << ")\n";
            return 1;
         }
      }
   }
   for (U i = 0; i < mod_range/2; ++i) {
      // We use std::vector to indirectly make a MontType array, since
      // MontType has no default constructor.
      std::vector<MontType> mf_vec;
      std::array<U, ARRAY_SIZE> exponent_arr;
      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         U mod = generate_random_value<U>(gen, distrib64);
         while (mod > maxMF)
            mod /= 2;
         if (mod % 2 == 0)
            ++mod;
         if (mod < 3)
            mod = 3;
         mf_vec.emplace_back(mod);
         exponent_arr[j] = generate_random_value<U>(gen, distrib64);
      }
      std::array<MontType, ARRAY_SIZE> mf_arr = vector_to_stdarray<ARRAY_SIZE>(mf_vec);

      std::array<U, ARRAY_SIZE> base_arr;
      std::array<V, ARRAY_SIZE> mont_base_arr;
      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         base_arr[j] = generate_random_value<U>(gen, distrib64);
         while (base_arr[j] >= mf_arr[j].getModulus())
            base_arr[j] /= 2;
         mont_base_arr[j] = mf_arr[j].convertIn(base_arr[j]);
      }

      auto mont_result_arr = hc::experimental::experimental_montgomery_pow_2kary::call
                  <MontType, U, ARRAY_SIZE, TABLE_BITS, CODE_SECTION, USE_SQUARING_VALUE_OPTIMIZATION>
                  (mf_arr, mont_base_arr, exponent_arr);

      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         U result = mf_arr[j].convertOut(mont_result_arr[j]);
         U standard_result = mf_arr[j].convertOut(mf_arr[j].pow(mont_base_arr[j], exponent_arr[j]));
         if (result != standard_result) {
            std::cout << "bug8 in array montgomery_pow_2kary found: got wrong result for ";
            std::cout << uint_to_string(base_arr[j]) << "^" << uint_to_string(exponent_arr[j]) << " (mod " <<
                  uint_to_string(mf_arr[j].getModulus()) << ")\n";
            return 1;
         }
      }
   }
   return 0;
}


template <class PTAG, size_t TABLE_BITS, size_t CODE_SECTION, size_t ARRAY_SIZE, class MontType,
          bool USE_SQUARING_VALUE_OPTIMIZATION, bool USE_SLIDING_WINDOW_OPTIMIZATION,
          typename ST>
int test_correctness_partial_array_pow(ST seed)
{
   namespace hc = hurchalla;
   using U = typename MontType::IntegerType;
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");
   using V = typename MontType::MontgomeryValue;

   constexpr U maxU = hc::ut_numeric_limits<U>::max();
   constexpr U range = static_cast<U>(100);

   constexpr U maxMF = MontType::max_modulus();
   auto mod_range = static_cast<U>(range);
   if (mod_range >= maxMF)
      mod_range = maxMF - 1;

   std::mt19937_64 gen(seed);
   std::uniform_int_distribution<uint64_t> distrib64;

   for (U i = 0; i < mod_range - 2; ++i) {
      U mod = (i % 2 == 0) ? maxMF - i : i + 2;
      MontType mf(mod);
      U exponent = (i % 3 == 0) ? i : maxU - i;

      std::array<U, ARRAY_SIZE> base_arr;
      std::array<V, ARRAY_SIZE> mont_base_arr;
      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         base_arr[j] = generate_random_value<U>(gen, distrib64);
         while (base_arr[j] >= mf.getModulus())
            base_arr[j] /= 2;
         mont_base_arr[j] = mf.convertIn(base_arr[j]);
      }

      auto mont_result_arr = hc::experimental::experimental_montgomery_pow_2kary::call<
                MontType, U, ARRAY_SIZE, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION, PTAG>(mf, mont_base_arr, exponent);
      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         U result = mf.convertOut(mont_result_arr[j]);
         U standard_result = mf.convertOut(mf.pow(mont_base_arr[j], exponent));
         if (result != standard_result) {
            std::cout << "bug9 in partial array montgomery_pow_2kary found: got wrong result for ";
            std::cout << uint_to_string(base_arr[j]) << "^" << uint_to_string(exponent) << " (mod " <<
                  uint_to_string(mf.getModulus()) << ")\n";
            return 1;
         }
      }
   }

   for (U i = 0; i < mod_range/2; ++i) {
      U mod = generate_random_value<U>(gen, distrib64);
      while (mod > maxMF)
         mod /= 2;
      if (mod % 2 == 0)
         ++mod;
      if (mod < 3)
         mod = 3;
      MontType mf(mod);
      U exponent = generate_random_value<U>(gen, distrib64);

      std::array<U, ARRAY_SIZE> base_arr;
      std::array<V, ARRAY_SIZE> mont_base_arr;
      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         base_arr[j] = generate_random_value<U>(gen, distrib64);
         mont_base_arr[j] = mf.convertIn(base_arr[j]);
      }

      auto mont_result_arr = hc::experimental::experimental_montgomery_pow_2kary::call<
                MontType, U, ARRAY_SIZE, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION, PTAG>(mf, mont_base_arr, exponent);

      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         U result = mf.convertOut(mont_result_arr[j]);
         U standard_result = mf.convertOut(mf.pow(mont_base_arr[j], exponent));
         if (result != standard_result) {
            std::cout << "bug10 in partial array montgomery_pow_2kary found: got wrong result for ";
            std::cout << uint_to_string(base_arr[j]) << "^" << uint_to_string(exponent) << " (mod " <<
                  uint_to_string(mf.getModulus()) << ")\n";
            return 1;
         }
      }
   }
   return 0;
}




struct TimingPA {
   bool is_low_uops;
   size_t table_bits;
   bool uses_sliding_window;
   size_t code_section;
   size_t array_size;
   std::chrono::duration<double>::rep time;
   bool uses_squaring_values;
   TimingPA(bool is_low_uops1, size_t table_bits1, bool uses_sliding_window1, size_t code_section1, size_t array_size1,
           std::chrono::duration<double>::rep time1, bool uses_squaring_values1)
      : is_low_uops(is_low_uops1), table_bits(table_bits1), uses_sliding_window(uses_sliding_window1),
        code_section(code_section1), array_size(array_size1),
        time(time1), uses_squaring_values(uses_squaring_values1) {}
   TimingPA() : is_low_uops(true), table_bits(0), uses_sliding_window(false), code_section(0),
               array_size(0), time(0.0), uses_squaring_values(false) {}
};


template <class PTAG, size_t TABLE_BITS, size_t CODE_SECTION, size_t ARRAY_SIZE,
          class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION, bool USE_SLIDING_WINDOW_OPTIMIZATION,
          typename U, typename ST>
TimingPA
bench_partial_array_pow(U min, U range, U& totalU, unsigned int max_modulus_bits_reduce, ST seed, unsigned int exponent_bits_reduce)
{
   HPBC_CLOCKWORK_PRECONDITION2(max_modulus_bits_reduce <
                     hurchalla::ut_numeric_limits<decltype(MontType::max_modulus())>::digits);
#if 1
   // run very short tests to hopefully catch a bugged experimental impl
   int tcpap_result = test_correctness_partial_array_pow<PTAG,
                                                     TABLE_BITS,
                                                     CODE_SECTION,
                                                     ARRAY_SIZE,
                                                     MontType,
                                                     USE_SQUARING_VALUE_OPTIMIZATION,
                                                     USE_SLIDING_WINDOW_OPTIMIZATION>(seed);
   if (tcpap_result != 0) {
      std::cout << "Failed on TABLE_BITS == " << TABLE_BITS;
      std::cout << ", CODE_SECTION == " << CODE_SECTION;
      if (USE_SQUARING_VALUE_OPTIMIZATION)
         std::cout << ", USE_SQUARING_VALUE_OPTIMIZATION == true";
      else
         std::cout << ", USE_SQUARING_VALUE_OPTIMIZATION == false";
      if (USE_SLIDING_WINDOW_OPTIMIZATION)
         std::cout << ", USE_SLIDING_WINDOW_OPTIMIZATION == true";
      else
         std::cout << ", USE_SLIDING_WINDOW_OPTIMIZATION == false";
      std::cout << ", ARRAY_SIZE == " << ARRAY_SIZE;

      using MontTag = typename MontType::MontType::MontyTag;
      if (std::is_same<MontTag, ::hurchalla::detail::TagMontyFullrangeMasked>::value)
         std::cout << ", MontTag == TagMontyFullrangeMasked";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyFullrange>::value)
         std::cout << ", MontTag == TagMontyFullrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyHalfrange>::value)
         std::cout << ", MontTag == TagMontyHalfrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyQuarterrange>::value)
         std::cout << ", MontTag == TagMontyQuarterrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyWrappedmath>::value)
         std::cout << ", MontTag == TagMontyWrappedmath";
      else
         std::cout << ", MontTag == UNKNOWN";

      if (std::is_same<PTAG, ::hurchalla::LowlatencyTag>::value)
         std::cout << ", PTAG == LowlatencyTag";
      else if (std::is_same<PTAG, ::hurchalla::LowuopsTag>::value)
         std::cout << ", PTAG == LowuopsTag";
      else
         std::cout << ", PTAG == UNKNOWN";
      std::cout << "\n";

      exit(1);
   }
# ifdef TEST_CORRECTNESS_ONLY
   return TimingPA();
# endif
#endif

   using V = typename MontType::MontgomeryValue;

   using namespace std::chrono;
   using dsec = duration<double>;

   constexpr auto max_modulus = MontType::max_modulus();

   constexpr bool randomizeModuli = true;
   unsigned int exponentreduction = exponent_bits_reduce;
   range *= 2;
   range = range / ARRAY_SIZE;

#if 1
   U maxMod = max_modulus >> max_modulus_bits_reduce;
#else
   if (range >= max_modulus/2)
      range = (max_modulus/2) - 1;
   U maxMod = 1 + range + (max_modulus >> max_modulus_bits_reduce);
#endif
   maxMod = maxMod - ((maxMod + 1) % 2);

   U max;
   if (range > maxMod) {
      min = 0;
      max = maxMod;
   } else {
      // if (min + range > maxMod) // this might possibly overflow
      if (min > maxMod - range)
         min = maxMod - range;
      max = min + range;
   }
   if (max % 2 == 0)
      --max;
   if (min % 2 == 0)
      ++min;
   while ((max - min) % 8 != 0)
      min += 2;

   {
      HPBC_CLOCKWORK_ASSERT(max > 0);
      int leading_zeros = hurchalla::count_leading_zeros(max);
      int numbits = hurchalla::ut_numeric_limits<U>::digits - leading_zeros;
      HPBC_CLOCKWORK_ASSERT(numbits > 0);
      U maxmask = static_cast<U>(1) << (numbits - 1);  //we need numbits-1 since numbits may be big enough to be UB to shift by.
      maxmask = 2 * maxmask - 1;

      std::mt19937_64 gen(seed);
      std::uniform_int_distribution<uint64_t> distrib64;

      std::vector<U> tmpvec;
      if (!randomizeModuli) {
         for (U x = max; x > min; x = x-2)
            tmpvec.push_back(x);
      } else {
         for (U x = max; x > min; x = x-2) {
            U val;
            // this while loop is somewhat of a hack, but should be good enough...
            do {
               val = generate_random_value<U>(gen, distrib64);
               val = val & maxmask;
            } while (val > max || val < max/2 || (val % 2 == 0));
            tmpvec.push_back(val);
         }
      }

      std::vector<MontType> mfvec;
      std::vector<V> randbaseV;
      for (size_t i=0; i<tmpvec.size(); ++i) {
         U x = tmpvec[i];
         MontType mf(x);
         mfvec.push_back(mf);

         for (size_t j=0; j < ARRAY_SIZE; ++j) {
            U base = generate_random_value<U>(gen, distrib64);
//            while (base >= mf.getModulus())
//               base /= 2;
            randbaseV.push_back(mf.convertIn(base));
         }
      }

      std::vector<U> randexpU;
      for (U x = max; x > min; x = x-2) {
         U val = generate_random_value<U>(gen, distrib64);
         uint64_t ranval2 = generate_random_value<uint64_t>(gen, distrib64);
         unsigned int extra_reduce = ranval2 & 7;
         U exponentmask = static_cast<U>(static_cast<U>(0) - static_cast<U>(1)) >> (exponentreduction + extra_reduce);
         val = val & exponentmask;
         if (val < exponentmask/2)
            val += exponentmask/2;
         randexpU.push_back(val);
      }

      auto t0_mfp = steady_clock::now();

      for (size_t i=0; i<mfvec.size(); ++i) {

#if 1   // Skip adding any cost of making the monts
         MontType mf = mfvec[i];
#elif 1 // Add the cost of making monts without RSquaredModN
         MontType mf = mfvec[i];
         {
            U x = tmpvec[i];
            U r_mod_n = ::hurchalla::get_R_mod_n(x);
            U inv_n = ::hurchalla::inverse_mod_R(x);
            totalU += r_mod_n + inv_n;
#  if 0 // Add the cost of making RSquaredModN
            U rsmn = ::hurchalla::get_Rsquared_mod_n(x, inv_n, r_mod_n);
            totalU += rsmn;
#  endif
         }
#else   // Add the full cost of making the monts
         MontType mf(tmpvec[i]);
#endif

         U exponent = randexpU[i];

         std::array<V, ARRAY_SIZE> mont_base_arr;
         for (size_t j=0; j < ARRAY_SIZE; ++j)
            mont_base_arr[j] = randbaseV[(i*ARRAY_SIZE) + j];

#if 0
// for comparison...
         std::array<V, ARRAY_SIZE> result = mf.pow(mont_base_arr, exponent);
#else
         std::array<V, ARRAY_SIZE> result = hurchalla::experimental::experimental_montgomery_pow_2kary::call<
                MontType, U, ARRAY_SIZE, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION, PTAG>(mf, mont_base_arr, exponent);
#endif

#if 0
         for (size_t j=0; j < ARRAY_SIZE; ++j)
            totalU += mf.convertOut(result[j]);
#else
         if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrangeMasked>::value) {
            struct OpenV : public V {
               HURCHALLA_FORCE_INLINE OpenV(V x) : V(x) {}
               HURCHALLA_FORCE_INLINE auto getbits() -> decltype(V::getbits()) { return V::getbits(); }
            };
            for (size_t j=0; j < ARRAY_SIZE; ++j)
               totalU += OpenV(result[j]).getbits();
         }
         else {
            struct OpenV : public V {
               HURCHALLA_FORCE_INLINE OpenV(V x) : V(x) {}
               HURCHALLA_FORCE_INLINE auto get() -> decltype(V::get()) { return V::get(); }
            };
            for (size_t j=0; j < ARRAY_SIZE; ++j)
               totalU += static_cast<U>(OpenV(result[j]).get());
         }
#endif
      }

      auto t1_mfp = steady_clock::now();
      dsec::rep mtp_time = dsec(t1_mfp-t0_mfp).count();

      return TimingPA(std::is_same<PTAG, ::hurchalla::LowuopsTag>::value,
                     TABLE_BITS, USE_SLIDING_WINDOW_OPTIMIZATION,
                     CODE_SECTION, ARRAY_SIZE, mtp_time, USE_SQUARING_VALUE_OPTIMIZATION);
   }
}




struct TimingA {
   size_t table_bits;
   size_t code_section;
   size_t array_size;
   std::chrono::duration<double>::rep time;
   bool uses_squaring_values;
   TimingA(size_t table_bits1, size_t code_section1, size_t array_size1,
           std::chrono::duration<double>::rep time1, bool uses_squaring_values1)
      : table_bits(table_bits1), code_section(code_section1), array_size(array_size1),
        time(time1), uses_squaring_values(uses_squaring_values1) {}
   TimingA() : table_bits(0), code_section(0), array_size(0), time(0.0), uses_squaring_values(false) {}
};


template <size_t TABLE_BITS, size_t CODE_SECTION, size_t ARRAY_SIZE,
          class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION,
          typename U, typename ST>
TimingA
bench_array_pow(U min, U range, U& totalU, unsigned int max_modulus_bits_reduce, ST seed, unsigned int exponent_bits_reduce)
{
   HPBC_CLOCKWORK_PRECONDITION2(max_modulus_bits_reduce <
                     hurchalla::ut_numeric_limits<decltype(MontType::max_modulus())>::digits);
#if 1
   // run very short tests to hopefully catch a bugged experimental impl
   int tcap_result = test_correctness_array_pow<TABLE_BITS,
                                                     CODE_SECTION,
                                                     ARRAY_SIZE,
                                                     MontType,
                                                     USE_SQUARING_VALUE_OPTIMIZATION>(seed);
   if (tcap_result != 0) {
      std::cout << "Failed on TABLE_BITS == " << TABLE_BITS;
      std::cout << ", CODE_SECTION == " << CODE_SECTION;
      if (USE_SQUARING_VALUE_OPTIMIZATION)
         std::cout << ", USE_SQUARING_VALUE_OPTIMIZATION == true";
      else
         std::cout << ", USE_SQUARING_VALUE_OPTIMIZATION == false";
      std::cout << ", ARRAY_SIZE == " << ARRAY_SIZE;
      using MontTag = typename MontType::MontType::MontyTag;
      if (std::is_same<MontTag, ::hurchalla::detail::TagMontyFullrangeMasked>::value)
         std::cout << ", MontTag == TagMontyFullrangeMasked";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyFullrange>::value)
         std::cout << ", MontTag == TagMontyFullrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyHalfrange>::value)
         std::cout << ", MontTag == TagMontyHalfrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyQuarterrange>::value)
         std::cout << ", MontTag == TagMontyQuarterrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyWrappedmath>::value)
         std::cout << ", MontTag == TagMontyWrappedmath";
      else
         std::cout << ", MontTag == UNKNOWN";
      std::cout << "\n";

      exit(1);
   }
# ifdef TEST_CORRECTNESS_ONLY
   return TimingA();
# endif
#endif

   using V = typename MontType::MontgomeryValue;

   using namespace std::chrono;
   using dsec = duration<double>;

   constexpr auto max_modulus = MontType::max_modulus();

   constexpr bool randomizeModuli = true;
   unsigned int exponentreduction = exponent_bits_reduce;
   range *= 2;

#if 1
   U maxMod = max_modulus >> max_modulus_bits_reduce;
#else
   if (range >= max_modulus/2)
      range = (max_modulus/2) - 1;
   U maxMod = 1 + range + (max_modulus >> max_modulus_bits_reduce);
#endif
   maxMod = maxMod - ((maxMod + 1) % 2);

   U max;
   if (range > maxMod) {
      min = 0;
      max = maxMod;
   } else {
      // if (min + range > maxMod) // this might possibly overflow
      if (min > maxMod - range)
         min = maxMod - range;
      max = min + range;
   }
   if (max % 2 == 0)
      --max;
   if (min % 2 == 0)
      ++min;
   while ((max - min) % 8 != 0)
      min += 2;


   {
      HPBC_CLOCKWORK_ASSERT(max > 0);
      int leading_zeros = hurchalla::count_leading_zeros(max);
      int numbits = hurchalla::ut_numeric_limits<U>::digits - leading_zeros;
      HPBC_CLOCKWORK_ASSERT(numbits > 0);
      U maxmask = static_cast<U>(1) << (numbits - 1);  //we need numbits-1 since numbits may be big enough to be UB to shift by.
      maxmask = 2 * maxmask - 1;

      std::mt19937_64 gen(seed);
      std::uniform_int_distribution<uint64_t> distrib64;

      std::vector<U> tmpvec;
      if (!randomizeModuli) {
         for (U x = max; x > min; x = x-2)
            tmpvec.push_back(x);
      } else {
         for (U x = max; x > min; x = x-2) {
            U val;
            // this while loop is somewhat of a hack, but should be good enough...
            do {
               val = generate_random_value<U>(gen, distrib64);
               val = val & maxmask;
            } while (val > max || val < max/2 || (val % 2 == 0));
            tmpvec.push_back(val);
         }
      }

      std::vector<MontType> mfvec;
      std::vector<V> randbaseV;
      for (size_t i=0; i<tmpvec.size(); ++i) {
         U x = tmpvec[i];
         MontType mf(x);
         mfvec.push_back(mf);

         U base = generate_random_value<U>(gen, distrib64);
//         while (base >= mf.getModulus())
//            base /= 2;
         randbaseV.push_back(mf.convertIn(base));
      }

      std::vector<U> randexpU;
      for (U x = max; x > min; x = x-2) {
         U val = generate_random_value<U>(gen, distrib64);
         uint64_t ranval2 = generate_random_value<uint64_t>(gen, distrib64);
         unsigned int extra_reduce = ranval2 & 7;
         U exponentmask = static_cast<U>(static_cast<U>(0) - static_cast<U>(1)) >> (exponentreduction + extra_reduce);
         val = val & exponentmask;
         if (val < exponentmask/2)
            val += exponentmask/2;
         randexpU.push_back(val);
      }

      std::vector<MontType> mfvec2;
      mfvec2.reserve(32);

      auto t0_mfp = steady_clock::now();

      for (size_t i=0; i + ARRAY_SIZE - 1 < mfvec.size(); i += ARRAY_SIZE) {
#if 1   // Skip adding any cost of making the monts
         std::array<MontType, ARRAY_SIZE> mfarr = array_to_stdarray<MontType, ARRAY_SIZE>(mfvec.data() + i);
#elif 1 // Add the cost of making monts without RSquaredModN
         std::array<MontType, ARRAY_SIZE> mfarr = array_to_stdarray<MontType, ARRAY_SIZE>(mfvec.data() + i);
         std::array<U, ARRAY_SIZE> r_mod_n;
         std::array<U, ARRAY_SIZE> inv_n;
         for (size_t j=0; j < ARRAY_SIZE; ++j) {
            U x = tmpvec[i+j];
            r_mod_n[j] = ::hurchalla::get_R_mod_n(x);
            inv_n[j] = ::hurchalla::inverse_mod_R(x);
            totalU += r_mod_n[j] + inv_n[j];
         }
#  if 1 // Add the cost of making RSquaredModN
#    if 1
         std::array<U, ARRAY_SIZE> n_arr;
         for (size_t j=0; j < ARRAY_SIZE; ++j) {
            n_arr[j] = tmpvec[i+j];
         }
         constexpr bool isQR = std::is_same<hurchalla::MontgomeryQuarter<U>, MontType>::value;
         std::array<U, ARRAY_SIZE> arr_rsmn = hurchalla::detail::impl_array_get_Rsquared_mod_n
                            <isQR, hurchalla::LowuopsTag>::call(n_arr, inv_n, r_mod_n);
         for (size_t j=0; j < ARRAY_SIZE; ++j) {
            totalU += arr_rsmn[j];
         }
#    elif 0
         for (size_t j=0; j < ARRAY_SIZE; ++j) {
            U modulus = tmpvec[i+j];
            U u_lo;
            U u_hi = hurchalla::unsigned_square_to_hilo_product(u_lo, r_mod_n[j]);
            U remainder;
            //U quotient = div_2U_by_1U(u_hi, u_lo, modulus, remainder);
            div_2U_by_1U(u_hi, u_lo, modulus, remainder);
            totalU += remainder;
         }
#    else
         for (size_t j=0; j < ARRAY_SIZE; ++j) {
            U x = tmpvec[i+j];
            U rsmn = ::hurchalla::get_Rsquared_mod_n(x,inv_n[j],r_mod_n[j]);
            totalU += rsmn;
         }
#    endif
#  endif
#else   // Add the full cost of making the monts
         mfvec2.clear();
         for (size_t j=0; j < ARRAY_SIZE; ++j) {
            U x = tmpvec[i+j];
            mfvec2.emplace_back(x);
         }
         std::array<MontType, ARRAY_SIZE> mfarr = vector_to_stdarray<ARRAY_SIZE>(mfvec2);
#endif

         std::array<U, ARRAY_SIZE> exparr;
         std::array<V, ARRAY_SIZE> mont_base_arr;
         for (size_t j=0; j < ARRAY_SIZE; ++j) {
            exparr[j] = randexpU[i + j];
            mont_base_arr[j] = randbaseV[i + j];
         }

         std::array<V, ARRAY_SIZE> result = hurchalla::experimental::experimental_montgomery_pow_2kary::call
               <MontType, U, ARRAY_SIZE, TABLE_BITS, CODE_SECTION, USE_SQUARING_VALUE_OPTIMIZATION>(mfarr, mont_base_arr, exparr);

#if 0
         for (size_t j=0; j < ARRAY_SIZE; ++j)
            totalU += mfarr[j].convertOut(result[j]);
#else
         if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrangeMasked>::value) {
            struct OpenV : public V {
               HURCHALLA_FORCE_INLINE OpenV(V x) : V(x) {}
               HURCHALLA_FORCE_INLINE auto getbits() -> decltype(V::getbits()) { return V::getbits(); }
            };
            for (size_t j=0; j < ARRAY_SIZE; ++j)
               totalU += OpenV(result[j]).getbits();
         }
         else {
            struct OpenV : public V {
               HURCHALLA_FORCE_INLINE OpenV(V x) : V(x) {}
               HURCHALLA_FORCE_INLINE auto get() -> decltype(V::get()) { return V::get(); }
            };
            for (size_t j=0; j < ARRAY_SIZE; ++j)
               totalU += static_cast<U>(OpenV(result[j]).get());
         }
#endif
      }

      auto t1_mfp = steady_clock::now();
      dsec::rep mtp_time = dsec(t1_mfp-t0_mfp).count();

      return TimingA(TABLE_BITS, CODE_SECTION, ARRAY_SIZE, mtp_time, USE_SQUARING_VALUE_OPTIMIZATION);
   }

}




struct Timing {
   size_t table_bits;
   bool uses_sliding_window;
   size_t code_section;
   std::chrono::duration<double>::rep time;
   bool uses_squaring_values;
   Timing(size_t table_bits1, bool uses_sliding_window1, size_t code_section1,
          std::chrono::duration<double>::rep time1, bool uses_squaring_values1)
      : table_bits(table_bits1), uses_sliding_window(uses_sliding_window1), code_section(code_section1),
        time(time1), uses_squaring_values(uses_squaring_values1) {}
   Timing() : table_bits(0), uses_sliding_window(false), code_section(0), time(0.0), uses_squaring_values(false) {}
};


template <size_t TABLE_BITS, bool USE_SLIDING_WINDOW_OPTIMIZATION,
          size_t CODE_SECTION, class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION,
          typename U, typename ST>
Timing
bench_range(U min, U range, U& totalU, unsigned int max_modulus_bits_reduce, ST seed, unsigned int exponent_bits_reduce)
{
   HPBC_CLOCKWORK_PRECONDITION2(max_modulus_bits_reduce <
                     hurchalla::ut_numeric_limits<decltype(MontType::max_modulus())>::digits);
#if 1
   // run very short tests to hopefully catch a bugged experimental impl
   int tcp_result = test_correctness_pow<TABLE_BITS,
                                             USE_SLIDING_WINDOW_OPTIMIZATION,
                                             CODE_SECTION,
                                             MontType,
                                             USE_SQUARING_VALUE_OPTIMIZATION>(seed);
   if (tcp_result != 0) {
      std::cout << "Failed on TABLE_BITS == " << TABLE_BITS;
      if (USE_SLIDING_WINDOW_OPTIMIZATION)
         std::cout << ", USE_SLIDING_WINDOW_OPTIMIZATION == true";
      else
         std::cout << ", USE_SLIDING_WINDOW_OPTIMIZATION == false";
      std::cout << ", CODE_SECTION == " << CODE_SECTION;
      if (USE_SQUARING_VALUE_OPTIMIZATION)
         std::cout << ", USE_SQUARING_VALUE_OPTIMIZATION == true";
      else
         std::cout << ", USE_SQUARING_VALUE_OPTIMIZATION == false";
      using MontTag = typename MontType::MontType::MontyTag;
      if (std::is_same<MontTag, ::hurchalla::detail::TagMontyFullrangeMasked>::value)
         std::cout << ", MontTag == TagMontyFullrangeMasked";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyFullrange>::value)
         std::cout << ", MontTag == TagMontyFullrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyHalfrange>::value)
         std::cout << ", MontTag == TagMontyHalfrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyQuarterrange>::value)
         std::cout << ", MontTag == TagMontyQuarterrange";
      else if (std::is_same<MontTag, ::hurchalla::detail::TagMontyWrappedmath>::value)
         std::cout << ", MontTag == TagMontyWrappedmath";
      else
         std::cout << ", MontTag == UNKNOWN";
      std::cout << "\n";

      exit(1);
   }
# ifdef TEST_CORRECTNESS_ONLY
   return Timing();
# endif
#endif

   using namespace std::chrono;
   using dsec = duration<double>;

   constexpr auto max_modulus = MontType::max_modulus();

   constexpr bool randomizeModuli = true;
   unsigned int exponentreduction = exponent_bits_reduce;
   range *= 2;

#if 1
   U maxMod = max_modulus >> max_modulus_bits_reduce;
#else
   if (range >= max_modulus/2)
      range = (max_modulus/2) - 1;
   U maxMod = 1 + range + max_modulus >> max_modulus_bits_reduce;
#endif
   maxMod = maxMod - ((maxMod + 1) % 2);


   U max;
   if (range > maxMod) {
      min = 0;
      max = maxMod;
   } else {
      // if (min + range > maxMod) // this might possibly overflow
      if (min > maxMod - range)
         min = maxMod - range;
      max = min + range;
   }
   if (max % 2 == 0)
      --max;
   if (min % 2 == 0)
      ++min;
   while ((max - min) % 8 != 0)
      min += 2;

//   std::cout << "max is " << uint_to_octal_string(max) << "\n";
//   std::cout << "min is " << uint_to_octal_string(min) << "\n";

   {
      HPBC_CLOCKWORK_ASSERT(max > 0);
      int leading_zeros = hurchalla::count_leading_zeros(max);
      int numbits = hurchalla::ut_numeric_limits<U>::digits - leading_zeros;
      HPBC_CLOCKWORK_ASSERT(numbits > 0);
      U maxmask = static_cast<U>(1) << (numbits - 1);  //we need numbits-1 since numbits may be big enough to be UB to shift by.
      maxmask += maxmask - 1;


      std::mt19937_64 gen(seed);
      std::uniform_int_distribution<uint64_t> distrib64;

      std::vector<U> tmpvec;
      if (!randomizeModuli) {
         for (U x = max; x > min; x = x-2)
            tmpvec.push_back(x);
      } else {
         for (U x = max; x > min; x = x-2) {
            U val;
            // this while loop is somewhat of a hack, but should be good enough...
            do {
               val = generate_random_value<U>(gen, distrib64);
               val = val & maxmask;
            } while (val > max || val < max/2 || (val % 2 == 0));
            tmpvec.push_back(val);
         }
      }

      using V = typename MontType::MontgomeryValue;

      std::vector<MontType> mfvec;
      std::vector<V> randbaseV;
      for (size_t i=0; i<tmpvec.size(); ++i) {
         U x = tmpvec[i];
         MontType mf(x);
         mfvec.push_back(mf);

         U base = generate_random_value<U>(gen, distrib64);
//         while (base >= mf.getModulus())
//            base /= 2;
         randbaseV.push_back(mf.convertIn(base));
      }

      std::vector<U> randexpU;
      {
         for (U x = max; x > min; x = x-2) {
            U val = generate_random_value<U>(gen, distrib64);
            uint64_t ranval2 = generate_random_value<uint64_t>(gen, distrib64);
            unsigned int extra_reduce = ranval2 & 7;
            U exponentmask = static_cast<U>(static_cast<U>(0) - static_cast<U>(1)) >> (exponentreduction + extra_reduce);
            val = val & exponentmask;
            if (val < exponentmask/2)
               val += exponentmask/2;
            randexpU.push_back(val);
         }
      }

      auto t0_mfp = steady_clock::now();

      for (size_t i=0; i<mfvec.size(); ++i) {

#if 1   // Skip adding any cost of making the monts
         MontType mf = mfvec[i];
#elif 1 // Add the cost of making monts without RSquaredModN
         MontType mf = mfvec[i];
         {
            U x = tmpvec[i];
            U r_mod_n = ::hurchalla::get_R_mod_n(x);
            U inv_n = ::hurchalla::inverse_mod_R(x);
            totalU += r_mod_n + inv_n;
#  if 0 // Add the cost of making RSquaredModN
            U rsmn = ::hurchalla::get_Rsquared_mod_n(x, inv_n, r_mod_n);
            totalU += rsmn;
#  endif
         }
#else   // Add the full cost of making the monts
         MontType mf(tmpvec[i]);
#endif

         U exponent = randexpU[i];
         V mont_base = randbaseV[i];

#if 0
// for comparison...
         auto val = mf.pow(mont_base, exponent);
#else
         auto val = hurchalla::experimental::experimental_montgomery_pow_2kary::call<
                MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION>(mf, mont_base, exponent);
#endif

#if 0
         totalU += mf.convertOut(val);
#else
         if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrangeMasked>::value) {
            using V = typename MontType::MontgomeryValue;
            struct OpenV : public V {
               HURCHALLA_FORCE_INLINE OpenV(V x) : V(x) {}
               HURCHALLA_FORCE_INLINE auto getbits() -> decltype(V::getbits()) { return V::getbits(); }
            };
            totalU += OpenV(val).getbits();
         }
         else {
            using V = typename MontType::MontgomeryValue;
            struct OpenV : public V {
               HURCHALLA_FORCE_INLINE OpenV(V x) : V(x) {}
               HURCHALLA_FORCE_INLINE auto get() -> decltype(V::get()) { return V::get(); }
            };
            totalU += static_cast<U>(OpenV(val).get());
         }
#endif
      }

      auto t1_mfp = steady_clock::now();
      dsec::rep mtp_time = dsec(t1_mfp-t0_mfp).count();

      return Timing(TABLE_BITS, USE_SLIDING_WINDOW_OPTIMIZATION, CODE_SECTION, mtp_time, USE_SQUARING_VALUE_OPTIMIZATION);
   }
}



#if 0
template <class PTAG, size_t TABLE_BITS, size_t CODE_SECTION,
          class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION, bool USE_SLIDING_WINDOW_OPTIMIZATION,
          typename U, typename ST>
void bench_PA(std::vector<TimingPA>& vecTimingPA,
   U maxU, U range, U& dummy, unsigned int mmbr, ST seed, unsigned int ebr)
{
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 2, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 3, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 4, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 5, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 6, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 7, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 8, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 10, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 12, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, TABLE_BITS, CODE_SECTION, 14, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
}
#endif



template <class PTAG, size_t ARRAY_SIZE,
          class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION, bool USE_SLIDING_WINDOW_OPTIMIZATION,
          typename U, typename ST>
void bench_PA_2(std::vector<TimingPA>& vecTimingPA,
   U maxU, U range, U& dummy, unsigned int mmbr, ST seed, unsigned int ebr)
{
//   vecTimingPA.push_back(bench_partial_array_pow
//      <PTAG, TABLE_BITS, CODE_SECTION, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 6, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 7, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 8, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 9, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 10, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 11, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 12, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 13, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 14, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 15, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 16, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));

   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 6, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 7, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 8, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 9, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 10, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 11, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 12, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));

   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 4, 6, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 4, 7, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 4, 8, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 4, 9, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 4, 10, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));


   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 2, 1, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));

   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 1, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 3, 2, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));

   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 4, 1, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 4, 2, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));

   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 5, 1, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
   vecTimingPA.push_back(bench_partial_array_pow
      <PTAG, 5, 2, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION, USE_SLIDING_WINDOW_OPTIMIZATION>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
}



template <class PTAG, size_t ARRAY_SIZE, class MontType, typename U, typename ST>
void bench_PA_PTAG(std::vector<TimingPA>& vecTimingPA,
   U maxU, U range, U& dummy, unsigned int mmbr, ST seed, unsigned int ebr)
{
   namespace hc = ::hurchalla;

//template <class PTAG, size_t TABLE_BITS, size_t CODE_SECTION, size_t ARRAY_SIZE,
//          class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION, bool USE_SLIDING_WINDOW_OPTIMIZATION,
//          typename U, typename ST>    TimingPA  bench_partial_array_pow(...)

      vecTimingPA.push_back(
         bench_partial_array_pow<PTAG, 2, 0, ARRAY_SIZE, MontType, false, false>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
      vecTimingPA.push_back(
         bench_partial_array_pow<PTAG, 2, 3, ARRAY_SIZE, MontType, false, false>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
      vecTimingPA.push_back(
         bench_partial_array_pow<PTAG, 2, 4, ARRAY_SIZE, MontType, false, false>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));
      vecTimingPA.push_back(
         bench_partial_array_pow<PTAG, 2, 5, ARRAY_SIZE, MontType, false, false>(static_cast<U>(maxU - range), range, dummy, mmbr, seed, ebr));


//template <class PTAG, size_t ARRAY_SIZE,
//          class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION, bool USE_SLIDING_WINDOW_OPTIMIZATION,
//          typename U, typename ST>bench_PA_2(...)

      bench_PA_2<PTAG, ARRAY_SIZE, MontType, false, false>(vecTimingPA, maxU, range, dummy, mmbr, seed, ebr);
      bench_PA_2<PTAG, ARRAY_SIZE, MontType, false, true>(vecTimingPA, maxU, range, dummy, mmbr, seed, ebr);


   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      bench_PA_2<PTAG, ARRAY_SIZE, MontType, true, false>(vecTimingPA, maxU, range, dummy, mmbr, seed, ebr);
      bench_PA_2<PTAG, ARRAY_SIZE, MontType, true, true>(vecTimingPA, maxU, range, dummy, mmbr, seed, ebr);
   }
}



template <size_t ARRAY_SIZE, class MontType, typename U, typename ST>
void bench_PA_all(std::vector<TimingPA>& vecTimingPA,
   U maxU, U range, U& dummy, unsigned int mmbr, ST seed, unsigned int ebr)
{
   namespace hc = ::hurchalla;

//template <class PTAG, size_t ARRAY_SIZE, class MontType, typename U, typename ST>
//void bench_PA_PTAG(...

      bench_PA_PTAG<hc::LowuopsTag, ARRAY_SIZE, MontType>(vecTimingPA, maxU, range, dummy, mmbr, seed, ebr);
      bench_PA_PTAG<hc::LowlatencyTag, ARRAY_SIZE, MontType>(vecTimingPA, maxU, range, dummy, mmbr, seed, ebr);
}





int main(int argc, char** argv)
{
   namespace hc = hurchalla;
   std::cout << "---Running Program---\n";


   unsigned int randomization_seed = 1;
   if (argc > 1)
      randomization_seed = string_to_uint<unsigned int>(std::string(argv[1]));

   unsigned int max_modulus_bits_reduce = 0;
   if (argc > 2)
      max_modulus_bits_reduce = string_to_uint<unsigned int>(std::string(argv[2]));

   unsigned int exponent_bits_reduce = 0;
   if (argc > 3)
      exponent_bits_reduce = string_to_uint<unsigned int>(std::string(argv[3]));


using namespace hurchalla;


#if !defined(TEST_ARRAY) && !defined(TEST_SCALAR) && !defined(TEST_PARTIAL_ARRAY)
#  error "You must define either TEST_ARRAY or TEST_SCALAR or TEST_PARTIAL_ARRAY"
#endif
#if defined(TEST_ARRAY) && defined(TEST_SCALAR)
#  error "You can only define one of TEST_ARRAY and TEST_SCALAR and TEST_PARTIAL_ARRAY"
#endif
#if defined(TEST_ARRAY) && defined(TEST_PARTIAL_ARRAY)
#  error "You can only define one of TEST_ARRAY and TEST_SCALAR and TEST_PARTIAL_ARRAY"
#endif
#if defined(TEST_SCALAR) && defined(TEST_PARTIAL_ARRAY)
#  error "You can only define one of TEST_ARRAY and TEST_SCALAR and TEST_PARTIAL_ARRAY"
#endif


#ifndef DEF_UINT_TYPE
#  error "DEF_UINT_TYPE was not defined"
#else
   using U = DEF_UINT_TYPE;  // uint64_t  __uint128_t  etc
#endif
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!(ut_numeric_limits<U>::is_signed), "");

#ifndef DEF_MONT_TYPE
#  error "DEF_MONT_TYPE was not defined"
#else
   using MontType = DEF_MONT_TYPE<U>; // MontgomeryQuarter<U>;
#endif
/*
# if 1
   // If you can guarantee your modulus will always be less than one quarter the
   // maximum value of type U, then use MontgomeryQuarter for speed.
   using MontType = hc::MontgomeryQuarter<U>;
# elif 0
   // If you can guarantee your modulus will always be less than one half the
   // maximum value of type U, then use MontgomeryHalf might be faster than
   // MontgomeryForm.
   using MontType = hc::MontgomeryHalf<U>;
# else
   // If you can't guarantee your modulus will always be small enough for
   // MontgomeryQuarter or MontgomeryHalf, then you must use MontgomeryForm.
   using MontType = hc::MontgomeryForm<U>;
# endif
*/

   constexpr U maxU = hc::ut_numeric_limits<U>::max();
   U range = static_cast<U>(100000);



// ------- Benchmarking --------

   std::cout << std::fixed;
   std::cout.precision(4);

   auto seed = randomization_seed;

   U dummy = 0; // dummy exists to prevent the compiler from optimizing away getting timings

   std::array<unsigned int, 4> mmbr = { 0, max_modulus_bits_reduce, 0, max_modulus_bits_reduce };
   unsigned int default_ebr = (std::is_same<MontType, MontgomeryQuarter<U>>::value) ? 2 :
                                 ((std::is_same<MontType, MontgomeryHalf<U>>::value) ? 1 : 0);
   std::array<unsigned int, 4> ebr = { default_ebr, default_ebr, exponent_bits_reduce, exponent_bits_reduce };


   constexpr int NUM_TEST_REPETITIONS = 10;


#if defined(TEST_PARTIAL_ARRAY)
   std::cout << "\nbegin benchmarks - partial array pow\n";

   // warm up call before benchmarking
//   bench_partial_array_pow<LowuopsTag, 4, 1, 4, MontType, false, false>(static_cast<U>(maxU - range), range, dummy,
//                                                            max_modulus_bits_reduce, seed, exponent_bits_reduce);

   // format is bench_PA<class PTAG, size_t TABLE_BITS, size_t CODE_SECTION,
   //       class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION, bool USE_SLIDING_WINDOW_OPTIMIZATION>(...)

   std::array<std::array<std::vector<TimingPA>, NUM_TEST_REPETITIONS>, 4> timingPA;

   for (size_t i=0; i<4; ++i) {
     for (size_t j=0; j<timingPA[i].size(); ++j) {

      bench_PA_all<2, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<3, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<4, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<5, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<6, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<7, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<8, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<10, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<12, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);
      bench_PA_all<14, MontType>(timingPA[i][j], maxU, range, dummy, mmbr[i], seed, ebr[i]);

     }
   }
#ifdef TEST_CORRECTNESS_ONLY
   std::cout << "no errors found" << "\n\n";
   return 0;
#endif

   std::array<std::vector<TimingPA>, 4> best_timingPA;
   for (size_t i=0; i<4; ++i) {
      best_timingPA[i] = timingPA[i][0];

      for (const auto& tvec : timingPA[i]) {
         for (size_t j=0; j<tvec.size(); ++j) {
            if (best_timingPA[i][j].time > tvec[j].time)
               best_timingPA[i][j].time = tvec[j].time;
         }
      }
   }

   std::vector<TimingPA> overall_bestPA = best_timingPA[0];
   for (size_t j=0; j<overall_bestPA.size(); ++j) {
      for (size_t i=1; i<4; ++i) {
         overall_bestPA[j].time += best_timingPA[i][j].time;
      }
   }
   std::sort(overall_bestPA.begin(), overall_bestPA.end(), [](TimingPA const& t1, TimingPA const& t2) { return t1.time < t2.time; });

   for (size_t i=0; i<4; ++i) {
      std::sort(best_timingPA[i].begin(), best_timingPA[i].end(), [](TimingPA const& t1, TimingPA const& t2) { return t1.time < t2.time; });
   }


   std::cout << "(ignore)" << uint_to_string(dummy) << "\n\n";

std::cout << "OVERALL BEST:\n";
   for (size_t j=0; j < overall_bestPA.size(); ++j) {
         const auto& t = overall_bestPA[j];
         std::cout << 10.0 * t.time << " " << t.table_bits << " ";
         if (t.code_section < 10)
            std::cout << "0";
         std::cout << t.code_section;
         if (t.uses_sliding_window)
            std::cout <<  " t";
         else
            std::cout <<  " x";
         if (t.uses_squaring_values)
            std::cout <<  " t";
         else
            std::cout <<  " x";
         if (t.array_size < 10)
            std::cout <<  " 0" << t.array_size;
         else
            std::cout <<  " " << t.array_size;
         if (t.is_low_uops)
            std::cout <<  " u";
         else
            std::cout <<  " y";
      std::cout << "\n";
   }
std::cout << "Timings By Test Type:\n";
   for (size_t j=0; j < best_timingPA[0].size(); ++j) {
      for (size_t i=0; i<4; ++i) {
         const auto& t = best_timingPA[i][j];
         std::cout << 10.0 * t.time << " " << t.table_bits << " ";
         if (t.code_section < 10)
            std::cout << "0";
         std::cout << t.code_section;
         if (t.uses_sliding_window)
            std::cout <<  " t";
         else
            std::cout <<  " x";
         if (t.uses_squaring_values)
            std::cout <<  "t";
         else
            std::cout <<  "x";
         if (t.is_low_uops)
            std::cout <<  "u";
         else
            std::cout <<  "y";
         if (t.array_size < 10)
            std::cout <<  " 0" << t.array_size;
         else
            std::cout <<  " " << t.array_size;
         if (i != 3)
            std:: cout << "   ";
      }
      std::cout << "\n";
   }
#endif





#if defined(TEST_ARRAY)
   std::cout << "\nbegin benchmarks - array pow\n";

   // warm up call before benchmarking
//   bench_array_pow<4, 0, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, max_modulus_bits_reduce, seed, exponent_bits_reduce);

      // format is bench_array_pow<TABLE_BITS, CODE_SECTION, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION>(...)

   std::array<std::array<std::vector<TimingA>, NUM_TEST_REPETITIONS>, 4> timingA;

   for (size_t i=0; i<4; ++i) {
     for (size_t j=0; j<timingA[i].size(); ++j) {

      timingA[i][j].push_back(
         bench_array_pow<2, 0, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<3, 0, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<4, 0, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<5, 0, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


      timingA[i][j].push_back(
         bench_array_pow<2, 1, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<3, 1, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<4, 1, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<5, 1, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 0, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<3, 0, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 0, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<4, 0, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 0, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<5, 0, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 0, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


      timingA[i][j].push_back(
         bench_array_pow<2, 1, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<2, 1, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<3, 1, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<3, 1, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<4, 1, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<4, 1, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_pow<5, 1, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_pow<5, 1, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
   }
     }
   }
#ifdef TEST_CORRECTNESS_ONLY
   std::cout << "no errors found" << "\n\n";
   return 0;
#endif

   std::array<std::vector<TimingA>, 4> best_timingA;
   for (size_t i=0; i<4; ++i) {
      best_timingA[i] = timingA[i][0];

      for (const auto& tvec : timingA[i]) {
         for (size_t j=0; j<tvec.size(); ++j) {
            if (best_timingA[i][j].time > tvec[j].time)
               best_timingA[i][j].time = tvec[j].time;
         }
      }
   }

   std::vector<TimingA> overall_bestA = best_timingA[0];
   for (size_t j=0; j<overall_bestA.size(); ++j) {
      for (size_t i=1; i<4; ++i) {
         overall_bestA[j].time += best_timingA[i][j].time;
      }
   }
   std::sort(overall_bestA.begin(), overall_bestA.end(), [](TimingA const& t1, TimingA const& t2) { return t1.time < t2.time; });

   for (size_t i=0; i<4; ++i) {
      std::sort(best_timingA[i].begin(), best_timingA[i].end(), [](TimingA const& t1, TimingA const& t2) { return t1.time < t2.time; });
   }


   std::cout << "(ignore)" << uint_to_string(dummy) << "\n\n";

std::cout << "OVERALL BEST:\n";
   for (size_t j=0; j < overall_bestA.size(); ++j) {
         const auto& t = overall_bestA[j];
         std::cout << 10.0 * t.time << "  " << t.table_bits << " ";
         if (t.code_section < 10)
            std::cout << "0";
         std::cout << t.code_section;
         if (t.uses_squaring_values)
            std::cout <<  " t";
         else
            std::cout <<  " x";
         if (t.array_size < 10)
            std::cout <<  " 0" << t.array_size;
         else
            std::cout <<  " " << t.array_size;
      std::cout << "\n";
   }
std::cout << "Timings By Test Type:\n";
   for (size_t j=0; j < best_timingA[0].size(); ++j) {
      for (size_t i=0; i<4; ++i) {
         const auto& t = best_timingA[i][j];
         std::cout << 10.0 * t.time << "  " << t.table_bits << " ";
         if (t.code_section < 10)
            std::cout << "0";
         std::cout << t.code_section;
         if (t.uses_squaring_values)
            std::cout <<  " t";
         else
            std::cout <<  " x";
         if (t.array_size < 10)
            std::cout <<  " 0" << t.array_size;
         else
            std::cout <<  " " << t.array_size;
         if (i != 3)
            std:: cout << "    ";
      }
      std::cout << "\n";
   }
#endif






#if defined(TEST_SCALAR)
   std::cout << "\nbegin benchmarks - scalar pow\n";

   // warm up call before benchmarking
//   bench_range<2, false, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, max_modulus_bits_reduce, seed, exponent_bits_reduce);

   // format is bench_range<TABLE_BITS, USE_SLIDING_WINDOW_OPTIMIZATION, CODE_SECTION, MontType, USE_SQUARING_VALUE_OPTIMIZATION>

   std::array<std::array<std::vector<Timing>, NUM_TEST_REPETITIONS>, 4> timings;

   for (size_t i=0; i<4; ++i) {
     for (size_t j=0; j<timings[i].size(); ++j) {

#if 1
      // Partial array pow using ARRAY_SIZE 1 is essentially a scalar pow.
      // We enter it into the main timing records with (code_section + 50) to distinguish it
      //
      // Note this is somewhat of a hack since we're assuming the tests and number of tests are
      // the same for both bench_range() and bench_partial_array_pow(), which they are for now.
      // If they weren't the same, then entering the PA timings into the main timings would
      // result in invalid timing comparisons/rankings.
      //
      std::vector<TimingPA> vecTimingPA;
      bench_PA_PTAG<hc::LowlatencyTag, 1, MontType>(vecTimingPA, maxU, range, dummy, mmbr[i], seed, ebr[i]);
      for (auto timingPA : vecTimingPA) {
         timings[i][j].push_back(Timing(timingPA.table_bits, timingPA.uses_sliding_window,
               timingPA.code_section + 50, timingPA.time, timingPA.uses_squaring_values));
      }
#endif

#if 1
      timings[i][j].push_back(
         bench_range<0, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<5, false, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5,  true, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<5, false, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5,  true, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 13, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 13, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 14, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 14, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 15, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 15, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 16, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 16, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 17, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 17, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<5, false, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5,  true, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 22, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 22, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 22, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 22, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 22, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 22, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 23, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 23, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 23, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 23, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 23, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 23, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 24, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 24, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 24, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 24, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 25, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 25, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 25, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 25, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 26, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 26, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 27, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 27, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 28, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 28, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 29, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 29, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 30, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 30, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      timings[i][j].push_back(
         bench_range<2, false, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<5, false, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5,  true, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<5, false, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5,  true, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 13, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 13, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 14, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 14, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 15, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 15, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 16, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 16, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 17, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 17, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<5, false, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5,  true, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 22, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 22, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 22, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 22, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 22, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 22, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 23, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 23, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 23, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 23, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, false, 23, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4,  true, 23, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 24, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 24, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 24, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 24, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 25, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 25, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, false, 25, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3,  true, 25, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 26, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 26, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 27, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 27, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 28, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 28, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 29, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 29, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, false, 30, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2,  true, 30, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
   }
#endif
     }
   }
#ifdef TEST_CORRECTNESS_ONLY
   std::cout << "no errors found" << "\n\n";
   return 0;
#endif

   std::array<std::vector<Timing>, 4> best_timings;
   for (size_t i=0; i<4; ++i) {
      best_timings[i] = timings[i][0];

      for (const auto& tvec : timings[i]) {
         for (size_t j=0; j<tvec.size(); ++j) {
            if (best_timings[i][j].time > tvec[j].time)
               best_timings[i][j].time = tvec[j].time;
         }
      }
   }

   std::vector<Timing> overall_best = best_timings[0];
   for (size_t j=0; j<overall_best.size(); ++j) {
      for (size_t i=1; i<4; ++i) {
         overall_best[j].time += best_timings[i][j].time;
      }
   }
   std::sort(overall_best.begin(), overall_best.end(), [](Timing const& t1, Timing const& t2) { return t1.time < t2.time; });

   for (size_t i=0; i<4; ++i) {
      std::sort(best_timings[i].begin(), best_timings[i].end(), [](Timing const& t1, Timing const& t2) { return t1.time < t2.time; });
   }


   std::cout << "(ignore)" << uint_to_string(dummy) << "\n\n";

std::cout << "OVERALL BEST:\n";
   for (size_t j=0; j < overall_best.size(); ++j) {
         const auto& t = overall_best[j];
         std::cout << 10.0 * t.time;
         if (t.uses_sliding_window)
            std::cout <<  "  t";
         else
            std::cout <<  "  x";
         if (t.uses_squaring_values)
            std::cout <<  " t ";
         else
            std::cout <<  " x ";
         std::cout << t.table_bits << " " << t.code_section;
         if (t.code_section < 10)
            std::cout <<  " ";
      std::cout << "\n";
   }
std::cout << "Timings By Test Type:\n";
   for (size_t j=0; j < best_timings[0].size(); ++j) {
      for (size_t i=0; i<4; ++i) {
         const auto& t = best_timings[i][j];
         std::cout << 10.0 * t.time;
         if (t.uses_sliding_window)
            std::cout <<  "  t";
         else
            std::cout <<  "  x";
         if (t.uses_squaring_values)
            std::cout <<  " t ";
         else
            std::cout <<  " x ";
         std::cout << t.table_bits << " " << t.code_section;
         if (t.code_section < 10)
            std::cout <<  " ";
         if (i != 3)
            std:: cout << "    ";
      }
      std::cout << "\n";
   }
#endif



   std::cout << "---Benchmark Program Finished---\n";
   return 0;
}

