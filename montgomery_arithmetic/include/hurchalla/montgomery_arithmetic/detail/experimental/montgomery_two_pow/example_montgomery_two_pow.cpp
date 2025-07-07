// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include "montgomery_two_pow.h"
#include "../montgomery_pow_kary.h"
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <utility>


#ifndef NDEBUG
#error "asserts are enabled and will slow performance"
#endif



// this utility function vector_to_array() is adapted from Mikhail's comment at
// https://stackoverflow.com/a/18497366
//
template<class T, std::size_t... Inds>
std::array<T, sizeof...(Inds)> vector_to_array_impl(std::vector<T>& vec, std::integer_sequence<std::size_t, Inds...>)
{
   HPBC_PRECONDITION(vec.size() >= sizeof...(Inds));
   return { vec[Inds]... };
}
template<std::size_t SIZE, class T>
std::array<T, SIZE> vector_to_array(std::vector<T>& vec)
{
   HPBC_PRECONDITION(vec.size() >= SIZE);
   return vector_to_array_impl(vec, std::make_index_sequence<SIZE>{});
}




// Note: uint_to_string() and string_to_uint() provide an easy way to
// do I/O with 128 bit (or larger) integer types.

template <typename U>
std::string uint_to_string(U number)
{
   namespace hc = hurchalla;
   static_assert(hc::ut_numeric_limits<U>::is_integer);
   static_assert(!hc::ut_numeric_limits<U>::is_signed);
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


struct STUException : public std::runtime_error
{
   STUException(std::string const& msg) : std::runtime_error(msg) {}
};

template <typename U>
U string_to_uint(const std::string& str)
{
   namespace hc = hurchalla;
   static_assert(hc::ut_numeric_limits<U>::is_integer);
   static_assert(!hc::ut_numeric_limits<U>::is_signed);
   constexpr U maxU = hc::ut_numeric_limits<U>::max();
   U number = 0;
   for (const auto& c : str) {
      char ch = static_cast<char>(c);
      if (ch < '0' || ch > '9')
         throw STUException("string_to_uint() called with invalid argument:"
               " non-digit character found in 'str'");
      U digit = ch - '0';
      if (number > (maxU - digit) / 10)
         throw STUException("string_to_uint() called with invalid argument:"
               " the contents of 'str' would convert to a value that is too"
               " large to fit in type 'U'");
      number = 10 * number + digit;
   }
   return number;
}





// benchmark a basic simulation of the pow calls in Fermat primality testing

template <class MontType, typename U>
void bench_range(U min, U range)
{
   constexpr auto max_modulus = MontType::max_modulus();

   U max;
   if (range > max_modulus) {
      min = 0;
      max = max_modulus;
   } else {
      // if (min + range > max_modulus) // this might possibly overflow
      if (min > max_modulus - range)
         min = max_modulus - range;
      max = min + range;
   }
   if (max % 2 == 0)
      --max;
   if (min % 2 == 0)
      ++min;
   while ((max - min) % 8 != 0)
      min += 2;

   using namespace std::chrono;
   using dsec = duration<double>;

   dsec::rep mtp_time = 0;
   {
      // the only purpose of total_zeros is to prevent the optimizer from
      // eliminating the function call we want to benchmark in the loop.
      U total_zeros = 0;
      auto t0 = steady_clock::now();

      for (U x = max; x > min; x = x-2) {
         MontType mf(x);
         auto val = hurchalla::montgomery_two_pow(mf, static_cast<U>(x-1));
         // the sole purpose of the next line is to prevent the optimizer from
         // being able to eliminate the function call above.
         if (mf.getCanonicalValue(val) == mf.getZeroValue())
            total_zeros++;
      }

      auto t1 = steady_clock::now();
      mtp_time = dsec(t1-t0).count();
      std::cout << "(ignore) " << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mkp_time = 0;
   {
      // the only purpose of total_zeros is to prevent the optimizer from
      // eliminating the function call we want to benchmark in the loop.
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      for (U x = max; x > min; x = x-2) {
         MontType mf(x);
         auto mont_two = mf.add(mf.getUnityValue(), mf.getUnityValue());
         auto val = hurchalla::montgomery_pow_kary(mf, mont_two, static_cast<U>(x-1));
         if (mf.getCanonicalValue(val) == mf.getZeroValue())
            total_zeros++;
      }
      auto t1 = steady_clock::now();
      mkp_time = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mfp_time = 0;
   {
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      for (U x = max; x > min; x = x-2) {
         MontType mf(x);
//         auto mont_two = mf.convertIn(2);
         auto mont_two = mf.add(mf.getUnityValue(), mf.getUnityValue());
         auto val = mf.pow(mont_two, x-1);
         if (mf.getCanonicalValue(val) == mf.getZeroValue())
            total_zeros++;
      }
      auto t1 = steady_clock::now();
      mfp_time = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mtp_time_2 = 0;
   {
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      for (U x = max; x > min && x >= 4; x = x-4) {
         std::array<MontType, 2> mf_arr { MontType(x), MontType(x - 2) };
         std::array<U, 2> exponent_arr { static_cast<U>(mf_arr[0].getModulus() - 1),  static_cast<U>(mf_arr[1].getModulus() - 1) };
         auto mont_result_arr = hurchalla::array_montgomery_two_pow(mf_arr, exponent_arr);
         if (mf_arr[0].getCanonicalValue(mont_result_arr[0]) == mf_arr[0].getZeroValue())
            total_zeros++;
         if (mf_arr[1].getCanonicalValue(mont_result_arr[1]) == mf_arr[1].getZeroValue())
            total_zeros++;
      }
      auto t1 = steady_clock::now();
      mtp_time_2 = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mtp_time_3 = 0;
   {
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      constexpr std::size_t ARRAY_SIZE = 3;
      for (U x = max; x > min && x >= (2*ARRAY_SIZE); x = x - (2*ARRAY_SIZE)) {
         std::array<MontType, ARRAY_SIZE> mf_arr {
            MontType(x), MontType(x - 2), MontType(x - 4) };
         std::array<U, ARRAY_SIZE> exponent_arr;
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            exponent_arr[j] = mf_arr[j].getModulus() - 1;
         auto mont_result_arr = hurchalla::array_montgomery_two_pow(mf_arr, exponent_arr);
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
            if (mf_arr[j].getCanonicalValue(mont_result_arr[j]) == mf_arr[j].getZeroValue())
               total_zeros++;
         }
      }
      auto t1 = steady_clock::now();
      mtp_time_3 = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mtp_time_4 = 0;
   {
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      constexpr std::size_t ARRAY_SIZE = 4;
      for (U x = max; x > min && x >= (2*ARRAY_SIZE); x = x - (2*ARRAY_SIZE)) {
         std::array<MontType, ARRAY_SIZE> mf_arr {
            MontType(x), MontType(x - 2), MontType(x - 4), MontType(x - 6) };
         std::array<U, ARRAY_SIZE> exponent_arr;
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            exponent_arr[j] = mf_arr[j].getModulus() - 1;
         auto mont_result_arr = hurchalla::array_montgomery_two_pow(mf_arr, exponent_arr);
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
            if (mf_arr[j].getCanonicalValue(mont_result_arr[j]) == mf_arr[j].getZeroValue())
               total_zeros++;
         }
      }
      auto t1 = steady_clock::now();
      mtp_time_4 = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mtp_time_5 = 0;
   {
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      constexpr std::size_t ARRAY_SIZE = 5;
      for (U x = max; x > min && x >= (2*ARRAY_SIZE); x = x - (2*ARRAY_SIZE)) {
         std::array<MontType, ARRAY_SIZE> mf_arr {
            MontType(x), MontType(x - 2), MontType(x - 4), MontType(x - 6), MontType(x - 8) };
         std::array<U, ARRAY_SIZE> exponent_arr;
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            exponent_arr[j] = mf_arr[j].getModulus() - 1;
         auto mont_result_arr = hurchalla::array_montgomery_two_pow(mf_arr, exponent_arr);
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
            if (mf_arr[j].getCanonicalValue(mont_result_arr[j]) == mf_arr[j].getZeroValue())
               total_zeros++;
         }
      }
      auto t1 = steady_clock::now();
      mtp_time_5 = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mtp_time_6 = 0;
   {
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      constexpr std::size_t ARRAY_SIZE = 6;
      for (U x = max; x > min && x >= (2*ARRAY_SIZE); x = x - (2*ARRAY_SIZE)) {
         std::array<MontType, ARRAY_SIZE> mf_arr {
            MontType(x), MontType(x - 2), MontType(x - 4), MontType(x - 6), MontType(x - 8), MontType(x - 10) };
         std::array<U, ARRAY_SIZE> exponent_arr;
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            exponent_arr[j] = mf_arr[j].getModulus() - 1;
         auto mont_result_arr = hurchalla::array_montgomery_two_pow(mf_arr, exponent_arr);
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
            if (mf_arr[j].getCanonicalValue(mont_result_arr[j]) == mf_arr[j].getZeroValue())
               total_zeros++;
         }
      }
      auto t1 = steady_clock::now();
      mtp_time_6 = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mtp_time_8 = 0;
   {
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      constexpr std::size_t ARRAY_SIZE = 8;
      for (U x = max; x > min && x >= (2*ARRAY_SIZE); x = x - (2*ARRAY_SIZE)) {
         std::array<MontType, ARRAY_SIZE> mf_arr {
            MontType(x), MontType(x - 2), MontType(x - 4), MontType(x - 6), MontType(x - 8), MontType(x - 10), MontType(x - 12), MontType(x - 14) };
         std::array<U, ARRAY_SIZE> exponent_arr;
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            exponent_arr[j] = mf_arr[j].getModulus() - 1;
         auto mont_result_arr = hurchalla::array_montgomery_two_pow(mf_arr, exponent_arr);
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
            if (mf_arr[j].getCanonicalValue(mont_result_arr[j]) == mf_arr[j].getZeroValue())
               total_zeros++;
         }
      }
      auto t1 = steady_clock::now();
      mtp_time_8 = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   dsec::rep mpkary_time = 0;
   dsec::rep mfpow_time = 0;
   {
      U total_zeros = 0;
      auto t0 = steady_clock::now();
      constexpr std::size_t ARRAY_SIZE = 4;
      for (U x = max; x > min && x >= (2*ARRAY_SIZE); x = x - (2*ARRAY_SIZE)) {
         MontType mf(x);
         U exponent = mf.getModulus() - 1;
         std::array<typename MontType::MontgomeryValue, ARRAY_SIZE> bases;
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            bases[j] = mf.convertIn(j + 5);
         auto mont_result_arr = hurchalla::array_montgomery_pow_kary(mf, bases, exponent);
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
            if (mf.getCanonicalValue(mont_result_arr[j]) == mf.getZeroValue())
               total_zeros++;
         }
      }
      auto t1 = steady_clock::now();
      mpkary_time = dsec(t1-t0).count();

      t0 = steady_clock::now();
      for (U x = max; x > min && x >= (2*ARRAY_SIZE); x = x - (2*ARRAY_SIZE)) {
         MontType mf(x);
         U exponent = mf.getModulus() - 1;
         std::array<typename MontType::MontgomeryValue, ARRAY_SIZE> bases;
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j)
            bases[j] = mf.convertIn(j + 5);
         auto mont_result_arr = mf.pow(bases, exponent);
         HURCHALLA_REQUEST_UNROLL_LOOP for (int j=0; j<ARRAY_SIZE; ++j) {
            if (mf.getCanonicalValue(mont_result_arr[j]) == mf.getZeroValue())
               total_zeros++;
         }
      }
      t1 = steady_clock::now();
      mfpow_time = dsec(t1-t0).count();
      std::cout << uint_to_string(total_zeros) << " ";
   }

   std::cout << "\n\n";

   std::cout << "montgomery_two_pow() time: " << mtp_time << "\n";
   std::cout << "montgomery_pow_kary() time: " << mkp_time << "\n";
   std::cout << "normal call mf.pow() time: " << mfp_time << "\n";
   std::cout << "array[2]_montgomery_two_pow() time: " << mtp_time_2 << "\n";
   std::cout << "array[3]_montgomery_two_pow() time: " << mtp_time_3 << "\n";
   std::cout << "array[4]_montgomery_two_pow() time: " << mtp_time_4 << "\n";
   std::cout << "array[5]_montgomery_two_pow() time: " << mtp_time_5 << "\n";
   std::cout << "array[6]_montgomery_two_pow() time: " << mtp_time_6 << "\n";
   std::cout << "array[8]_montgomery_two_pow() time: " << mtp_time_8 << "\n";
   std::cout << "performance ratio = " << mfp_time / mtp_time << "\n";
   std::cout << "array2 performance ratio = " << mtp_time / mtp_time_2 << "\n";
   std::cout << "array3 performance ratio = " << mtp_time / mtp_time_3 << "\n";
   std::cout << "array4 performance ratio = " << mtp_time / mtp_time_4 << "\n";
   std::cout << "array5 performance ratio = " << mtp_time / mtp_time_5 << "\n";
   std::cout << "array6 performance ratio = " << mtp_time / mtp_time_6 << "\n";
   std::cout << "array8 performance ratio = " << mtp_time / mtp_time_8 << "\n";
   std::cout << "\narraykary performance ratio = " << mfpow_time / mpkary_time << "\n";

   std::cout << '\n';
}





int main()
{
   namespace hc = hurchalla;
   std::cout << "---Running Example Program---\n\n";

// These are types and values that you may wish to change:
   using U = __uint128_t;
//   using U = uint64_t;

      constexpr int UDIGITS = hc::ut_numeric_limits<U>::digits;
      // Note you're not required to use string_to_uint().  I just used it as a way to set values greater than 2^64 without getting a compile error.
   U exponent = string_to_uint<U>("8");
   U modulus;
   if constexpr (UDIGITS >= 128)
      modulus = string_to_uint<U>("1234567890123456789012345678901");
   else if constexpr (UDIGITS >= 64)
      modulus = string_to_uint<U>("1234567890123456789");
   else if constexpr (UDIGITS >= 32)
      modulus = string_to_uint<U>("123456789");
   else if constexpr (UDIGITS >= 16)
      modulus = string_to_uint<U>("12345");
   else
      modulus = string_to_uint<U>("63");
   if (modulus % 2 == 0) {
      std::cout << "Error: modulus must be odd to use Montgomery arithmetic\n";
      return 1;
   }

#if 1
   // If you can guarantee your modulus will always be less than one quarter the
   // maximum value of type U, then use MontgomeryQuarter for speed.
   using MontType = hc::MontgomeryQuarter<U>;
#else
   // If you can't guarantee your modulus will always be small enough for
   // MontgomeryQuarter, then you must use MontgomeryForm.
   using MontType = hc::MontgomeryForm<U>;
#endif

// demonstration of montgomery_two_pow()
   MontType mf(modulus);
   auto mont_result = hc::montgomery_two_pow(mf, exponent);
   U result = mf.convertOut(mont_result);
   std::cout << "2^" << uint_to_string(exponent) << " (mod " <<
      uint_to_string(modulus) << ") == " << uint_to_string(result) << '\n';

// demonstration of array_montgomery_two_pow(), with an array size of 2.
   // (array_montgomery_two_pow allows you to use any array size > 0.)
   // On my M2 macbook with U = __uint128_t, an array size of 4 benchamrked
   // as fastest per exponentiation, at roughly 1.9x the speed of the plain
   // (non-array) function montgomery_two_pow.
   std::array<MontType, 2> mf_arr { MontType(modulus), MontType(modulus + 2) };  // modulus + 2 is just an arbitrary second modulus value
   std::array<U, 2> exponent_arr { exponent, static_cast<U>(exponent + 3) };  // exponent + 3 is just an arbitrary second exponent value
   auto mont_result_arr = hc::array_montgomery_two_pow(mf_arr, exponent_arr);
   std::array<U, 2> result_arr { mf_arr[0].convertOut(mont_result_arr[0]),
                                 mf_arr[1].convertOut(mont_result_arr[1]) };
   for (int j=0; j<2; ++j) {
      std::cout << "2^" << uint_to_string(exponent_arr[j]) << " (mod "
            << uint_to_string(mf_arr[j].getModulus()) << ") == "
            << uint_to_string(result_arr[j]) << '\n';
   }

   std::cout << '\n';

// ------ End of example portion -------

// Nothing beyong this point is interesting for purposes of an example.
// (You probably don't want to copy anything beyond here, and you don't need to
// read anything beyond here.)





// ------ Tests for correctneess ------

   // test for correctness with a range of exponents
   U range = static_cast<U>(100000);
   constexpr U maxU = hc::ut_numeric_limits<U>::max();
   auto mont_two = mf.add(mf.getUnityValue(), mf.getUnityValue());
   for (exponent = maxU; exponent > maxU-range; exponent = exponent-2) {
      mont_result = hc::montgomery_two_pow(mf, exponent);
      result = mf.convertOut(mont_result);
      U standard_result = mf.convertOut(mf.pow(mont_two, exponent));
      if (result != standard_result) {
         std::cout << "bug in montgomery_two_pow found: got wrong result for ";
         std::cout << "2^" << uint_to_string(exponent) << " (mod " <<
               uint_to_string(modulus) << ")\n";
         return 1;
      }
   }
   for (exponent = maxU; exponent > maxU-range; exponent = exponent-2) {
      constexpr size_t ARRAY_SIZE = 5;
      // We use std::vector to indirectly make a MontType array, since
      // MontType has no default constructor.
      std::vector<MontType> mf_vec;
      std::array<U, ARRAY_SIZE> exponent_arr;
      for (int j=0; j<ARRAY_SIZE; ++j) {
         mf_vec.push_back(mf);
         exponent_arr[j] = exponent + j * 1000000;   // overflow is ok here
      }
      std::array<MontType, ARRAY_SIZE> mf_arr = vector_to_array<ARRAY_SIZE>(mf_vec);

      auto mont_result_arr = hc::array_montgomery_two_pow(mf_arr, exponent_arr);
      for (int j=0; j<ARRAY_SIZE; ++j) {
         result = mf.convertOut(mont_result_arr[j]);
         U standard_result = mf.convertOut(mf.pow(mont_two, exponent_arr[j]));
         if (result != standard_result) {
            std::cout << "bug2 in array_montgomery_two_pow found: got wrong result for ";
            std::cout << "2^" << uint_to_string(exponent_arr[j]) << " (mod " <<
                  uint_to_string(mf.getModulus()) << ")\n";
            return 1;
         }
      }
   }

   // test for correctness with a range of moduli.
   // simulates fermat primality tests
   constexpr auto maxMF = MontType::max_modulus();
   auto mod_range = static_cast<decltype(maxMF)>(range);
   if (mod_range >= maxMF)
      mod_range = maxMF - 1;
   for (auto mod = maxMF; mod > maxMF-mod_range; mod = mod-2) {
      MontType mt(mod);
      auto mont_two = mt.add(mt.getUnityValue(), mt.getUnityValue());
      mont_result = hc::montgomery_two_pow(mt, static_cast<decltype(mod)>(mod-1));
      result = mt.convertOut(mont_result);
      U standard_result = mt.convertOut(mt.pow(mont_two, mod-1));
      if (result != standard_result) {
         std::cout << "bug3 in montgomery_two_pow found: got wrong result for ";
         std::cout << "2^" << uint_to_string(static_cast<decltype(mod)>(mod-1)) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }
   }

   mod_range -= 16;
   for (auto mod = maxMF; mod > maxMF-mod_range; mod = mod-2) {
      constexpr size_t ARRAY_SIZE = 3;
      // We use std::vector to indirectly make a MontType array, since
      // MontType has no default constructor.
      std::vector<MontType> mf_vec;
      std::array<U, ARRAY_SIZE> exponent_arr;
      for (int j=0; j<ARRAY_SIZE; ++j) {
         mf_vec.emplace_back(mod - 2*j);
         exponent_arr[j] = mod + j * 100000;   // overflow is ok here
      }
      std::array<MontType, ARRAY_SIZE> mf_arr = vector_to_array<ARRAY_SIZE>(mf_vec);

      auto mont_result_arr = hc::array_montgomery_two_pow(mf_arr, exponent_arr);
      for (int j=0; j<ARRAY_SIZE; ++j) {
         result = mf_arr[j].convertOut(mont_result_arr[j]);
         auto mont_two = mf_arr[j].add(mf_arr[j].getUnityValue(), mf_arr[j].getUnityValue());
         U standard_result = mf_arr[j].convertOut(mf_arr[j].pow(mont_two, exponent_arr[j]));
         if (result != standard_result) {
            std::cout << "bug4 in array_montgomery_two_pow found: got wrong result for ";
            std::cout << "2^" << uint_to_string(exponent_arr[j]) << " (mod " <<
                  uint_to_string(mf_arr[j].getModulus()) << ")\n";
            return 1;
         }
      }
   }

   std::cout << "All tests succeeded.\n\n";




// ------- Benchmarking --------

   bench_range<MontType>(static_cast<U>(maxU - range), range);

   std::cout << "---Example Program Finished---\n";
   return 0;
}

