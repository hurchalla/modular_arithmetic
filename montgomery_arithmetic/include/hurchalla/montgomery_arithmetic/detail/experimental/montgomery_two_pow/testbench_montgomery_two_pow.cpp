// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/impl_array_get_Rsquared_mod_n.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/count_leading_zeros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"

#include "experimental_montgomery_two_pow.h"

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



// The function divides a 2U-width dividend by a 1U-width divisor, and
//   produces the quotient and remainder.
// As input precondition, it requires dividend_hi < divisor.
// U can be any unsigned integer type.
//
// This function is adapted from Hacker's Delight 2nd edition, by Henry Warren.
// It is his algorithm of Figure 9-3.
//
// This compiles, superficial results appear correct, but it's basically untested.
// My purpose was possibly to use it to get RSquaredModN, but on M2 it seems to
// have almost exactly the same speed as getRSquaredModN() for __uint128_t.
// I haven't tried any other sizes yet.
// It's unlikely other platforms than M2 will do better, given that M2 has
// great division; so far I expect this function won't be useful to me.
template <typename U>
U div_2U_by_1U(U dividend_hi, U dividend_lo, U divisor, U& remainder)
{
    namespace hc = hurchalla;

    static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
    static_assert(hc::ut_numeric_limits<U>::is_integer, "");

    HPBC_CLOCKWORK_PRECONDITION2(dividend_hi < divisor);

    constexpr int bitsU = hc::ut_numeric_limits<U>::digits;

    const U b = static_cast<U>(1) << (bitsU/2);  // Number base.
    U un1, un0,                     // Norm. dividend LSD’s.
            vn1, vn0,               // Norm. divisor digits.
            q1, q0,                 // Quotient digits.
            un32, un21, un10,       // Dividend digit pairs.
            rhat;                   // A remainder.
    int s;                          // Shift amount for norm.

    HPBC_CLOCKWORK_ASSERT2(divisor > 0);
    s = hc::count_leading_zeros(divisor);
    divisor = divisor << s;         // Normalize divisor.

// note: assuming U is 128bit, vn1 and vn0 can fit in uint64_t
    U mask = (static_cast<U>(1) << (bitsU/2)) - 1;

    vn1 = divisor >> (bitsU/2);     // Break divisor up into
    vn0 = divisor & mask;           // into hi and lo parts

    HPBC_CLOCKWORK_ASSERT2(s < bitsU);
    un32 = (dividend_hi << s) |
           ((dividend_lo >> (bitsU - s - 1)) >> 1);
    un10 = dividend_lo << s;        // Shift dividend left.

    un1 = un10 >> (bitsU/2);        // Break low half of
    un0 = un10 & mask;              // dividend into two parts.

    q1 = un32/vn1;                  // Compute the first
    rhat = un32 - q1*vn1;           // quotient digit, q1.

    while (q1 >= b || q1*vn0 > b*rhat + un1) {
        q1 = q1 - 1;
        rhat = rhat + vn1;
        if (rhat >= b) break;
    }

    un21 = un32*b + un1 - q1*divisor; // Multiply and subtract.

    q0 = un21/vn1;                  // Compute the second
    rhat = un21 - q0*vn1;           // quotient digit, q0.

    while (q0 >= b || q0*vn0 > b*rhat + un0) {
        q0 = q0 - 1;
        rhat = rhat + vn1;
        if (rhat >= b) break;
    }

    remainder = (un21*b + un0 - q0*divisor) >> s;
    return q1*b + q0;
}




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
          size_t CODE_SECTION, class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION>
int test_correctness_two_pow()
{
   namespace hc = hurchalla;
   using U = typename MontType::IntegerType;
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");

   constexpr U maxU = hc::ut_numeric_limits<U>::max();
   U range = static_cast<U>(100);

   constexpr U maxMF = MontType::max_modulus();
   auto mod_range = static_cast<U>(range);
   if (mod_range >= maxMF)
      mod_range = maxMF - 1;

   for (U i = 0; i < mod_range - 2; ++i) {
      U mod = (i % 2 == 0) ? maxMF - i : i + 2;
      MontType mf(mod);
      auto mont_two = mf.add(mf.getUnityValue(), mf.getUnityValue());
      U exponent = (i % 3 == 0) ? i : maxU - i;
      auto mont_result = hc::experimental::experimental_montgomery_two_pow::call<
                MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION>(mf, exponent);
      U result = mf.convertOut(mont_result);
      U standard_result = mf.convertOut(mf.pow(mont_two, exponent));
      if (result != standard_result) {
         std::cout << "bug in montgomery_two_pow found: got wrong result for ";
         std::cout << "2^" << uint_to_string(exponent) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }
   }
   // simulate fermat primality tests
   for (U i = 0; i < mod_range - 2; ++i) {
      U mod = (i % 2 == 0) ? maxMF - i : i + 2;
      MontType mt(mod);
      auto mont_two = mt.add(mt.getUnityValue(), mt.getUnityValue());
      auto mont_result = hc::experimental::experimental_montgomery_two_pow::call<
                MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION>(mt, static_cast<decltype(mod)>(mod-1));
      U result = mt.convertOut(mont_result);
      U standard_result = mt.convertOut(mt.pow(mont_two, mod-1));
      if (result != standard_result) {
         std::cout << "bug3 in montgomery_two_pow found: got wrong result for ";
         std::cout << "2^" << uint_to_string(static_cast<decltype(mod)>(mod-1)) << " (mod " <<
               uint_to_string(mod) << ")\n";
         return 1;
      }
   }
   return 0;
}


template <size_t TABLE_BITS, size_t CODE_SECTION, size_t ARRAY_SIZE, class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION>
int test_correctness_array_two_pow()
{
   namespace hc = hurchalla;
   using U = typename MontType::IntegerType;
   static_assert(!hc::ut_numeric_limits<U>::is_signed, "");
   static_assert(hc::ut_numeric_limits<U>::is_integer, "");

   constexpr U range = static_cast<U>(100);

   constexpr U maxMF = MontType::max_modulus();
   auto mod_range = static_cast<U>(range);
   if (mod_range >= maxMF)
      mod_range = maxMF - 1;

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

      auto mont_result_arr = hc::experimental::experimental_montgomery_two_pow::call
                  <MontType, U, ARRAY_SIZE, TABLE_BITS, CODE_SECTION, USE_SQUARING_VALUE_OPTIMIZATION>
                  (mf_arr, exponent_arr);

      for (size_t j=0; j<ARRAY_SIZE; ++j) {
         U result = mf_arr[j].convertOut(mont_result_arr[j]);
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
   return 0;
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
bench_array_two_pow(U min, U range, U& totalU, unsigned int max_modulus_bits_reduce, ST seed, unsigned int exponent_bits_reduce)
{
   HPBC_CLOCKWORK_PRECONDITION2(max_modulus_bits_reduce <
                     hurchalla::ut_numeric_limits<decltype(MontType::max_modulus())>::digits);

   // run very short tests to hopefully catch a bugged experimental impl
   int tcatp_result = test_correctness_array_two_pow<TABLE_BITS,
                                                     CODE_SECTION,
                                                     ARRAY_SIZE,
                                                     MontType,
                                                     USE_SQUARING_VALUE_OPTIMIZATION>();
   if (tcatp_result != 0) {
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

   using namespace std::chrono;
   using dsec = duration<double>;

   constexpr auto max_modulus = MontType::max_modulus();

   // I initialize mf using x (which can be as high as maxMod), so
   // x (and therefore maxMod) will need to be <= max_modulus.

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
      for (size_t i=0; i<tmpvec.size(); ++i) {
         U x = tmpvec[i];
         MontType mf(x);
         mfvec.push_back(mf);
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



      using V = typename MontType::MontgomeryValue;

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
            U u_hi = hurchalla::unsigned_multiply_to_hilo_product(u_lo, r_mod_n[j], r_mod_n[j]);
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
         for (size_t j=0; j < ARRAY_SIZE; ++j)
            exparr[j] = randexpU[i + j];

         std::array<V, ARRAY_SIZE> result = hurchalla::experimental::experimental_montgomery_two_pow::call
               <MontType, U, ARRAY_SIZE, TABLE_BITS, CODE_SECTION, USE_SQUARING_VALUE_OPTIMIZATION>(mfarr, exparr);

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
      //std::cout << mtp_time << "     (ignore " << uint_to_string(totalU) << ")\n";

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


// benchmark a basic simulation of the pow calls in Fermat primality testing

template <size_t TABLE_BITS, bool USE_SLIDING_WINDOW_OPTIMIZATION,
          size_t CODE_SECTION, class MontType, bool USE_SQUARING_VALUE_OPTIMIZATION,
          typename U, typename ST>
Timing
bench_range(U min, U range, U& totalU, unsigned int max_modulus_bits_reduce, ST seed, unsigned int exponent_bits_reduce)
{
   HPBC_CLOCKWORK_PRECONDITION2(max_modulus_bits_reduce <
                     hurchalla::ut_numeric_limits<decltype(MontType::max_modulus())>::digits);

   // run very short tests to hopefully catch a bugged experimental impl
   int tctp_result = test_correctness_two_pow<TABLE_BITS,
                                             USE_SLIDING_WINDOW_OPTIMIZATION,
                                             CODE_SECTION,
                                             MontType,
                                             USE_SQUARING_VALUE_OPTIMIZATION>();
   if (tctp_result != 0) {
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

   using namespace std::chrono;
   using dsec = duration<double>;

   constexpr auto max_modulus = MontType::max_modulus();

   // I initialize mf using x (which can be as high as maxMod), so
   // x (and therefore maxMod) will need to be <= max_modulus.

   constexpr bool randomizeModuli = true;
   unsigned int exponentreduction = exponent_bits_reduce; // 2; // 2 // 50
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
      std::vector<MontType> mfvec;
      for (size_t i=0; i<tmpvec.size(); ++i) {
         U x = tmpvec[i];
         MontType mf(x);
         mfvec.push_back(mf);
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
#if 0
// for comparison...
         auto two = mf.two_times(mf.getUnityValue());
         auto val = mf.pow(two, exponent);
#else
         auto val = hurchalla::experimental::experimental_montgomery_two_pow::call<
                MontType, U, USE_SLIDING_WINDOW_OPTIMIZATION, TABLE_BITS, CODE_SECTION,
                USE_SQUARING_VALUE_OPTIMIZATION>(mf, static_cast<U>(exponent));
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
      //std::cout << mtp_time << "     (ignore " << uint_to_string(totalU) << ")\n";

//      return mtp_time;
      return Timing(TABLE_BITS, USE_SLIDING_WINDOW_OPTIMIZATION, CODE_SECTION, mtp_time, USE_SQUARING_VALUE_OPTIMIZATION);
   }
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



#ifndef DEF_UINT_TYPE
#  error "DEF_UINT_TYPE was not defined"
#else
   using U = DEF_UINT_TYPE;
#endif
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!(ut_numeric_limits<U>::is_signed), "");

//#define XSTR(x) STR(x)
//#define STR(x) #x
//#pragma message "The value of DEF_UINT_TYPE: " XSTR(DEF_UINT_TYPE)


/*
#if 1
   using U = __uint128_t;
#else
   using U = uint64_t;
#endif
*/

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



   constexpr int UDIGITS = hc::ut_numeric_limits<U>::digits;
      // Note you're not required to use string_to_uint().  I just used it as a way to set values greater than 2^64 without getting a compile error.

   U modulus;
   if (UDIGITS >= 128)
      modulus = string_to_uint<U>("1234567890123456789012345678901");
   else if (UDIGITS >= 64)
      modulus = string_to_uint<U>("1234567890123456789");
   else if (UDIGITS >= 32)
      modulus = string_to_uint<U>("123456789");
   else if (UDIGITS >= 16)
      modulus = string_to_uint<U>("12345");
   else
      modulus = string_to_uint<U>("63");

   if (modulus % 2 == 0) {
      std::cout << "Error: modulus must be odd to use Montgomery arithmetic\n";
      return 1;
   }


   constexpr U maxU = hc::ut_numeric_limits<U>::max();
   U range = static_cast<U>(100000);




// ------- Benchmarking --------

//   constexpr size_t TABLE_BITS = 0; // 1; // 0; // 3;
//   constexpr bool USE_SLIDING_WINDOW_OPTIMIZATION = true;
//   constexpr size_t CODE_SECTION = 0; // 1; // 0; // 1;

   std::cout << std::fixed;
   std::cout.precision(4);

   auto seed = randomization_seed;

   U dummy = 0; // dummy exists to prevent the compiler from optimizing away getting timings

   std::array<unsigned int, 4> mmbr = { 0, max_modulus_bits_reduce, 0, max_modulus_bits_reduce };
   unsigned int default_ebr = (std::is_same<MontType, MontgomeryQuarter<U>>::value) ? 2 :
                                 (std::is_same<MontType, MontgomeryHalf<U>>::value) ? 1 : 0;
   std::array<unsigned int, 4> ebr = { default_ebr, default_ebr, exponent_bits_reduce, exponent_bits_reduce };


   constexpr int NUM_TEST_REPETITIONS = 2;


#if 1
   std::cout << "\nbegin benchmarks - array two_pow\n";

   // warm up call
   bench_array_two_pow<5, 8, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, max_modulus_bits_reduce, seed, exponent_bits_reduce);

      // format is bench_array_two_pow<TABLE_BITS, CODE_SECTION, ARRAY_SIZE, MontType, USE_SQUARING_VALUE_OPTIMIZATION>(...)

   std::array<std::array<std::vector<TimingA>, NUM_TEST_REPETITIONS>, 4> timingA;

   for (size_t i=0; i<4; ++i) {
     for (size_t j=0; j<timingA[i].size(); ++j) {

#if 1
      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));



      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));



   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 27, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 28, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 29, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
   }
# endif

# if 1
      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
   }
# endif

# if 1
      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 2, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 3, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 10, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 11, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<0, 0, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 1, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 4, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 5, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 6, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<0, 7, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<3, 8, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<4, 8, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<5, 8, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<6, 8, 12, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timingA[i][j].push_back(
         bench_array_two_pow<1, 8, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<2, 8, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timingA[i][j].push_back(
         bench_array_two_pow<7, 8, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
   }
# endif
     }
   }


   std::array<std::vector<TimingA>, 4> best_timingA;
   for (size_t i=0; i<4; ++i) {
      best_timingA[i] = timingA[i][0];

      for (const auto& tvec : timingA[i]) {
         for (size_t j=0; j<tvec.size(); ++j) {
            if (best_timingA[i][j].time > tvec[j].time)
               best_timingA[i][j].time = tvec[j].time;
         }
      }
      std::sort(best_timingA[i].begin(), best_timingA[i].end(), [](TimingA const& t1, TimingA const& t2) { return t1.time < t2.time; });
   }

   std::cout << "(ignore)" << uint_to_string(dummy) << "\n\n";

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






#if 1
   std::cout << "\nbegin benchmarks - scalar two_pow\n";

   //  warm up to get cpu boost (or throttle) going
   for (size_t i=0; i<4; ++i)
      bench_range<0, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, max_modulus_bits_reduce, seed, exponent_bits_reduce);

//   std::array<std::vector<Timing>, 4> timings;

   std::array<std::array<std::vector<Timing>, NUM_TEST_REPETITIONS>, 4> timings;

   for (size_t i=0; i<4; ++i) {
     for (size_t j=0; j<timings[i].size(); ++j) {

#if 0
      timings[i][j].push_back(
         bench_range<0, true , 3, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 3, MontType, false >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 3, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 3, MontType, false >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
#endif
#if 0
   // this repeats calls below - it is the best performers on clang for MontyFull __uint128_t when
   // considering only squaring value optimization
   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      timings[i][j].push_back(
         bench_range<0, true , 0, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 0, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 1, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 1, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 2, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 2, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
//      timings[i][j].push_back(
//         bench_range<0, true , 3, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, false, 15, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 16, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 15, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 16, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 17, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
   }
#endif

#if 1
      timings[i][j].push_back(
         bench_range<0, true , 17, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 17, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 18, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 18, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 22, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 22, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 23, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 23, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 24, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 24, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 25, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 25, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 26, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 26, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 27, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 27, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, false, 28, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 29, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
         timings[i][j].push_back(
         bench_range<0, false, 30, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 19, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 19, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 20, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 20, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 21, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, false, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, false, 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 13, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 10, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 11, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 12, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 13, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, false, 14, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 15, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 16, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 14, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 15, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 16, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      timings[i][j].push_back(
         bench_range<0, true , 17, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 17, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 18, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 18, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 22, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 22, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 23, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 23, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 24, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 24, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 25, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 25, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 26, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 26, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 27, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 27, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, false, 28, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 29, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 19, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 19, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 20, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 20, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 21, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 0, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 0, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 1, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 1, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, true , 2, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 2, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 3, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 3, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, false, 10, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 11, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 12, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 13, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 10, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 11, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 12, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 13, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<0, false, 14, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 15, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, false, 16, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 14, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 15, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<0, true , 16, MontType, true >(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
   }

      timings[i][j].push_back(
         bench_range<1, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<1, true , 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<1, true , 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<1, true , 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2, true , 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3, true , 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
#endif
# if 1
      timings[i][j].push_back(
         bench_range<1, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<1, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, true , 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, true , 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, false, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, false, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, false, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<5, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 3, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 4, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 5, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 6, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 7, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 8, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 9, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<6, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, true , 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, false, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, false, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<7, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, true , 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, true , 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, false, 1, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, false, 2, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<8, true , 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<8, false, 0, MontType, false>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));


   if constexpr (std::is_same<typename MontType::MontType::MontyTag,
                                    ::hurchalla::detail::TagMontyFullrange>::value) {
      timings[i][j].push_back(
         bench_range<1, true , 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<1, false, 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<2, true , 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<2, false, 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<3, true , 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<3, false, 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<4, true , 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, true , 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, true , 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, true , 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, false, 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, false, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, false, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<4, false, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<5, true , 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, true , 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 3, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 4, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 5, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 6, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 7, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 8, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<5, false, 9, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<6, true , 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, true , 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, true , 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, false, 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, false, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<6, false, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<7, true , 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, true , 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, true , 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, false, 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, false, 1, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<7, false, 2, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));

      timings[i][j].push_back(
         bench_range<8, true , 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
      timings[i][j].push_back(
         bench_range<8, false, 0, MontType, true>(static_cast<U>(maxU - range), range, dummy, mmbr[i], seed, ebr[i]));
   }
# endif
     }
   }


   std::array<std::vector<Timing>, 4> best_timings;
   for (size_t i=0; i<4; ++i) {
      best_timings[i] = timings[i][0];

      for (const auto& tvec : timings[i]) {
         for (size_t j=0; j<tvec.size(); ++j) {
            if (best_timings[i][j].time > tvec[j].time)
               best_timings[i][j].time = tvec[j].time;
         }
      }
      std::sort(best_timings[i].begin(), best_timings[i].end(), [](Timing const& t1, Timing const& t2) { return t1.time < t2.time; });
   }

   std::cout << "(ignore)" << uint_to_string(dummy) << "\n\n";

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



   std::cout << "---Example Program Finished---\n";
   return 0;
}

