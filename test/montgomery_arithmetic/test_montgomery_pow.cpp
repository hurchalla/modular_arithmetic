// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---


// For more complete testing we'll define HURCHALLA_ALLOW_INLINE_ASM_ALL,
// which will cause MontgomeryForm to use any helper inline asm functions that
// are available.  Internally, these inline asm functions will also call their
// corresponding generic template helper functions inside a postcondition, in
// order to make sure that the asm result is correct.  Of course postcondition
// checks must be enabled for this check to occur - the easiest way to ensure
// postconditions are enabled is to undefine NDEBUG, which is why we undef
// NDEBUG here too.
#undef HURCHALLA_ALLOW_INLINE_ASM_ALL
#define HURCHALLA_ALLOW_INLINE_ASM_ALL 1

#undef NDEBUG


//#include "hurchalla/modular_arithmetic/modular_multiplication.h"
//#include "hurchalla/modular_arithmetic/modular_addition.h"
//#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/modular_pow.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/montgomery_pow.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <type_traits>
#include <array>


template <std::size_t NUM_BASES, typename M>
void test_pow_array(M& mf, typename M::T_type base, typename M::T_type exponent)
{
    namespace hc = hurchalla;
    using T = typename M::T_type;
    using V = typename M::MontgomeryValue;

    T modulus = mf.getModulus();

    std::array<T, NUM_BASES> bases;
    std::array<V, NUM_BASES> mv_bases;
    for (std::size_t i=0; i<bases.size(); ++i) {
        bases[i] = static_cast<T>((base + i) % modulus);
        mv_bases[i] = mf.convertIn(bases[i]);
    }

    std::array<V,NUM_BASES> mv_result = mf.pow(mv_bases, exponent);
    for (std::size_t i=0; i<bases.size(); ++i) {
        EXPECT_TRUE(mf.convertOut(mv_result[i]) ==
                           hc::modular_pow<T>(bases[i], exponent, modulus));
    }
}

template <typename M>
void test_pow(M& mf, typename M::T_type base, typename M::T_type exponent)
{
    namespace hc = hurchalla;
    using T = typename M::T_type;

    T modulus = mf.getModulus();

    // first try the non-array overload of pow
    T result = mf.convertOut(mf.pow(mf.convertIn(base), exponent));
    EXPECT_TRUE(result == hc::modular_pow<T>(base, exponent, modulus));

    // then try the array template overload of pow using different array sizes
    test_pow_array<1>(mf, base, exponent);
    test_pow_array<2>(mf, base, exponent);
    test_pow_array<3>(mf, base, exponent);
    test_pow_array<4>(mf, base, exponent);
    test_pow_array<5>(mf, base, exponent);
    test_pow_array<6>(mf, base, exponent);
    test_pow_array<7>(mf, base, exponent);
    test_pow_array<8>(mf, base, exponent);
    test_pow_array<9>(mf, base, exponent);
    test_pow_array<14>(mf, base, exponent);
    test_pow_array<29>(mf, base, exponent);
    test_pow_array<61>(mf, base, exponent);
    test_pow_array<120>(mf, base, exponent);
}


template <typename M>
void run_pow_tests()
{
    using T = typename M::T_type;

    // Try a basic test case first that is valid for all possible Monty types
    // (including even type M == MontySqrtRange<std::uint8_t>).
    {
        T modulus = 13;
        M mf(modulus);
        T base = 6;
        T exponent = 11;
        test_pow(mf, base, exponent);
    }
    // Try a test with the smallest possible modulus
    {
        T modulus = 3;
        M mf(modulus);
        T base = 2;
        T exponent = 5;
        test_pow(mf, base, exponent);
    }
    // Try the largest possible modulus
    {
        T modulus = M::max_modulus();
        M mf(modulus);
        T base = static_cast<T>(modulus - 1);
        T exponent = 179;
        test_pow(mf, base, exponent);
    }

    // Try a bunch of general tests...
    {
        M mf(113);
        T base = 5; T exponent = 6;
        test_pow(mf, base, exponent);
        base = 10;  exponent = 0;
        test_pow(mf, base, exponent);
        base = 0;  exponent = 0;
        test_pow(mf, base, exponent);
        base = 0;  exponent = static_cast<T>(1356);
        test_pow(mf, base, exponent);
        base = 1;  exponent = static_cast<T>(541);
        test_pow(mf, base, exponent);
        base = 67;  exponent = 1;
        test_pow(mf, base, exponent);
        base = 71;  exponent = static_cast<T>(934);
        test_pow(mf, base, exponent);
    }
    {
        T max = M::max_modulus();
        M mf(static_cast<T>(max-2));
        T base = static_cast<T>(max-3); T exponent = 24;
        test_pow(mf, base, exponent);
        base = static_cast<T>(max/2);  exponent = 43;
        test_pow(mf, base, exponent);
        base = static_cast<T>(max/2-1);  exponent = 253;
        test_pow(mf, base, exponent);
        base = 1;  exponent = 135;
        test_pow(mf, base, exponent);
    }
    {
        M mf((M::max_modulus()/4)*2 + 1);
        T base = 5; T exponent = 89;
        test_pow(mf, base, exponent);
        base=static_cast<T>(mf.getModulus()-4); exponent=3;
        test_pow(mf, base, exponent);
        base=static_cast<T>(mf.getModulus()/2); exponent=2;
        test_pow(mf, base, exponent);
        base=static_cast<T>(mf.getModulus()/2-1); exponent=4;
        test_pow(mf, base, exponent);
        base=0; exponent=123;
        test_pow(mf, base, exponent);
    }
}


namespace {
    TEST(MontgomeryArithmetic, montgomery_pow) {
        namespace hc = hurchalla;
        run_pow_tests<hc::MontgomeryFull<std::uint8_t>>();
        run_pow_tests<hc::MontgomeryQuarter<std::uint8_t>>();
        run_pow_tests<hc::MontgomeryStandardMathWrapper<std::uint8_t>>();

        run_pow_tests<hc::MontgomeryFull<std::uint16_t>>();
        run_pow_tests<hc::MontgomeryQuarter<std::uint16_t>>();
        run_pow_tests<hc::MontgomeryStandardMathWrapper<std::uint16_t>>();

        run_pow_tests<hc::MontgomeryFull<std::uint32_t>>();
        run_pow_tests<hc::MontgomeryQuarter<std::uint32_t>>();
        run_pow_tests<hc::MontgomeryStandardMathWrapper<std::uint32_t>>();

        run_pow_tests<hc::MontgomeryFull<std::uint64_t>>();
        run_pow_tests<hc::MontgomeryQuarter<std::uint64_t>>();
        run_pow_tests<hc::MontgomeryStandardMathWrapper<std::uint64_t>>();

#if HURCHALLA_COMPILER_HAS_UINT128_T()
        run_pow_tests<hc::MontgomeryFull<__uint128_t>>();
        run_pow_tests<hc::MontgomeryQuarter<__uint128_t>>();
        run_pow_tests<hc::MontgomeryStandardMathWrapper<__uint128_t>>();
#endif
    }
}
