// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

// Test the REDC function versions that contain inline asm.

#undef HURCHALLA_ALLOW_INLINE_ASM_ALL
#define HURCHALLA_ALLOW_INLINE_ASM_ALL

// For extra coverage, we also enable the asserts, so that the internal REDC
// function postconditions call corresponding non-inline asm functions to
// check their results.
#undef HURCHALLA_CLOCKWORK_ENABLE_ASSERTS
#define HURCHALLA_CLOCKWORK_ENABLE_ASSERTS
#undef HURCHALLA_CLOCKWORK_ASSERT_LEVEL
#define HURCHALLA_CLOCKWORK_ASSERT_LEVEL 3


#include "test_REDC.h"


TEST(MontgomeryArithmetic, REDC8_inline_asm) {
    std::vector<uint8_t> moduli { 3, 255, 19, 21, 211, 23, 171 };
    for (auto n : moduli)
        REDC_test_all(n);
}
TEST(MontgomeryArithmetic, REDC16_inline_asm) {
    std::vector<uint16_t> moduli { 3, 17, UINT16_C(65535),
                          UINT16_C(65533), UINT16_C(357), UINT16_C(32253),
                          UINT16_C(11111) };
    for (auto n : moduli)
        REDC_test_all(n);
}
TEST(MontgomeryArithmetic, REDC32_inline_asm) {
    std::vector<uint32_t> moduli { 3, 13, UINT32_C(4294967295),
                          UINT32_C(4294967293), UINT32_C(2147483347),
                          UINT32_C(246098243), UINT32_C(1111111) };
    for (auto n : moduli)
        REDC_test_all(n);
}
TEST(MontgomeryArithmetic, REDC64_inline_asm) {
    std::vector<uint64_t> moduli { 3, 11, UINT64_C(18446744073709551615),
                          UINT64_C(18446744073709551613),
                          UINT64_C(4294967295),
                          UINT64_C(3194806714689), UINT64_C(11111111311) };
    for (auto n : moduli)
        REDC_test_all(n);
}

#if !defined(__GNUC__) || __GNUC__ >= 11 || defined(__INTEL_COMPILER) || \
                                            defined(__clang__)
// Older versions of GCC (most of them prior to v11) have a compiler bug that
// causes an incorrect value of n to be produced and thus results in one of my
// google test assertions failing. See
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98474 .  The bug appears to have
// been introduced as a regression to gcc in v5.1.  It exists up to the latest
// released version (v10.2) of gcc at the time of this writing.  It's unclear
// at the moment whether __uint128_t is safe to use with any version of gcc
// between 5.1 and 10.2.  The patch appears to fix the bug, and it is scheduled
// to be in the gcc 11 release.
// The #if above disables the following tests on gcc prior to gcc v11, since
// they will fail at optimization level -O1 or higher due to the compiler bug.
# if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(MontgomeryArithmetic, REDC128_inline_asm) {
    __uint128_t zero = 0;
    std::vector<__uint128_t> moduli { 3, 11, zero-1, zero-3,
                  static_cast<__uint128_t>(UINT64_C(18446744073709551613)) *
                                             UINT64_C(18446744073709551611),
                  static_cast<__uint128_t>(UINT64_C(35698723439051265)) *
                                                UINT64_C(70945870135873583),
                  static_cast<__uint128_t>(UINT64_C(34069834503)) *
                                              UINT64_C(895835939) };
    for (auto n : moduli)
        REDC_test_all(n);
}
# endif
#endif
