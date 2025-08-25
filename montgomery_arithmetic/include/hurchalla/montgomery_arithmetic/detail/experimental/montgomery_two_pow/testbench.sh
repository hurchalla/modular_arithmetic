#!/bin/bash

# Copyright (c) 2025 Jeffrey Hurchalla.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


exit_on_failure () {
  if [ $? -ne 0 ]; then
    exit 1
  fi
}





# Array two_pow  -O3  M2,  no asm,  uint128_t,  MontQuarter, full range modulus (and some noted reduced ranges)
# ------
#  using getRSquaredModN - with gcc appears to add about 6% penalty to array two_pow
#  using getRSquaredModN - with clang appears to add about 3% penalty to array two_pow
# clang: with the above considered (and even without it on most tests), 0 code2 is the definite winner
# gcc:   with the above considered... 0 code0 is the winner (even with the Rsquare penalty) at the very high exponents I'm testing.
#        this remains true with exponent reduced by >> 63.  (presumably the table init overhead is a big hit, just like getRSquaredModN)
#        it remains true when using modulus just above a power of 2 also.
#    ! 0code0 wins, 0 code2 bad
#
# uint128_t,  MontQuarter, modulus reduced
# gcc:   with the above considered... by reducing modulus range >> 2,  0 code0  and 0 code2 are about break-even.
#        with it *not* considered it reaches break-even at about modulus range >> 3 to 5 (up to 5 is needed to break even with exponent? reduced >> 50)
#
# uint128_t,  MontQuarter, using ARM64 asm for two_times_restricted
#   0 code2  wins for gcc and especially clang
#
# CONCLUSION:   For two_pow  use 0 code2, except when using gcc and not using asm (then use 0 code0).
#               In theory for gcc no-asm I could switch to 0 code2 for cases of modulus < (max >> 3),
#                  but because I don't understand why 0 code2 is so slow at full modulus, I generally don't trust 0 code2 for gcc no-asm.


# Array two_pow  -O3  M2,  no asm,  uint128_t,  MontHalfrange, full range modulus (and some noted reduced ranges)
# ------
# clang:  0 code2 pretty clearly wins, regardless of modulus reduction or not.
#         note:  0 code0 wins at array size 2 when ignoring rsquaredmodN, otherwise breaks-even at size 2.
# gcc:    0 code0 wins big-time, regardless of modulus reduction or not.
#         although  asm two_times_restricted doesn't help 0 code2 (which is terrible), that means very little -
#               Halfrange would depend on having HalfrangeCanonical asm to (maybe) get gains.
#
# CONCLUSION:   0 code2 for clang.   0 code0 for gcc.   I don't expect that RsquaredModN will change this.
#               Reevaluate when there's fully enabled asm.


# Array two_pow  -O3  M2,  no asm,  uint128_t,  MontFullrange, full range modulus (and some noted reduced ranges)
# ------
# clang:  0 code0 wins (breaks-even with 0 code2 when considering getRSquaredModN).  5 code0 comes close for array size 7+
#         when considering getRSquaredModN, 5 code0 wins at array size 4+ (array size3 is break-even).
#         asm two_times_restricted doesn't seem to make a difference, but enabling other asm might help for MontFullrange REDC?
# gcc:    0 code0 wins big-time over all others, even considering getRSquaredModN.   Even the best timings are 1.6x worse than clang!!
#         asm two_times_restricted doesn't help - 0 code2 timings are still terrible, presumably caused by no-asm REDC?
#
# uint128_t,  MontQuarter, modulus reduced
# clang:  reducing modulus range >> 2,  situation seems roughly the same as above.
#         with exponent reducted >> 50, relative situation seems roughly unchanged.
# gcc:    reducing modulus range >> 2,  relative winners are roughly the same as above.  note: all timings other than 0 code2 are much improved.
#         reducing modulus range >> 8,  relative winners roughly same:  0 code0 clear winner,  0 code2 still terrible.
#         with exponent reducted >> 50, relative situation unchanged, regardless of modulus range reduction of 0 or 8.
#
# CONCLUSION:   Use 0 code0 for clang and (especially) gcc.  With no-asm and Mont that skips RSquare, reevaluate clang timings.
#               Reevaluate all when there's fully enabled asm.






# Standard two_pow  -O3  M2,  no asm,  uint128_t,  MontQuarter, full range modulus (and some noted reduced ranges)
#  using getRSquaredModN - with gcc appears to add about 5.4% penalty to standard two_pow
#  using getRSquaredModN - with clang appears to add about 4.3% penalty to standard two_pow
#
# clang:  0 false code3  wins,  regardless of rsquaredmodn.  modulus reduction makes no difference.  exponent reduction makes no difference.
#               0.07328  0  false  3
#               0.07397  0  true   2
#               0.07448  0  false  0
#               0.07450  0  true   3
#               0.07479  0  false  1
#               0.07492  0  true   0
#               0.07496  0  false  2
#               0.07616  0  true   1
#               0.07971  5  true   2
#           with modulus reduction (8)
#               0.07346  0  false  3
#               0.07411  0  true   2
#               0.07458  0  false  0
#               0.07468  0  true   3
#               0.07477  0  false  2
#               0.07501  0  false  1
#               0.07518  0  true   0
#               0.07600  0  true   1
#               0.07998  5  true   1
#           with  full modulus and exponent reduction (50)
#               0.04397  0  false  3
#               0.04441  0  true   2
#               0.04480  0  true   3
#               0.04489  0  false  1
#               0.04508  0  true   0
#               0.04518  0  false  2
#               0.04535  0  false  0
#               0.04587  0  true   1
#               0.04924  5  true   2
#           with modulus reduction (8) and exponent reduction (50)
#               0.04394  0  false  3
#               0.04461  0  true   2
#               0.04481  0  true   3
#               0.04490  0  false  1
#               0.04525  0  true   0
#               0.04527  0  false  2
#               0.04531  0  false  0
#               0.04571  0  true   1
#               0.04937  5  true   2
# gcc:  0 false code2  wins,  regardless of rsquaredmodn.
#           with  full modulus and full exponent
#               0.07989  0  false  2
#               0.07999  0  false  0
#               0.08064  0  false  3
#               0.08071  0  true   1
#               0.08160  0  false  1
#               0.08226  0  true   3
#               0.08558  0  true   2
#               0.08689  0  true   0
#               0.09374  4  false  0
#           with modulus reduction (8) and full exponent
#               0.07992  0  false  2
#               0.07995  0  false  0
#               0.08057  0  false  3
#               0.08061  0  true   1
#               0.08153  0  false  1
#               0.08202  0  true   3
#               0.08573  0  true   2
#               0.08714  0  true   0
#               0.09104  1  true   1
#               0.09364  4  false  2
#           with  full modulus and exponent reduction (50)
#               0.04843  0  false  3
#               0.04848  0  true   1
#               0.04854  0  false  2
#               0.04870  0  false  0
#               0.04895  0  false  1
#               0.04945  0  true   3
#               0.05167  0  true   2
#               0.05255  0  true   0
#               0.05864  4  false  1
#           with modulus reduction (8) and exponent reduction (50)
#               0.04808  0  false  3
#               0.04839  0  true   1
#               0.04843  0  false  2
#               0.04863  0  false  0
#               0.04886  0  false  1
#               0.04919  0  true   3
#               0.05166  0  true   2
#               0.05258  0  true   0
#               0.05636  1  true   1
#               0.05844  4  false  2


# Standard two_pow  -O3  M2,  no asm,  uint128_t,  MontHalfrange
# clang:  0 true 2  wins...  but really,  0 false 3  is equally as good.
#           with  full modulus and full exponent
#               0.07593  0  true   2
#               0.07600  0  false  3
#               0.07620  0  true   3
#               0.07624  0  true   0
#               0.07658  0  false  1
#               0.07722  0  false  0
#               0.07761  0  false  2
#               0.07776  0  true   1
#               0.08090  5  true   1
#               0.08100  5  true   2
#           with modulus reduction (8) and full exponent
#               0.07571  0  true   2
#               0.07587  0  false  3
#               0.07614  0  true   0
#               0.07627  0  true   3
#               0.07650  0  false  1
#               0.07704  0  false  0
#               0.07720  0  false  2
#               0.07762  0  true   1
#               0.08063  5  true   1
#               0.08070  5  true   2
#           with  full modulus and exponent reduction (50)
#               0.04525  0  false  3
#               0.04529  0  true   2
#               0.04536  0  true   0
#               0.04556  0  true   3
#               0.04569  0  false  1
#               0.04597  0  false  0
#               0.04608  0  false  2
#               0.04634  0  true   1
#               0.04956  5  true   2
#               0.04963  5  true   1
#           with modulus reduction (8) and exponent reduction (50)
#               0.04529  0  false  3
#               0.04530  0  true   2
#               0.04541  0  true   3
#               0.04545  0  true   0
#               0.04553  0  false  1
#               0.04596  0  false  0
#               0.04610  0  false  2
#               0.04637  0  true   1
#               0.04943  5  true   2
#               0.04951  5  true   1
# gcc:  0 false 2  wins.  0 false 3  is close to just as good.
#           with  full modulus and full exponent
#               0.08142  0  false  3
#               0.08147  0  false  2
#               0.08287  0  false  1
#               0.08404  0  true   1
#               0.08516  0  true   3
#               0.08659  2  true   1
#               0.08661  1  true   1
#               0.08687  3  true   1
#               0.08731  0  false  0
#               0.08807  0  true   2
#               0.09063  0  true   0
#           with modulus reduction (8) and full exponent
#               0.07998  0  false  2
#               0.08067  0  false  3
#               0.08206  0  false  1
#               0.08288  3  true   1
#               0.08290  1  true   1
#               0.08300  0  true   1
#               0.08341  2  true   1
#               0.08409  0  true   3
#               0.08519  0  false  0
#               0.08560  1  true   2
#               0.08685  0  true   2
#               0.08954  0  true   0
#           with  full modulus and exponent reduction (50)
#               0.04877  0  false  2
#               0.04881  0  false  3
#               0.04962  0  false  1
#               0.04997  0  true   1
#               0.05072  0  true   3
#               0.05197  0  false  0
#               0.05249  0  true   2
#               0.05313  2  true   1
#               0.05332  3  true   1
#               0.05333  1  true   1
#               0.05401  0  true   0
#               0.05979  1  true   2
#           with modulus reduction (8) and exponent reduction (50)
#               0.04802  0  false  2
#               0.04841  0  false  3
#               0.04887  0  false  1
#               0.04937  0  true   1
#               0.05023  0  true   3
#               0.05107  0  false  0
#               0.05192  0  true   2
#               0.05198  2  true   1
#               0.05206  1  true   1
#               0.05244  3  true   1
#               0.05338  0  true   0
#               0.05376  1  true   2
#               0.06426  4  false  2


# Standard two_pow  -O3  M2,  no asm,  uint128_t,  MontFullrange
# ------
# clang:   0 false 3  wins.   0 false 2 is definitely not as good.
#           with  full modulus and full exponent
#               0.08671  0  false  3
#               0.08727  0  true   3
#               0.08746  0  true   2
#               0.08752  0  true   0
#               0.08755  0  false  1
#               0.08874  0  true   1
#               0.09084  0  false  0
#               0.09085  0  false  2
#               0.09505  5  true   1
#               0.09530  5  true   2
#               0.09550  5  true   6
#               0.09565  5  true   9
#           with modulus reduction (8) and full exponent
#               0.08677  0  false  3
#               0.08725  0  true   2
#               0.08734  0  true   3
#               0.08748  0  true   0
#               0.08752  0  false  1
#               0.08879  0  true   1
#               0.09067  0  false  2
#               0.09078  0  false  0
#               0.09512  5  true   1
#               0.09537  5  true   6
#               0.09537  5  true   2
#               0.09583  5  true   7
#           with  full modulus and exponent reduction (50)
#               0.05119  0  false  3
#               0.05147  0  true   3
#               0.05156  0  true   0
#               0.05163  0  true   2
#               0.05173  0  false  1
#               0.05203  0  true   1
#               0.05291  0  false  2
#               0.05305  0  false  0
#               0.05683  5  true   1
#               0.05689  5  true   2
#               0.05731  5  true   3
#               0.05731  5  true   4
#           with modulus reduction (8) and exponent reduction (50)
#               0.05106  0  false  3
#               0.05151  0  true   2
#               0.05161  0  true   3
#               0.05162  0  false  1
#               0.05171  0  true   0
#               0.05216  0  true   1
#               0.05294  0  false  0
#               0.05299  0  false  2
#               0.05684  5  true   2
#               0.05697  5  true   1
#               0.05728  5  true   3
# gcc:  0 false 0 wins. (code_section0 ???).  Use 0 false 2 instead.
#           with  full modulus and full exponent
#               0.09623  0  false  0
#               0.09756  0  false  2
#               0.09905  0  false  3
#               0.09911  0  false  1
#               0.10250  0  true   2
#               0.10252  0  true   0
#               0.10402  0  true   3
#               0.10570  0  true   1
#               0.10985  4  false  3
#               0.11014  5  false  0
#           with modulus reduction (8) and full exponent
#               0.08306  0  false  1
#               0.08318  0  false  3
#               0.08318  0  false  0
#               0.08396  0  false  2
#               0.08620  0  true   2
#               0.08636  0  true   0
#               0.08730  0  true   3
#               0.08892  0  true   1
#               0.09564  4  true   3
#               0.09580  4  true   0
#           with  full modulus and exponent reduction (50)
#               0.05734  0  false  0
#               0.05820  0  false  2
#               0.05849  0  false  3
#               0.05882  0  false  1
#               0.06089  0  true   2
#               0.06098  0  true   0
#               0.06151  0  true   3
#               0.06272  0  true   1
#               0.06788  4  false  3
#               0.06792  4  false  2
#           with modulus reduction (8) and exponent reduction (50)
#               0.04879  0  false  3
#               0.04892  0  false  1
#               0.04951  0  false  0
#               0.04968  0  false  2
#               0.05126  0  true   2
#               0.05147  0  true   0
#               0.05151  0  true   3
#               0.05255  0  true   1
#               0.05880  4  false  0
#               0.05881  4  true   2









# Array two_pow  -O3  M2,  no asm,  uint128_t,  MontQuarter
# clang: [update no change]   0 code2 wins, regardless getRSquaredModN
# gcc:   [update: with asm 0 code 2 wins slightly over 0 code0 - but presumably 0 code0 still far better without asm]  0 code0 wins at full mod, 0 code2 bad.  with mod reduce by 7, 0 code 2 10% faster with rsquaredmodn considered
#
# Array two_pow  -O3  M2,  no asm,  uint128_t,  MontHalfrange
# clang:  [update no change]  0 code2 pretty clearly wins
# gcc:    [update no change]  0 code0 wins big-time
#
# Array two_pow  -O3  M2,  no asm,  uint128_t,  MontFullrange
# clang:  [update no change]  0 code0 wins.  when considering getRSquaredModN,  5 code0 wins.
# gcc:    [update no change]  0 code0 wins clearly.  (note: gcc is being pretty badly hobbled all around at full mod - branch prediction I'm guessing)


# Array two_pow  -O2  M2,  no asm,  uint128_t,  MontQuarter
# clang: 0 code2 wins, regardless getRSquaredModN.   actually 0 code0 wins by ~2% when not considering r2modn.
# gcc:   0 code0 wins at full mod, 0 code2 bad.  with mod reduce by 7, 0 code 2 10% faster with rsquaredmodn considered
#
# Array two_pow  -O2  M2,  no asm,  uint128_t,  MontHalfrange
# clang:  0 code2 pretty clearly wins
# gcc:    0 code0 wins big-time
#
# Array two_pow  -O2  M2,  no asm,  uint128_t,  MontFullrange
# clang:  0 code0 wins.  when considering getRSquaredModN,  5 code0 wins.
# gcc:    0 code0 wins clearly.


#gcc   uint64_t cost to make array r2modN is  ~0.00037
#clang uint64_t cost to make array r2modN is  ~0.0003525

# Array two_pow  -O3  M2,  no asm,  uint64_t,  MontQuarter
# clang:  [update 0 code 0 wins, 0 code2 roughly even with it considering r2modn]  0 code0 wins when not considering r2modn,  0 code2 wins when considering it.  pretty even advantage to both, depending on the criteria.
# gcc:    [update no change]   0 code0 wins when not considering r2modn,  0 code2 definitely wins when considering it.
#
# Array two_pow  -O3  M2,  no asm,  uint64_t,  MontHalfrange
# clang:  [update 0 code 0 wins, 0 code2 maybe wins considering r2modn]  0 code0 slight advantage when not considering r2modn; otherwise 0 code2 clear winner.
# gcc:    [update no change]  0 code0 wins when not considering r2modn,  0 code2 wins when considering it.  pretty even advantage to both, depending on the criteria.
#
# Array two_pow  -O3  M2,  no asm,  uint64_t,  MontFullrange
# clang:  [update no change]  0 code0 wins, regardless getRSquaredModN
# gcc:    [update no change]  0 code0 wins, regardless getRSquaredModN
#
# using -O2, the relative standings seem to stay mostly the same.





# Standard two_pow  -O3  M2,  no asm,  uint128_t,  MontQuarter, full range modulus (and some noted reduced ranges)
# clang:  [update 0 true 3 wins, asm or not]  0 false code3  wins,  regardless of rsquaredmodn.  modulus reduction makes no difference.  exponent reduction makes no difference.
# gcc:  [update 0 true 1 wins, asm or not]  0 false code2  wins,  regardless of rsquaredmodn.  0 false 3 not too far behind, though not as good.
#
# Standard two_pow  -O3  M2,  no asm,  uint128_t,  MontHalfrange
# clang:  [update 0 true 3 wins, asm or not]  0 true 2  -or-  0 false 3 , equally as good.
# gcc:  [update 0 false 2 wins, asm or not]   0 false 2  wins.  0 false 3  is close to just as good.
#
# Standard two_pow  -O3  M2,  no asm,  uint128_t,  MontFullrange
# clang:   [update 0 true 3 wins, asm or not]  0 false 3  wins.   0 false 2 is definitely not as good.
# gcc:  [update 0 false 2 without asm, 0 true 3 wins with asm]   0 false 0 wins. (code section 0???).  Use 0 false 2 instead.


#gcc   uint64_t  getRSquaredModN penalty is 0.00061
#clang uint64_t  getRSquaredModN penalty is 0.00055

# Standard two_pow  -O3  M2,  no asm,  uint64_t,  MontQuarter
# clang: [update 0 true 3 wins with asm, 0 true 1 without]   0 true  3 wins, regardless.  0 false 3 is about 1.5% slower.
# gcc:   [update 0 true 3 wins, asm or not]   0 false 3 wins, regardless.
#
# Standard two_pow  -O3  M2,  no asm,  uint64_t,  MontHalf
# clang: [update 0 true 3 wins with asm, 0 true 1 without]   0 true  3 wins, regardless.  0 false 3 is about 2% slower.
# gcc:   [update 0 true 3 wins, with or without asm]   0 false 3 wins, regardless.
#
# Standard two_pow  -O3  M2,  no asm,  uint64_t,  MontFull
# clang: [update 0 true 1 wins, asm or not]   0 true 3 wins, regardless.  0 false 3 about 2.5% slower
# gcc:   [update 0 true 3 wins, asm or not]   0 true 3 wins, regardless.  0 false 3 about 1.3% slower
#
# using -O2, the relative standings seem to stay mostly the same.






#cppcompiler=g++
#cppcompiler=clang++
cppcompiler=$1

#optimization_level=O2
#optimization_level=O3
optimization_level=$2

#define_mont_type=-DPREDEF_MONT_TYPE=MontgomeryQuarter<U>
define_mont_type=-DPREDEF_MONT_TYPE=$3
define_uint_type=-DPREDEF_UINT_TYPE=$4

define_use_asm=$8


cpp_standard=c++17
ndebug=-DNDEBUG


# You need to clone the util, factoring, and modular_arithmetic repos
# from https://github.com/hurchalla


# SET repo_directory TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT
# REPOSITORIES.  (or otherwise ensure the compiler /I flags correctly specify
# the needed hurchalla include directories)
repo_directory=/Users/jeffreyhurchalla/Desktop



if [[ $cppcompiler == "g++" ]]; then
  error_limit=-fmax-errors=3
else
  error_limit=-ferror-limit=3
fi


# argument $8 (if present), should be -DHURCHALLA_ALLOW_INLINE_ASM_ALL



$cppcompiler  \
        $error_limit   -$optimization_level  $ndebug \
        $define_mont_type  $define_uint_type  $define_use_asm \
        -Wall -Wextra -Wpedantic \
        -std=$cpp_standard \
        -I${repo_directory}/modular_arithmetic/modular_arithmetic/include \
        -I${repo_directory}/modular_arithmetic/montgomery_arithmetic/include \
        -I${repo_directory}/util/include \
        -c testbench_montgomery_two_pow.cpp

$cppcompiler  -$optimization_level  -std=$cpp_standard  -o testbench_montgomery_two_pow  testbench_montgomery_two_pow.o -lm

#$cppcompiler  -$optimization_level  -std=$cpp_standard  -o testbench_montgomery_two_pow  sample.o -lm



exit_on_failure

echo "compilation finished, now executing:"


# argument $5 (if present), is the randomization seed for std::mt19937_64
# argument $6 (if present), is max_modulus_bits_reduce
# argument $7 (if present), is exponent_bits_reduce

./testbench_montgomery_two_pow $5 $6 $7

# To give you an example of invoking this script at the command line:
#   ./testbench.sh clang++ O3 MontgomeryFull __uint128_t 191 8 50


