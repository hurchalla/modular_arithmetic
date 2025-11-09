
// -------
// This file is intended to be included multiple times, using different macro definitions
// -------

// The macros to set:

// I was not using this
// HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLESIZE_INIT
//    Set this as blank, or to HURCHALLA_REQUEST_UNROLL_LOOP

// I was using this
// HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE
//    Set this as blank, or to HURCHALLA_REQUEST_UNROLL_LOOP

// for 128bit I was NOT unrolling on NUM_TABLES in table init, but I was for 64 bit.
// HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_NUM_TABLES_INIT
//    Set this as blank, or to HURCHALLA_REQUEST_UNROLL_LOOP

// I was using this
// HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLE_BITS
//    Set this as blank, or to HURCHALLA_REQUEST_UNROLL_LOOP

// I was not using this
// HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_NUM_TABLES_MAINLOOP
//    Set this as blank, or to HURCHALLA_REQUEST_UNROLL_LOOP



// For now, I'm not going to utlize these macros, but I could:
// I was using this
// HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLE_BITS_ENDING
//    Set this as blank, or to HURCHALLA_REQUEST_UNROLL_LOOP





    std::array<std::array<std::array<V, ARRAY_SIZE>, TABLESIZE>, NUM_TABLES> table;

    V mont_one = mf.getUnityValue();

    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q) {
        table[0][0][q] = mont_one;
        table[0][1][q] = x[q];
    }
    if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
        HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
            table[0][2][q] = mf.square(x[q]);
        HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
            table[0][3][q] = mf.template multiply<PTAG>(table[0][2][q], x[q]);
    }
    if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
        HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLESIZE_INIT for (size_t i=4; i<TABLESIZE; i+=2) {
            size_t j = i/2;
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                table[0][i][q] = mf.template square<LowuopsTag>(table[0][j][q]);
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                table[0][i+1][q] = mf.template multiply<LowuopsTag>(table[0][j+1][q], table[0][j][q]);
        }
    }


    U n_orig = n;
    (void)n_orig;   // silence potential unsed var warnings
    int shift;
    size_t tmp;
    if (n > MASKBIG) {
        HPBC_CLOCKWORK_ASSERT2(n > 0);
        int leading_zeros = count_leading_zeros(n);
        int numbits = ut_numeric_limits<decltype(n)>::digits - leading_zeros;
        HPBC_CLOCKWORK_ASSERT2(numbits > NUMBITS_MASKBIG);
        shift = numbits - NUMBITS_MASKBIG;
        HPBC_CLOCKWORK_ASSERT2(shift > 0);
        tmp = static_cast<size_t>(branchless_shift_right(n, shift));
        // this preps n ahead of time for the main loop
        n = branchless_shift_left(n, leading_zeros + NUMBITS_MASKBIG);
    }
    else {
        shift = 0;
        tmp = static_cast<size_t>(n);
    }
    HPBC_CLOCKWORK_ASSERT2(shift >= 0);

    HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);


    std::array<V, ARRAY_SIZE> result;
    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q) {
        result[q] = table[0][tmp & MASK][q];
    }


//    constexpr int digitsRU = hurchalla::ut_numeric_limits<typename MFE_LU::RU>::digits;

    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_NUM_TABLES_INIT for (size_t k=1; k < NUM_TABLES; ++k) {
        if HURCHALLA_CPP17_CONSTEXPR (UseEarlyExitInInit) {
            // this part could be removed - it provides fast return when n_orig is small.
            size_t limit_in_progress = static_cast<size_t>(1) << (k * TABLE_BITS);
            if (n_orig < limit_in_progress)
                return result;
        }
        HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q) {
            table[k][0][q] = mont_one;
            table[k][1][q] = mf.square(table[k-1][TABLESIZE/2][q]);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE >= 4) {
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                table[k][2][q] = mf.square(table[k][1][q]);
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                table[k][3][q] = mf.template multiply<PTAG>(table[k][2][q], table[k][1][q]);
        }
        if HURCHALLA_CPP17_CONSTEXPR (TABLESIZE > 4) {
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLESIZE_INIT for (size_t i=4; i<TABLESIZE; i+=2) {
                size_t j = i/2;
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    table[k][i][q] = mf.template square<LowuopsTag>(table[k][j][q]);
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    table[k][i+1][q] = mf.template multiply<LowuopsTag>(table[k][j+1][q], table[k][j][q]);
            }
        }

        size_t index = (tmp >> (k * TABLE_BITS)) & MASK;
        HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
            result[q] = mf.template multiply<LowuopsTag>(table[k][index][q], result[q]);
    }
    int bits_remaining = shift;


    while (bits_remaining >= NUMBITS_MASKBIG) {
        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            SV sv[ARRAY_SIZE];
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                sv[q] = MFE_LU::getSquaringValue(mf, result[q]);
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (bits_remaining > NUMBITS_MASKBIG &&
                                (static_cast<size_t>(n >> high_word_shift) &
                                     (static_cast<size_t>(1) << (digits_smaller - 1))) == 0) {
                    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                        sv[q] = MFE_LU::squareSV(mf, sv[q]);
                    n = static_cast<U>(n << 1);
                    --bits_remaining;
                }
            }
            HPBC_CLOCKWORK_ASSERT2(bits_remaining >= NUMBITS_MASKBIG);

            tmp = static_cast<size_t>(n >> high_word_shift) >> small_shift;
            n = static_cast<U>(n << NUMBITS_MASKBIG);
            bits_remaining -= NUMBITS_MASKBIG;

            V val1[ARRAY_SIZE];
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                val1[q] = table[0][tmp & MASK][q];

            static_assert(TABLE_BITS >= 1, "");
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLE_BITS for (size_t i=0; i<TABLE_BITS - 1; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    sv[q] = MFE_LU::squareSV(mf, sv[q]);
            }

            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_NUM_TABLES_MAINLOOP for (size_t k=1; k<NUM_TABLES; ++k) {
                tmp = tmp >> TABLE_BITS;
                size_t index = tmp & MASK;
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    val1[q] = mf.template multiply<LowuopsTag>(val1[q], table[k][index][q]);

                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLE_BITS for (size_t i=0; i<TABLE_BITS; ++i) {
                    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                        sv[q] = MFE_LU::squareSV(mf, sv[q]);
                }
            }
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                result[q] = MFE_LU::squareToMontgomeryValue(mf, sv[q]);

            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                result[q] = mf.template multiply<PTAG>(result[q], val1[q]);
        }
        else {
            if HURCHALLA_CPP17_CONSTEXPR (USE_SLIDING_WINDOW_OPTIMIZATION) {
                while (bits_remaining > NUMBITS_MASKBIG &&
                                (static_cast<size_t>(n >> high_word_shift) &
                                     (static_cast<size_t>(1) << (digits_smaller - 1))) == 0) {
                    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                        result[q] = mf.template square<PTAG>(result[q]);
                    n = static_cast<U>(n << 1);
                    --bits_remaining;
                }
            }
            HPBC_CLOCKWORK_ASSERT2(bits_remaining >= NUMBITS_MASKBIG);

            tmp = static_cast<size_t>(n >> high_word_shift) >> small_shift;
            n = static_cast<U>(n << NUMBITS_MASKBIG);
            bits_remaining -= NUMBITS_MASKBIG;

            V val1[ARRAY_SIZE];
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                val1[q] = table[0][tmp & MASK][q];

            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLE_BITS for (size_t i=0; i<TABLE_BITS; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    result[q] = mf.template square<PTAG>(result[q]);
            }

            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_NUM_TABLES_MAINLOOP for (size_t k=1; k<NUM_TABLES; ++k) {
                tmp = tmp >> TABLE_BITS;
                size_t index = tmp & MASK;
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    val1[q] = mf.template multiply<LowuopsTag>(val1[q], table[k][index][q]);

                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLE_BITS for (size_t i=0; i<TABLE_BITS; ++i)
                    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                        result[q] = mf.template square<PTAG>(result[q]);
            }

            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                result[q] = mf.template multiply<PTAG>(result[q], val1[q]);
        }
    }
    if (bits_remaining == 0)
        return result;

    HPBC_CLOCKWORK_ASSERT2(0 < bits_remaining && bits_remaining < NUMBITS_MASKBIG);

    tmp = static_cast<size_t>(n >> high_word_shift) >> (digits_smaller - bits_remaining);
    HPBC_CLOCKWORK_ASSERT2(tmp <= MASKBIG);

    V val1[ARRAY_SIZE];
    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
        val1[q] = table[0][tmp & MASK][q];


    if HURCHALLA_CPP17_CONSTEXPR (NUM_TABLES <= 2) {
        // here we only handle NUM_TABLES <= 2 because when we have larger
        // numbers of tables we optimize for that below
        HURCHALLA_REQUEST_UNROLL_LOOP for (size_t k=1; k<NUM_TABLES; ++k) {
            size_t index = (tmp >> (k * TABLE_BITS)) & MASK;
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                val1[q] = mf.template multiply<PTAG>(val1[q], table[k][index][q]);
        }
        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            SV sv[ARRAY_SIZE];
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                sv[q] = MFE_LU::getSquaringValue(mf, result[q]);
            HPBC_CLOCKWORK_ASSERT2(bits_remaining >= 1);
            for (int i=0; i<bits_remaining-1; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    sv[q] = MFE_LU::squareSV(mf, sv[q]);
            }
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                result[q] = MFE_LU::squareToMontgomeryValue(mf, sv[q]);
        }
        else {
            for (int i=0; i<bits_remaining; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    result[q] = mf.template square<PTAG>(result[q]);
            }
        }
    }
    else {
        if HURCHALLA_CPP17_CONSTEXPR (USE_SQUARING_VALUE_OPTIMIZATION) {
            SV sv[ARRAY_SIZE];
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                sv[q] = MFE_LU::getSquaringValue(mf, result[q]);
            int i=0;
            for (size_t k=1; i + static_cast<int>(TABLE_BITS) < bits_remaining;
                             i += static_cast<int>(TABLE_BITS), ++k) {
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLE_BITS for (size_t h=0; h<TABLE_BITS; ++h) {
                    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                        sv[q] = MFE_LU::squareSV(mf, sv[q]);
                }
                size_t index = (tmp >> (k * TABLE_BITS)) & MASK;
                HPBC_CLOCKWORK_ASSERT2(k < NUM_TABLES);
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    val1[q] = mf.template multiply<PTAG>(val1[q], table[k][index][q]);
            }
            HPBC_CLOCKWORK_ASSERT2(bits_remaining >= 1);
            HPBC_CLOCKWORK_ASSERT2(i < bits_remaining);
            for (; i<bits_remaining-1; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    sv[q] = MFE_LU::squareSV(mf, sv[q]);
            }
            HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                result[q] = MFE_LU::squareToMontgomeryValue(mf, sv[q]);
        }
        else {
            int i=0;
            for (size_t k=1; i + static_cast<int>(TABLE_BITS) < bits_remaining;
                             i += static_cast<int>(TABLE_BITS), ++k) {
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_TABLE_BITS for (size_t h=0; h<TABLE_BITS; ++h) {
                    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                        result[q] = mf.template square<PTAG>(result[q]);
                }
                size_t index = (tmp >> (k * TABLE_BITS)) & MASK;
                HPBC_CLOCKWORK_ASSERT2(k < NUM_TABLES);
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    val1[q] = mf.template multiply<PTAG>(val1[q], table[k][index][q]);
            }
            for (; i<bits_remaining; ++i) {
                HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
                    result[q] = mf.template square<PTAG>(result[q]);
            }
        }
    }

    HURCHALLA_REQUEST_UNROLL_LOOP_2KARY_ARRAY_SIZE for (size_t q=0; q<ARRAY_SIZE; ++q)
        result[q] = mf.template multiply<PTAG>(result[q], val1[q]);
    return result;
