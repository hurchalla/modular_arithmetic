This file supplements the document [README_REDC.md](README_REDC.md).

We can improve upon the inline assembly we saw for traditional REDC, though the code becomes harder to understand.  The improvements also can't be implemented well in standard C; none of the major compilers (gcc, clang, MSVC, icc) are able to compile standard C versions of these functions without adding significant extra latency and uops, even with idiomatic use of the ternary operator for conditional move.

Since the alternate REDC function from [README_REDC.md](README_REDC.md) does better on uops and equals or betters the latency, all while being easier to understand, and friendlier for compilers if written in standard C, we should certainly prefer the alternate REDC to the functions that follow.  Nevertheless the functions below do improve the traditional REDC inline asm, so they could be useful as an easy drop-in replacement of an existing REDC function (which will almost certainly be traditional REDC with the negative inverse), or they might be interesting for anyone curious.

The improved functions below are correct and produce output equivalent to the previous inline asm we saw for the traditional REDC.  You can find a rough proof of correctness in comments of the C++ function ["REDC(T u_hi, T u_lo, T n, T neg_inv_n, FullrangeTag, InplaceLowlatencyTag)" of an old git commit](https://github.com/hurchalla/modular_arithmetic/blob/66281af1639031b04bdaf9b916e5d5638d3ded25/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/platform_specific/RedcLargeR.h#L365).

The first function below is optimized for lowest latency(REDC_traditional_improved1).  The second is optimized for lowest uops (REDC_traditional_improved2):</br>


<pre>
// On Intel Skylake: 9 cycles latency, 11 fused uops
inline uint64_t REDC_traditional_improved1(uint64_t T_hi, uint64_t T_lo,
                                                   uint64_t N, uint64_t negInvN)
{
    assert(T_hi < N);   // REDC requires T < NR, and this enforces it.
    uint64_t rrax = T_lo;
    uint64_t Thi = T_hi;
    uint64_t tmp;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"
        "imulq %[inv], %%rax \n\t"    /* m = T_lo * negInvN */
        "mulq %[N] \n\t"              /* mN = m * N */
        "movq %[Thi], %%rax \n\t"
        "subq %[N], %%rax \n\t"       /* diff = T_hi - N */
        "negq %[tmp] \n\t"            /* Sets carry to (T_lo != 0) */
        "adcq %%rdx, %[Thi] \n\t"     /* sum1 = addcarry(T_hi, mN_hi) */
        "negq %[tmp] \n\t"            /* Sets carry to (T_lo != 0) */
        "adcq %%rdx, %%rax \n\t"      /* sum2 = addcarry(diff, mN_hi) */
        "cmovaeq %[Thi], %%rax \n\t"  /* rrax = (sum2 >= mN_hi) ? sum1 : sum2 */
        : [Thi]"+r"(Thi), "+&a"(rrax), [tmp]"=&r"(tmp)
        : [N]"r"(N), [inv]"r"(negInvN)
        : "rdx", "cc");
    return rrax;
}
</pre>
<i>Improved Traditional REDC (low latency version)</i>

</br>

<pre>
// On Intel Skylake: 10 cycles latency, 9 fused uops.
inline uint64_t REDC_traditional_improved2(uint64_t T_hi, uint64_t T_lo,
                                                   uint64_t N, uint64_t negInvN)
{
    assert(T_hi < N);   // REDC requires T < NR, and this enforces it.
    uint64_t rrax = T_lo;
    uint64_t Thi = T_hi;
    uint64_t tmp;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"
        "imulq %[inv], %%rax \n\t"        /* m = T_lo * negInvN */
        "mulq %[N] \n\t"                  /* mN = m * N */
        "subq %[N], %[Thi] \n\t"          /* diff = T_hi - N */
        "negq %[tmp] \n\t"                /* Sets carry to (T_lo != 0) */
        "adcq %[Thi], %%rdx \n\t"         /* rdx = addcarry(diff, mN_hi) */
        "leaq (%%rdx, %[N]), %%rax \n\t"  /* rax = rdx + N */
        "cmovbq %%rdx, %%rax \n\t"        /* rrax = (rdx &lt; mN_hi) ? rdx : rax */
        : [Thi]"+&r"(Thi), "+&a"(rrax), [tmp]"=&r"(tmp)
        : [N]"r"(N), [inv]"r"(negInvN)
        : "rdx", "cc");
    return rrax;
}
</pre>
<i>Improved Traditional REDC (low uops version)</i>

</br>

All code in this file is licensed under the MIT Open Source License:

Copyright (c) 2022 by Jeffrey Hurchalla.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
