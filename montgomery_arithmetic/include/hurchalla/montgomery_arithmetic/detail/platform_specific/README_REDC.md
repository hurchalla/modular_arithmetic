# Montgomery REDC using the positive inverse (mod R)

Traditionally the Montgomery REDC algorithm uses the negative inverse (mod R) of the modulus - see Peter Montgomery's famous paper, ["Modular multiplication without trial division"](https://www.ams.org/journals/mcom/1985-44-170/S0025-5718-1985-0777282-X/home.html), which has been cited over 3000 times.  However, there exists an alternate form of this algorithm with significant performance advantages, and no apparent disadvantages.  It's remarkable that nearly all the research community has overlooked this obvious (once you see it) improvement.

I've had limited luck tracing its origins - I first encountered it in Ernst W Mayer's paper ["Efficient long division via Montgomery multiply"](https://arxiv.org/abs/1303.0328), where he presents it in introductory material without remark.  Via email, Ernst gave credit to Peter Montgomery for teaching him Montgomery multiplication using this approach in the 1990s.  After I scanned through all of Montgomery's publications, I found that Montgomery used this approach in one paper he co-authored, ["Montgomery Multiplication Using Vector Instructions"](https://www.researchgate.net/publication/274651516_Montgomery_Multiplication_Using_Vector_Instructions).  The method is treated as an incidental minor detail without citation; Mayer's paper predates this paper.  Otherwise Montgomery used the original 1985 method in all his papers and in his 1992 dissertation.  Aside from one web site's description of [Montgomery multiplication](https://cp-algorithms.com/algebra/montgomery_multiplication.html) (it's a translation of a Russian source that I haven't tracked down) and the aforementioned papers by Mayer and by Montgomery, every book, paper, and website I have been able to find uses the 1985 approach.  This includes the books [*Handbook of Applied Cryptography*](http://cacr.uwaterloo.ca/hac/), [*Prime Numbers: A Computational Perspective*](https://www.springer.com/gp/book/9780387252827), [*Modern Computer Arithmetic*](https://www.cambridge.org/gb/academic/subjects/mathematics/computational-science/modern-computer-arithmetic?format=HB&isbn=9780521194693), [Wikipedia](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication), and 20+ papers on Montgomery multiplication/reduction (a partial list: [1](https://ieeexplore.ieee.org/document/502403) [2](https://web.archive.org/web/20190327184738/http://www.hackersdelight.org/MontgomeryMultiplication.pdf) [3](https://www.researchgate.net/publication/2437595_Montgomery_Exponentiation_with_no_Final_Subtractions_Improved_Results) [4](https://eprint.iacr.org/2017/1115) [5](https://link.springer.com/article/10.1007/s13389-012-0031-5) [6](https://www.semanticscholar.org/paper/Public-key-Cryptography-on-SIMD-Mobile-Devices-Martins/9ea026c6b9baaabcf7889541c268187f2f1676ea) [7](https://www.semanticscholar.org/paper/Montgomery-Multiplication-on-the-Cell-Bos-Kaihara/91570a8ef6afebfcd8f4dcd64acc36e728e90675) [8](https://www.semanticscholar.org/paper/High-Performance-Modular-Multiplication-on-the-Cell-Bos/633d2dd5f85b1b06f1a5579a120645d55a770c37) [9](https://www.semanticscholar.org/paper/An-RNS-Montgomery-Modular-Multiplication-Algorithm-Bajard-Didier/ee9189804f8caf9b4122290348a2e0d4d2834234) [10](https://www.semanticscholar.org/paper/Systolic-Arrays-for-Modular-Exponentiation-Using-Iwamura-Matsumoto/3eb225b269cdcee1e5e14bdc5e2766f421b1fd11) [11](https://www.semanticscholar.org/paper/Montgomery-Arithmetic-from-a-Software-Perspective-Bos-Montgomery/2cc0b6b03f2186e29a471e52f9751be2b9841961) [12](https://www.semanticscholar.org/paper/On-Software-Parallel-Implementation-of-Pairings-Grabher-Gro%C3%9Fsch%C3%A4dl/f7ff0236036870af0b372961b9b8ab0efffd63d4) [13](https://link.springer.com/chapter/10.1007/11545262_6) [14](https://www.semanticscholar.org/paper/Architectural-Enhancements-for-Montgomery-on-RISC-Gro%C3%9Fsch%C3%A4dl-Kamendje/f0b34b1d821bce2490d8809ffa5c7399120dae53) [15](https://link.springer.com/chapter/10.1007/3-540-48059-5_9) [16](https://www.researchgate.net/publication/220949520_Hardware_Implementation_of_a_Montgomery_Modular_Multiplier_in_a_Systolic_Array) [17](https://www.researchgate.net/publication/224133204_Faster_Interleaved_Modular_Multiplication_Based_on_Barrett_and_Montgomery_Reduction_Methods) [18](https://www.semanticscholar.org/paper/New-Speed-Records-for-Montgomery-Modular-on-8-Bit-Liu-Gro%C3%9Fsch%C3%A4dl/a5c643394e2fc2c0aa7e940c4bf7cae16d26e068) [19](https://www.semanticscholar.org/paper/Montgomery-Modular-Multiplication-on-ARM-NEON-Seo-Liu/44604c49fc8642d99c44c83e54f4ac4e411cd4c9) [20](https://www.semanticscholar.org/paper/Montgomery-Modular-Multiplication-Algorithm-on-Fan-Sakiyama/29a1eaa8d9ac3e04940831bfdead45067c048931) [21](https://www.semanticscholar.org/paper/New-frameworks-for-Montgomery's-modular-method-McLaughlin/e25cca2c65c69aed668f30c08792519e8fba121e) [22](https://www.semanticscholar.org/paper/Defeating-modexp-side-channel-attacks-with-traces-Granlund/3139b32f21fafb334498cff777ad11d77596db35) [23](https://www.semanticscholar.org/paper/Speeding-the-Pollard-and-elliptic-curve-methods-of-Montgomery/a377c34ff3252c94ae864b7e0cfa64c05d01ef6d) [24](https://www.semanticscholar.org/paper/Montgomery-reduction-within-the-context-of-residue-Bajard-Eynard/05e309ec18fd150df416f9629597e68f0e78459c))

#### Background

It's perhaps plausible that Peter Montgomery discovered the 1985 algorithm after asking himself the question, "could I add some multiple of N to an integer T, to make the sum divisible by R (assuming R is a power of 2 and N is odd)?"  The answer is yes, and his 1985 paper shows how to calculate that multiple m and use it to in turn calculate TR<sup>-1</sup> (mod N), which is the basis for Montgomery modular multiplication.  [Wikipedia](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm) gives a fairly easy and thorough description of the algorithm and proof.

#### The Alternate Version

We'll get a slightly different algorithm if we ask a slightly different question - "could we subtract some multiple of N from an integer T, to make the difference divisible by R?"  We can show that the answer is yes.  Some integer 'm' always exists such that <span class="nowrap">T - mN &equiv; 0 (mod R)</span>.  And we can provide a way to calculate m, and <span class="nowrap">TR<sup>-1</sup> (mod N)</span>.

Given that R is a power of 2 and N is odd, we know the greatest common divisor of R and N is 1, and thus by [B&#233;zout's identity](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_identity) there exist integers x and y such that <span class="nowrap">xN + yR = 1</span>.  Taking this mod R, we have <span class="nowrap">xN &equiv; 1 (mod R)</span>, and so x is the [modular multiplicative inverse](https://en.wikipedia.org/wiki/Modular_multiplicative_inverse) of N (mod R), which we denote as <span class="nowrap">x &equiv; N<sup>-1</sup> (mod R)</span>.  Multiplying both sides by T, we get  <span class="nowrap">xTN &equiv; T (mod R)</span>, which after rearranging is  <span class="nowrap">T - xTN &equiv; 0 (mod R)</span>.  Therefore our desired integer 'm' equals xT.  More precisely, <span class="nowrap">m &equiv; N<sup>-1</sup>T (mod R)</span>, and <span class="nowrap">T - mN &equiv; 0 (mod R)</span>.

A different way to get the same result is to follow a line of reasoning similar to that used in the [Wikipedia REDC article](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm).  Let <span class="nowrap">m = ((T mod R)N<sup>-1</sup>) mod R</span>.  Then <span class="nowrap">m &equiv; N<sup>-1</sup>T (mod R)</span> and </br>
<span class="nowrap">T - mN &equiv; T - (TN<sup>-1</sup>)N &equiv; T - T &equiv; 0 (mod R).</span>

Continuing to follow reasoning similar to the Wikipedia article, we can see that the above result proves that <span class="nowrap">T - mN</span> is divisible by R, and therefore <span class="nowrap">t = (T - mN)/R</span>  must be an integer.  This means that</br>
<span class="nowrap">t &equiv; (T - mN)R<sup>-1</sup> &equiv; TR<sup>-1</sup> - (mR<sup>-1</sup>)N &equiv; TR<sup>-1</sup> (mod N).</span>

So we have a way to calculate both m and <span class="nowrap">TR<sup>-1</sup> (mod N)</span>.

Finally, let's deduce bounds for t, assuming we have bounds on the inputs of <span class="nowrap">0 < N < R</span>, and <span class="nowrap">0 <= T < RN</span>.  Since m is calculated mod R, we know <span class="nowrap">0 <= m < R</span>, and so <span class="nowrap">-R < -m <= 0</span> and <span class="nowrap">-RN < -mN <= 0</span>.  Thus, <span class="nowrap">-RN <= T - RN < T - mN <= T < RN</span>.  This gives us  <span class="nowrap">-RN < T - mN < RN</span>, and since all parts are divisible by R, we have  <span class="nowrap">-N < (T - mN)/R < N</span>.  And thus we have the bounds for t:</br>
<span class="nowrap">-N < t < N</span>.

We can now write our alternate REDC algorithm:

<pre><b>function</b> REDC2 <b>is</b>
    <b>input:</b> Integers <i>R</i> and <i>N</i> with <span class="nowrap">gcd(<i>R</i>, <i>N</i>) = 1</span>, and <i>N</i> in <span class="nowrap">[0, <i>R</i> &#8722; 1]</span>,
           Integer <i>N</i><sup>&#8722;1</sup> in <span class="nowrap">[0, <i>R</i> &#8722; 1]</span> such that <span class="nowrap"><i>NN</i><sup>&#8722;1</sup> &equiv; 1 (mod <i>R</i>)</span>,
           Integer <i>T</i> in the range <span class="nowrap">[0, <i>RN</i> &#8722; 1]</span>
    <b>output:</b> Integer <i>S</i> in the range <span class="nowrap">[0, <i>N</i> &#8722; 1]</span> such that <span class="nowrap"><i>S</i> &equiv; <i>TR</i><sup>&#8722;1</sup> (mod <i>N</i>)</span>

    <i>m</i> &lArr; ((<i>T</i> mod <i>R</i>)<i>N</i><sup>&#8722;1</sup>) mod <i>R</i>
    <i>t</i> &lArr; (<i>T</i> &#8722; <i>mN</i>) / <i>R</i>
    <b>if</b> <i>t</i> &lt; 0 <b>then</b>
        <b>return</b> <span class="nowrap"><i>t</i> + <i>N</i></span>
    <b>else</b>
        <b>return</b> <i>t</i>
    <b>end if</b>
<b>end function</b>
</pre>

</br>

#### Why It Matters
On the surface, the alternate REDC algorithm and the traditional 1985 version look very similar.  However, they differ in important ways when implemented, with the result that the alternate version is significantly simpler and more efficient.  Primarily the implementation advantage arises from the alternate version's use of <span class="nowrap">(<i>T</i> &#8722; <i>mN</i>)</span> verses the traditional version's use of <span class="nowrap">(<i>T</i> + <i>mN</i>).</span>  &nbsp;&nbsp;<i>mN</i> is calculated differently in the two versions: for the alternate version, <span class="nowrap"><i>T</i> &#8722; <i>mN</i> &equiv; 0 (mod <i>R</i>),</span> and for the traditional version, <span class="nowrap"><i>T</i> + <i>mN</i> &equiv; 0 (mod <i>R</i>).</span>  So for the former implementation, the low half bits of the subtraction will always equal 0 and will never generate a borrow/carry.  For the latter implementation, the low half bits of the addition will always equal 0, and will always generate a carry unless both <i>T</i> and <i>mN</i> equal 0.  This exception for <i>T</i> and <i>mN</i> equaling 0 means that the carry must always be calculated for the traditional version (usually by performing the low part addition), and then an add-with-carry instruction (or its equivalent) must be used to sum the high part of the addition.  In contrast, the alternate algorithm's low part subtraction <span class="nowrap"><i>T</i> &#8722; <i>mN</i></span> can be completely ignored and omitted from its implementation, because we can prove it never has any effect on the rest of the algorithm - the result of the low part subtraction is always zero and it never generates a borrow/carry.  A further difference between the two algorithm versions' implementations is that the traditional version usually has two conditionals that must be checked: possible overflow on <span class="nowrap"><i>t &lArr; T</i> + <i>mN</i>,</span> and <span class="nowrap"><i>t</i> &ge; <i>N</i>;</span>  in contrast, the alternate version has only one conditional check: <span class="nowrap"><i>t</i> &lt; 0,</span> which is itself an overflow check.

It may be easiest to see the differences by comparing assembly code.  We'll look at two x86_64 inline assembly C functions, implementing the traditional REDC and the alternate REDC for a 64 bit unsigned integer type.</br></br>

#### x86_64 Assembly Implementations

We'll first consider the traditional REDC x86_64 inline asm implementation:

<pre>
// ~11 cycles latency, 11 fused uops.
inline uint64_t REDC_traditional(uint64_t T_hi, uint64_t T_lo, uint64_t N,
                                                               uint64_t negInvN)
{
    assert(T_hi < N);   // REDC requires T < NR, and this enforces it.
    uint64_t rrax = T_lo;
    uint64_t rrdx, tmp;
    __asm__ (
        "movq %%rax, %[tmp] \n\t"
        "imulq %[inv], %%rax \n\t"    /* m = T_lo * negInvN */
        "mulq %[N] \n\t"              /* mN = m * N */
        "addq %[tmp], %%rax \n\t"     /* rax = T_lo+mN_lo. Sets carry, rax==0 */
        "adcq %[Thi], %%rdx \n\t"     /* t_hi = addcarry(T_hi, mN_hi) */
        "cmovaeq %[N], %%rax \n\t"    /* rax = (t_hi >= T_hi) ? N : 0 */
        "xorl %k[tmp], %k[tmp] \n\t"  /* tmp = 0 */
        "subq %[N], %%rdx \n\t"       /* rdx = t_hi - N */
        "cmovaeq %[tmp], %%rax \n\t"  /* rax = (t_hi >= N) ? 0 : rax */
        : "+&a"(rrax), "=&d"(rrdx), [tmp]"=&r"(tmp)
        : [Thi]"r"(T_hi), [N]"r"(N), [inv]"r"(negInvN)
        : "cc");
    return rrax + rrdx;   // let compiler choose between add/lea
}
</pre>
<i>Traditional REDC</i>

</br>

The code above is decent but not completely ideal.  We can see that it incorporates an add-carry and two conditional moves, as we discussed.  We will see below that the alternate REDC does not need an add-carry (or a subtract-borrow), and it also needs only one conditional move instead of two.

Another argument against this code (as given) is simply that it is inline asm.  Our absolute ideal would be code written in only standard C (using no inline asm at all) that compiles to assembly that is nearly as optimal as our inline asm functions.  This is easy to wish for, though not always easy to achieve.  If we write a standard C version of the traditional REDC, we can see from our inline asm's use of conditional move and add-carry instructions that we may have trouble coaxing the C compiler into generating optimal assembly; those instructions don't exist in C.  We could use C compiler idioms to help the compiler, such as the ternary operator (which typically translates into a conditional move), but we still will find none of the major compilers (gcc, clang, MSVC, icc) produce optimal assembly from a standard C version of the traditional REDC.  This isn't necessarily a major problem, since we can of course use the given inline asm, but we should note that some compilers don't support inline asm (MSVC 64bit) and also that inline asm code tends to be unusually bug prone.  I can give you comfortable assurances that I'm skilled at writing high quality inline asm, but you would be wise to be skeptical of this and of any inline asm you see - including our inline asm functions here.  Unit tests of inline asm are far less helpful than you might think, potentially producing false negatives caused by the particular register choices your compiler happens to make for code surrounding your inline asm under test.  There are reasonable arguments to be made for [banning inline asm](https://gcc.gnu.org/wiki/DontUseInlineAsm) from your projects.

We'll briefly note that we could rewrite the inline asm above to [improve the traditional REDC](https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/experimental/README_REDC_supplement.md).  However, in nearly all respects the alternate REDC will provide an even greater improvement, so this is mentioned for curiosity sake, or for the case that you happen to want a drop-in replacement for an existing traditional REDC function.

Below we have the alternate REDC algorithm inline asm C function, which comes close to our ideal.  Note that it uses the positive inverse as a parameter, whereas the traditional algorithm uses the negative inverse as a parameter.  The function consists of a relatively straightforward translation of the alternate algorithm into inline asm, and thus it's fairly easy to understand.  It can be expressed as standard C code almost as effectively as with inline asm; clang produces close to optimal assembly from a standard C version of it (other major compilers produce assembly that is acceptable but not as good).  For x86_64, it provides us with the best performance we've seen, with both lowest latency and fewest used instructions/uops/registers.  I would expect the alternate REDC to provide similar benefits for ARM and other architectures.</br></br>

<pre>
// 9 cycles latency, 7 fused uops.
inline uint64_t REDC_alternate(uint64_t T_hi, uint64_t T_lo, uint64_t N,
                                                                  uint64_t invN)
{
    assert(T_hi < N);   // REDC requires T < NR, and this enforces it.
    uint64_t rrax = T_lo;
    uint64_t Thi = T_hi;
    __asm__ (
        "imulq %[inv], %%rax \n\t"        /* m = T_lo * invN */
        "mulq %[N] \n\t"                  /* mN = m * N */
        "leaq (%[Thi], %[N]), %%rax \n\t" /* rax = T_hi + N */
        "subq %%rdx, %%rax \n\t"          /* rax = rax - mN_hi */
        "subq %%rdx, %[Thi] \n\t"         /* t_hi = T_hi - mN_hi */
        "cmovbq %%rax, %[Thi] \n\t"       /* t_hi = (T_hi&lt;mN_hi) ? rax : t_hi */
        : [Thi]"+&bcSD"(Thi), "+&a"(rrax)
        : [N]"r"(N), [inv]"r"(invN)
        : "rdx", "cc");
    return Thi;
}
</pre>
<i>Alternate REDC (Positive Inverse)</i>

</br>

All code in this file is licensed under the MIT Open Source License:

Copyright (c) 2020 by Jeffrey Hurchalla.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
