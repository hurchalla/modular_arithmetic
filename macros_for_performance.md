
Optional macros you can define to tune performance
--------------------------------------------------
There are a number of macros you can optionally define when you compile, to tune the
performance on your system for the modular arithmetic functions.  It is
generally recommended not to do so, but in some cases you may find it useful.
You would define one or more of these macros when compiling *your* sources,
given that Clockwork is a header-only library.

For example, if you are compiling using clang or gcc from the command line, you would
specify the -D compilation flag, similarly to the following: 
clang++ -DHURCHALLA_ALLOW_INLINE_ASM_ALL  ...more arguments...
As another example, if you are using CMake you would add the command "target_compile_definitions"
to your CMakeLists.txt, similarly to the following: 
target_compile_definitions(&lt;your_target_name&gt;  PRIVATE  HURCHALLA_ALLOW_INLINE_ASM_ALL) 
\
\
HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE - define this macro if your target system's
instruction set does not include division.  Although it is unusual, some
microcontrollers do not have division, and defining this macro might improve
performance in such a case.

HURCHALLA_AVOID_CSELECT - you may wish to define this macro if your target
system's instruction set does not include conditional move or conditional
select.  It may improve performance in such a case.  This macro is normally
already defined for RISC-V.

HURCHALLA_ALLOW_INLINE_ASM_ALL - defining this macro will enable all
available inline asm functions.  Although this is the easiest macro to use, you
can more selectively enable inline asm for particular functions, using macros
listed below.  In some cases HURCHALLA_ALLOW_INLINE_ASM_ALL may improve
performance up to 20% (gcc often benefits), and in other cases it may make
essentially no difference or harm performance (clang does not seem to benefit).
It is not enabled by default because inline asm is extremely difficult to verify
for correctness.  While I believe I'm skilled at writing high quality inline
asm, I advise you to be skeptical of this and of any inline asm you see.
Unit tests of inline asm are far less helpful than you might think - the ability
of a unit test to detect a bug in inline asm often depends upon the register
allocation choices the compiler makes for surrounding test code, which is mostly
outside a programmer's control.  Generally speaking, it is [difficult to
recommend inline asm](https://gcc.gnu.org/wiki/DontUseInlineAsm) unless there is
a large performance benefit or performance is critical.

HURCHALLA_ALLOW_INLINE_ASM_REDC  
HURCHALLA_ALLOW_INLINE_ASM_ABSDIFF  
HURCHALLA_ALLOW_INLINE_ASM_MODADD  
HURCHALLA_ALLOW_INLINE_ASM_MODSUB  
HURCHALLA_ALLOW_INLINE_ASM_QUARTERRANGE_GET_CANONICAL  
HURCHALLA_ALLOW_INLINE_ASM_HALFRANGE_GET_CANONICAL  
- these macros selectively enable inline asm for functions.  They may or may not
improve performance, and the warnings above for HURCHALLA_ALLOW_INLINE_ASM_ALL
apply here too.  To determine if they are even useful, you would need to
compare performance with different ASM macros defined/not defined.  Generally
you would want to start with HURCHALLA_ALLOW_INLINE_ASM_REDC.
