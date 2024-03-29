; Copyright (c) 2020-2022 Jeffrey Hurchalla.
; This Source Code Form is subject to the terms of the Mozilla Public
; License, v. 2.0. If a copy of the MPL was not distributed with this
; file, You can obtain one at https://mozilla.org/MPL/2.0/.

.code

; This uses Microsoft x64 calling convention


; extern "C" uint64_t modular_multiply_uint64_asm_UID7b5f83fc983(uint64_t a,
;                                              uint64_t b, uint64_t modulus);
; Preconditions: 0 <= a < modulus,  0 <= b < modulus,  modulus > 0
; Postconditions: returns (a*b)%modulus
;
; rcx == a, rdx == b, r8 == modulus
; return register is rax
PUBLIC  modular_multiply_uint64_asm_UID7b5f83fc983
modular_multiply_uint64_asm_UID7b5f83fc983  PROC
    mov rax, rcx
    mul rdx         ; RDX:RAX = RAX*RDX; high-order bits of the product in RDX
    div r8          ; (quotient RAX, remainder RDX) = RDX:RAX/R8
    mov rax, rdx    ; return the remainder
    ret 0
modular_multiply_uint64_asm_UID7b5f83fc983  ENDP


End
