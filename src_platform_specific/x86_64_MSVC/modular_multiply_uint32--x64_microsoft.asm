
.code

; This uses Microsoft calling convention


; extern "C" uint32_t modular_multiply_uint32_asm_UID7b5f83fc983(uint32_t a,
;                                              uint32_t b, uint32_t modulus);
; Preconditions: 0 <= a < modulus,  0 <= b < modulus,  modulus > 0
; Postconditions: returns (a*b)%modulus
;
; rcx == a, rdx == b, r8 == modulus
; return register is rax
PUBLIC  modular_multiply_uint32_asm_UID7b5f83fc983
modular_multiply_uint32_asm_UID7b5f83fc983  PROC
    mov eax, ecx
    mul edx         ; EDX:EAX = EAX*EDX; high-order bits of the product in EDX
    div r8d         ; (quotient EAX, remainder EDX) = EDX:EAX/R8D
    mov eax, edx    ; return the remainder
    ret 0
modular_multiply_uint32_asm_UID7b5f83fc983  ENDP


End
