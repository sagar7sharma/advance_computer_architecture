;Group 11 Members: Sagar Sharma (368058), Teodora Anitoaei (368014), Varun Gowtham  (368053)
.data
N_COEFFS: .word 3
coeff: .double 0.5,1.0,0.5
N_SAMPLES: .word 5
sample: .double 1.0,2.0,1.0,2.0,1.0
result: .double 0.0,0.0,0.0,0.0,0.0
base: .double 1.0,0.0
.text
Main:
    l.d f4,coeff(r0)        ;loading coeff(0)
    daddi r4,r0,8           ;increment reference reg, point to next word
    l.d f0,base(r4)         ;load f0 (zero comparison reg) with zero, to be sure 
    l.d f1,base(r4)         ;load f1 (norm reg) with zero 
    l.d f5,coeff(r4)        ;loading coeff(1)
    daddi r5,r4,8           ;increment reference reg, to point to next word
    l.d f6,coeff(r5)        ;loading coeff(2)
    c.lt.d f4,f0            ;checking if coeff(0) is less than 0
    bc1t L1                 ;branching to subtract
    nop                     ;
    add.d f1,f1,f4          ;add coeff(0) to norm reg
    j G1                    ;
    l.d f7,base(r0)         ;load numerator reg with 1 to compute 1/norm
L1: sub.d f1,f1,f4          ;make coeff(0) +ve and add to norm
    l.d f7,base(r0)         ;
G1: c.lt.d f5,f0            ;checking if coeff(1) is less than 0
    bc1t L2                 ;
    nop                     ;
    add.d f1,f1,f5          ;add coeff(1) to norm reg
    j G2                    ;
    daddi r4,r0,8           ;make r4 as 8 to compute words later
L2: sub.d f1,f1,f5          ;make coeff(1) +ve and add to norm
    daddi r4,r0,8           ;
G2: c.lt.d f6,f0            ;checking if coeff(2) is less than 0
    bc1t L3                 ;
    nop                     ;
    add.d f1,f1,f6          ;add coeff(2) to norm reg
    j OUT                   ;
    ld r2,N_SAMPLES(r0)     ;
L3: sub.d f1,f1,f6          ;make coeff(2) +ve and add to norm
    ld r2,N_SAMPLES(r0)     ;r2 contains the number of samples
OUT: div.d f1,f7,f1         ;compute 1/norm and store it in norm reg

     daddi r7,r2,-1         ;r7 contains make N_SAMPLES-1
     l.d f2,sample(r0)      ;load the first sample to f2
     dmulu r9,r7,r4         ;increment r9 to point to last sample
     s.d f2,result(r0)      ;result[0] = sample[0]
     daddi r7,r7,-1         ;r7 pointing to N_SAMPLES-2
     l.d f2,sample(r9)      ;To make
     s.d f2,result(r9)      ;result[N_SAMPLES-1] = sample[N_SAMPLES-1]
     mul.d f4,f4,f1         ;compute coeff(0)/norm
     mul.d f5,f5,f1         ;compute coeff(1)/norm
     mul.d f6,f6,f1         ;compute coeff(2)/norm

FOR: l.d f13,sample(r4)     ;load f13 with sample(n)
     mul.d f13,f13,f5       ;f13=sample(n)*(coeff(1)/norm)
     daddi  r14,r4,-8       ;r14 points to sample(n-1)
     l.d f14,sample(r14)    ;load f14 with sample(n-1)
     mul.d f14,f14,f4       ;f14=sample(n-1)*(coeff(0)/norm)
     daddi  r15,r4,8        ;r15 points to sample(n+1)
     l.d f15,sample(r15)    ;load f15 with sample(n+1)
     add.d f16,f13,f14      ;f16=f13+f14
     mul.d f15,f15,f6       ;f115=sample(n+1)*(coeff(2)/norm)
     add.d f17,f16,f15      ;f17=f15+f16
     daddi r7,r7,-1         ;decrement loop reference r7
     s.d f17,result(r4)     ;store result 
     bnez r7,FOR            ;loop back to for if r7!=0
     daddi r4,r4,8          ;point r7 to next word in the sample
     halt




