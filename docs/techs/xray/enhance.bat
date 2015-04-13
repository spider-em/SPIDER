;
; Calculate 1D pw, 
; then create scattering curve in the same format,
; search for close frequences -- intensity can be normalized/ multiplied
; or divided to get it closer to pw values

FR G
<vol>bpr004    ; input volume to correct/enhance

FR G
<scatter>scattering282     ; input X-ray scattering power spectrum

FR G
<pow>rojo     ; output 1D power spectrum of input volume

FR G
<out>fen      ; output enhancement curve


x80 = 2.82    ; pixel size (in A)

x35 = 35     ; filter limit (in pixels). Filter the data out to point P,
             ; where P < (volume_diameter / 2)
             ; E.g., resolution cutoff

; -------------------------------------------------
; create a doc file with a 1D power spectrum of the input volume
PW
<vol>
_1

SQ
_1
_2

RO
_2
_3

DE
<pow>

LI D
_3
<pow>
r
1

UD N,x66
<pow>

; create the output enhancement curve
x55=0.0
x56=100
x21=1

DE
<out>

DO LB1 x21=2,x35   ; curve goes out to filter limit

UD IC,x21,x51
<pow>

x55=(x21-1)/(2*(x66-1))
x56=x80/x55

@pwsc[x56,x77]
<scatter>

; Xray/EM
x77=SQR(x77/X51)
x78=LON(X77)

SD x21,x77,x55,x56,X78
<out>

LB1

EN D

