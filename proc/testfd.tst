;
; PURPOSE: Test FD (filter using doc file) ArDean Leith May 00
; (INDEXED)

FD                ; Filter according to a doc. file
sav0020           ; Image                         (input)
jnk321            ; Filtered image                (output)
savdoc006         ; Doc file containing coeffs.   (input)

FS [max],[min]    ; Get volume statistics
jnk321            ; Filtered image                (input)

VM
echo ' Range:     ({%f8.4%[min]} ... {%f8.4%[max]})'
VM
echo ' Should be: ( -0.0350 ...  -0.0038) (for 80 x 1000)'

; FMIN=-0.350E-01 FMAX=-0.379E-02 AV=-0.16118E-01  SIG= 0.86503E-02
; FMIN=-0.238E-01 FMAX= 0.245E-02 AV=-0.80590E-02  SIG= 0.70794E-02
; last line is june 2000 results

RE
