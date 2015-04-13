; apply the enhancement/correction to the volume,
; Then filter the volume at the resolution level.

; Fit parameters
x31 = 164.471    ; A
x32 = 32521.7    ; B
x33 = -2.04786   ; C   - set to zero if only used A,B

x50 = 0.32      ; resolution level (low pass filter cutoff)

x35 = 35     ; filter limit (in pixels). see enhance.bat

FR G
<vol>bpr004     ; input volume to be enhanced

FR G
<filter>filter  ; output document file with filter values

FR G
<fvol>vol_filt  ; output corrected, filtered volume

; ---------------------------------------------------------

DE
filter

DO LB1 x21=2,x35
   x77 = (x33*x21*x21*x21) + (x31*x21*x21) + x32
   SD x21,x77
   <filter>
LB1

SD E
<filter>

; Filter using a document file
FD
<vol>     ; input volume
_1        ; output corrected volume
<filter>  ; document file to use as filter values

; Low pass filter to the resolution level
FQ
_1
<fvol>
7
x50-0.05,x50+0.05

EN D
