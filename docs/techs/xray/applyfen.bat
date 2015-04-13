; Apply the enhancement/correction to the volume,
; Then filter the volume at the resolution level.


x50 = 0.32      ; resolution level (low pass filter cutoff)

FR G
<vol>bpr004     ; input volume to be enhanced

FR G
<fen>fen      ; input enhancement curve

FR G
<envol>vol_en       ; output corrected volume

FR G
<fvol>vol_en_filt  ; output corrected, filtered volume

; ------------------------------------------


; Filter using the enhancement file
FD      
<vol>
<envol>
<fen>

FQ 
<envol>
<fvol>
7
x50-0.05,x50+0.05  

EN D
