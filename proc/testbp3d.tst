;
; !!!!  FS FAILS if weighting used !!!!!!!!!!!!!!!!!

; create inline stack for 2 images
MS
_8@         ; Inline stack
75,75,1     ; size
2           ; 20 images

mo
jnk_bp3d
75 75
t

cp 
jnk_bp3d
_8@1

; Run  backprojection now
BP 3D         
_8@****       ; Template for 2-d input images 
1             ; Projection files number to be used
savangbprp    ; Angles doc file
75,75,75      ; Output volume dimensions
1,5           ; First, last slice to be reconstructed
1             ; SNR / DIAMETER
jnkbp3dout    ; Output volume

FS
_8@1
;;jnkbp3dout

RE
