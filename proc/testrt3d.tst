; PURPOSE TEST  'RT 3D'  Bimal Mar 2000
;               'MO 3'   Helix  al Aug 2002

 MO 3                    ; Create model volume
   _1                    ; Inline volume  (output)
   100  100 100          ; Dimensions
   W                     ; Density wedge

 RT 3D                   ; Rotate the volume
   _1                    ; Inline volume  (input)
   _2                    ; Rotated file  (output)
   10.0,20.0,30          ; Angles
   30.0 

 ; Find statistics
 FS [max],[min],[avg]    ; Get volume statistics
   _2                    ; File       (input)

 SYS
   echo ' Range:    ({%f8.2%[min]} ... {%f8.2%[max]})'
 SYS
   echo ' Shud be:  (    0.02 ...    35.00)'

 ; Average pixel value should be: 12.107
 [avgt] = INT([avg]*1000)
 [avgt]

 IQ REG
   [avgt] 12107

 MO 3                 ; Test 'MO 3' helix
   jnk000             ; Output volume      (output)
   100 100 100        ; Volume dimensions
   H                  ; Helix
   2                  ; Sphere density
   3 40               ; Sphere radius, helix radius
   60  3              ; 60 spheres in 3 turns

 RE 
