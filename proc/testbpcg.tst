; SOURCE:   testbpcg.tst  INDEXED
; PURPOSE:  Test 'BP CG' 

MO                ; Create sine wave image 
jnk001            ; Image file    (output)
75 75             ; Size
T                 ; Sine wave

DO [inum]=2,12 ; Copy sine wave to 12 images
   CP             ; SPIDER Copy
   jnk001         ; File           (input)
   jnk{***[inum]} ; File           (output)
ENDDO

VM
echo ' Testing BP CG with 12 images -----'

BP CG          ; Run  backprojection
jnk***         ; Template for 2-d images  (input)
1-12           ; File numbers
35             ; Radius of object
savangbprp     ; Angles doc file          (input)
N              ; No symmmetries
jnkbpcg_out    ; Volume                   (output)
1.0E5, 0.0     ; Error limit, Chi^2 limit
25,1           ; Iteration limit, mode
2000           ; Lambda

FS [max],[min] ; Get volume statistics
jnkbpcg_out    ; File                        (input)

VM             ; Echo volume statistics
echo ' BP CG    Range: {%g13.5%[min]}...{%g13.5%[max]}'
VM
echo ' Correct  Range:  -0.10567E-01...  0.95824E-02'
VM            
echo '  '

RE


