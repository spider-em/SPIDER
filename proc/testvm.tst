; PURPOSE:  Loop for renaming/renumbering series of SPIDER images
; ArDean Leith Sept. 99  (INDEXED)

MO                  ; Create test file
jnk001              ; Image          (output)
63 56               ; Size
T                   ; Sine waves

MO                  ; Create test file
jnk002              ; Image          (output)
63 56               ; Size
T                   ; Sine waves

MO                  ; Create test file
jnk004              ; Image          (output)
63 56               ; Size
T                   ; Sine waves

; ----------------------------------------------------------------------

; Loop for renaming/renumbering series of SPIDER images 

[num]=0              ; Initialize output file number counter

DO [i]=1,2           ; Loop

   [num]=[num]+1     ; Increment output file number 

   VM                ; System call for renaming/renumbering 
   cp jnk{***[i]}.$DATEXT jnkout{***[num]}.$DATEXT 

ENDDO                ; End loop

VM
ls -l jnkout*.$DATEXT

; ----------------------------------------------------------------------

; Loop for consecutively renaming/renumbering SPIDER image series 
;    ( skipping any nonexistent image numbers) 

[num]=0                 ; Initialize output file number
DO  [n]=1,4             ; Loop

   IQ FI [got]          ; Use "IQ FI" to see if file exists
   jnk{***[n]}          ; Filename                    (input)

   IF([got].LE.0)CYCLE  ; Skip if file not found

   [num]=[num]+1        ; Increment output file number

   VM                   ; System call for renaming/renumbering 
   cp jnk{***[n]}.$DATEXT jnkout{***[num]}.$DATEXT 
ENDDO                   ; End loop

RE
