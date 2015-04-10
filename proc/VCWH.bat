[x24,x25,x26,x27,x77,x78,x79,x80,x11,x12,x16,x17,x18,x87,x88,x91]

;VolCorWH--slave		
; this is the file that will be run by VolCor to find correlation
; values for 2 volumes rough draft 11/24/03 jamie lebarron
; accepted large volume values from VolCorr for explorer, 12/18/03
; added size of contol points box and changed ouput to vector 12/23/03

;---ONE TIME RUN CALCULATIONS-----------
;---calculate values to remove center of window after
;---cross correlation function
x46=INT(x80/2)+1		;x value =Bx33
x47=INT(x11/2)+1
x48=INT(x12/2)+1

;---line counter for file print out, must increment every new peak
x50=0

;---calculate size of half cube, for finding begining 
;---cube for loops
x81=(x80+1)/2		;x half cube
x82=(x11+1)/2		;y half cube
x83=(x12+1)/2		;z half cube, altered for numbers sake

;---error catch if cube is larger than partition
IF((x18-x83).LE.0)THEN
 VM
 echo "cube larger than Z slab"
 END 
ENDIF

IF((x17-x82).LE.0)THEN
 VM
 echo "cube larger than Y slab"
 END 
ENDIF

IF((x16-x81).LE.0)THEN
 VM
 echo "cube larger than X slab"
 END 
ENDIF
;---end error catch

;---total number of loops =X*Y*Z
x98=x77*x78*x79

;---start codon for explorer/datascribe
x69=666
;---calculate the total size of the volume
x21=(x77*x24)+(x26*(x77-1))	;x
x22=(x78*x25)+(x27*(x78-1))	;y
x23=(x79*x87)+(x88*(x79-1))	;z

;---send the size of the volume to docfile as a comment
;---it will always be preseeded by "666"
SD -1,x69,x77,x78,x79,x21,x22,x23
[doc_out]

;---clear values for their use as loop counters
x21=0
x22=0
x23=0

;---END ONE TIME RUN CALCULATIONS-----------

;---BEGIN GIGUNDO XYZ LOOPS
DO LB3 x23=1, x79		;Z loop
     VM
     echo "Entering Z loop number " {****x23}"of"{****x79}
  DO LB2 x22=1, x78		;Y loop
    DO LB1 x21=1, x77		;X loop
     
;---loop counter
x99=(x23-1)*x78*x77+(x22-1)*x77+x21
VM
echo "entering loop" {***x99} "of " {***x98}  

;---values to find upper left corner of larger cube in primary
x96=(x23-1)*x87+1
IF(x23.GT.1) x96=x96+(x88*(x23-1))	;upper left Z value

x28=(x22-1)*x25+1
IF(x22.GT.1) x28=x28+(x27*(x22-1))	;upper left Y value

x29=(x21-1)*x24+1
IF(x21.GT.1) x29=x29+(x26*(x21-1))	;upperleft X value

;---find the upper left corner of smaller, search cube, in secondary
x97=x96-1+x18-x83			;z value
x74=x28-1+x17-x82			;y value
x75=x29-1+x16-x81			;x value

;---both the value to print out what the center of the search
;---cube SHOULD be, but also used to realign the CC results
;---so that it jives with primary space, not partition
x58=x16+(x24+x26)*(x21-1)	;center x in primary  
x59=x17+(x25+x27)*(x22-1)	;y value
x60=(x87+x88)*(x23-1)+x18	;centered z

;---cut out searching in, larger region
WI
[primary]	;cutting from
[slab]		;now called
x24,x25,x87	;dimentions of portion removed
x29,x28,x96	;where taken from

;---set x71-73 to size of  slab
FI x71, x72, x73
[slab]
(12,2,1)

SQ		;square volume searching IN
[slab]
_1

FT		;fourier transform of (primary^2)
_1
[ft_slab_sq]
DE
_1

FT		;fourier transform of slab image
[slab]
[ft_slab]

;---calculate size difference between partition and primary
;---used for windowing after WUrtzel
x34=INT(x71-x80)+1		;x difference, =Bx86
x35=INT(x72-x11)+1		;y difference
x36=INT(x73-x12)+1		;z difference

;---calculate value to move partition in so that it is centered 
;---in primary volume, used before cross correlation
x37=INT(((x71-x80)/2)+0.5)	;x centered, =Bx36 
x38=INT(((x72-x11)/2)+0.5)	;y centered   
x39=INT(((x73-x12)/2)+0.5)	;z centered

;---partition searching for slab
WI			
[secondary]
_1			;output
x80,x11,x12		;dimentions of portion removed
x75,x74,x97		;where to take from

;---create mask using seperated partition
TH M			
_1			;input
_2			;output
B			;below threshold, zeros
x91			;threshold value

;---find number of nonzero pixels in mask
;---my x19 = Bimal's x50
FS x30,x30,x19,x30	;returning avg of pixels
_2

x20=x80*x11*x12 	;size of partition in 3D =Bx26
x31=x19*x20		;# non-0 pixels =Bx50

;--error message
IF (x31.LE.0) THEN
 VM
 echo 'no usable pixels in mask, change x91'
 EN
ENDIF

;---my x32 = Bimal's x51
x32=1.0/x31		;for speed, i dont know why

;---create volume same size as slab, with only search partition in it 
PD
_2			;mask created with TH M, from WI
_3			;output--PADDED MASK---
x71,x72,x73		;size
N			;background not average
(0)			;value of background
(1,1,1)			;where to place mask, shouldnt matter where

;---FT on padded mask 
FT
_3			;input
_4			;output FT(padded mask)

;---multiply FT of primary with FT of padded mask
MU M
[ft_slab]		;input
_4			;input
_3			;output
*			;end multiplication

;--find inverse of FT[primary]* FT(padded mask)
FT
_3			;input, FT[primary]* FT(padded mask)
_6			;ouptut INVERSE FT[primary]* FT(padded mask)

;---normalize the INV(FT[primary]* FT(padded mask))
x33=x32*x32		;used below, my x33 = Bimal's 52
AR
_6 			;input INV (FT[primary]* FT(padded mask)
_6			;output NORM(inv(FT[primary]* FT(padded mask)))
P1*P1*x33		;function to do to each pixel of input

;---ft(primary^2)*complex conjugate of padded mask ?=FTpadded mask?
MU M
[ft_slab_sq]		;input
_4			;input FT of padded mask
_3			;output
*			;end multiplication

;---find inverse of above function
FT
_3			;input
_4			;output INV(ft(primary^2)*FT(padded mask)

;---normalize, with factors
AD F
_4			;input INV(ft(primary^2)*FTpadded mask)
_6			;input NORM ((INV(FT[primary]* FT(padded mask))))
x32, -1.0		;factors to multiply the images by
_4			;output

;---local std. dev.
WU
_4			;input
_4			;output

;---window out local std. dev. info
WI
_4			;input
_9			;output
x34,x35,x36		;amount of space to be removed
(1,1,1)			;where to begin cutting

;--------------begin step 2----------------
;---prepare partition so that avg inside mask =0, and std.
;---dev. inside mask = 1
MM
_2			;mask made from partition, in=BMaskdRotdRef 
_1			;part of secondary we're looking for=BPaddedRefVol
(0)			;background for <0.5

;get avg of masked partition avg=x40, x30 used at garbage
FS x30,x30,x40,x30
_1

;---find sum of all non-0 pixels in masked partition
x41=x40*x20		;non-0 pixels = avg.* how many pixels

;---square masked partition
SQ
_1			;input
_4			;output

;---find average of MASKED(partition))^2
FS x30,x30,x40,x30
_4			;input

;---find sum of all non-0 pixels in MSKD(partition))^2
x42=x20*x40		

;---find std.dev. inside masked (?)
x43=SQRT((x42-(x41*x41*x32))/(x31-1))

;---avg inside mask (?)
x44=x41/x31		;=Bx47
x45=1.0/x43		;=Bx48

;---normalize the masked partition
AR
_1			;input, masked partition
_6			;output
(P1-x44)*x45		;function that defines _11 values

;---trying to streamline the processing
DE 
_3

DE 
_1

;---multiply mask * normalized(masked(Partition))
MM
_2			;mask, from thresholding partition
_6			;input/output file
(0)			;background value

;---create empty volume, in which to paste the MSKD(NORMD(MSKD(partition)))
BL
_4			;output
x71,x72,x73		;size, of the slab volume
N			;not avg background
(0)			;value of background

;---paste in MSKD(NORMD(MSKD(partition)))
IN
_6			;what is getting pasted 
_4			;empty volume
x37,x38,x39		;where to paste, =Bx36

;---take cross corr. of pasted volume with primary volume
CC
[slab]
_4
_6

;---remove pertinent volume
WI
_6			;input, where removing from
_4			;output
x34,x35,x36		;size 
x46,x47,x48		;where to remove from

;---(pertinent CC volume)/(# of non-0 pixels)
AR
_4			:input
_6			;output
P1*x32			;function

;---above/coresponding local std. dev.
MU D
_6			;above volume, input
_9			;local std. dev. volume
_2			;output

;---debugging, used to "focus" search area 
x70 = 1			;x corner
x76 = 1			;y corner
x92 = 1			;z corner 

x93 = x34		;x size 
x94 = x35		;y size 
x95 = x36		;z size 

;---decide if X small enough 
IF (x34.GT.21) THEN
 x93 = 21
 x70 = 1+ INT((x34-x93)/2)
ENDIF 

;---decide if Y small
IF (x35.GT.21) THEN
 x94 = 21
 x76 = 1+ INT((x35-x94)/2)
ENDIF

;---decide if Z small enough
IF (x36.GT.21) THEN
 x95 = 21
 x92 = 1+ INT((x36-x95)/2)
ENDIF

;---look at only a cube of 21 pixels
WI
_2		;input
_6		;output
x93,x94,x95 	;size 
x70,x76,x92	;upper corner to cut from 


;---remove results from loop before
DE
[doc_del]

;---peak search
PK 3D
_6			;input file
+			;search for maxima
(5,0)			;# of peaks, regular center
N			;no gravity
N			;no search box
[doc_del]		;output file

;---want only highedst peak,x52-3 x,y,z values of maxima in 
;---relation to partition, x30 is garbage, x61 value of maxima

  UD S,1,x52,x53,x54,x30,x30,x30,x61
  [doc_del]

  x55=x52		;convert x
  x56=x53		;convert y
  x57=x54		;convert z
  x50=x50+1		;increment write-to line

;---x50 is line of doc out to print to, x58-60 is center  of 
;---partition searched for, x55-7 is center of where match maxima 
;---occurs relative to primary, x61 is maxima value
  SD x50,x58,x59,x60,x55,x56,x57,x61
  [doc_out]

  UD E			;end document send

;---delete the temp file for each partiton
DE
[doc_del]

;---close the loops and the program
    LB1 		;end x loop     
   MY FL		;fluch results
   LB2			;end y loop 

  LB3			;end z loop

;---delete slab files
DE 
[slab]
DE
[ft_slab_sq]
DE
[ft_slab]

;---delete output folder
VM M
cd output
rm -f *.*
cd ..
rmdir output
.

DE



RE 			;return to master program
