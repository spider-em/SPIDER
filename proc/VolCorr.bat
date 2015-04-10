; VolCorr---master
; find correlation numbers for 2 different tilt volumes, 11/4/03 jamie lebarron
; thanks to Ardean Leith for help, and Bimal Rath for help and core loop
; passed large volume values to VCWH for explorer, 12/18/03

; To Test: use same volume as both primary and secondary, identical xyz #'s

; This program assumes: no rotation, that the two volumes are already 
; rotated so they match,both volumes at same resolution, all border white noize 
; removed, that VCWH in same directory

; outputs : 
; x24,x25,x87 =size of searched in partitions, x26,x27,x88 =space between  
; partitions, x77-9= number of partitons/search cubes, x80,x11,x12= total size 
; ofsearch cubes, x16-8 =center of partition/search cubes, x91 = pixel masking 
; value, pixels less that this will be set to zero

;-----------------------------ALL INPUTS----------------------------
;---input volume file names, must be orienteated the same way---
;---if not, must use if statements, and rotations to fix

FR G 
[primary]wth5sin001	;the volume that will be searched in
FR G
[secondary]wth5sec001		;volume that will be searched for 
;---end input file name---

;---x91 pixel mask 
x91 = -15000

;---- number of search volumes, must be odd, i'll check
x77 = 5		;number of divisions in x
x78 = 5		;number of divisions in y
x79 = 3		;number of divisions in z

;--- size of cubes, must be odd
x80 = 51		;x
x11 = 51		;y
x12 = 31		;z
;------------------------END ALL INPUTS--------------------------

;---test to make sure user inputs are odd--------
x95=x77
x96=INT(x95/2)
 IF((x96*2).EQ.x95) THEN
 VM
 echo "something was not odd, change input values"
 EN
 ENDIF
x95=x78
x96=INT(x95/2)
 IF((x96*2).EQ.x95) THEN
 VM
 echo "something was not odd, change input values"
 EN
 ENDIF
x95=x79
x96=INT(x95/2)
 IF((x96*2).EQ.x95) THEN
 VM
 echo "something was not odd, change input values"
 EN
 ENDIF
x95=x80
x96=INT(x95/2)
 IF((x96*2).EQ.x95) THEN
 VM
 echo "something was not odd, change input values"
 EN
 ENDIF
x95=x11
x96=INT(x95/2)
 IF((x96*2).EQ.x95) THEN
 VM
 echo "something was not odd, change input values"
 EN
 ENDIF
x95=x12
x96=INT(x95/2)
 IF((x96*2).EQ.x95) THEN
 VM
 echo "something was not odd, change input values"
 EN
 ENDIF
;----end odd test for inputs-------------------------

;-------determine size of primary and secondary volumes----
FI x93,x94,x95	;x, y, z max for primary 
[primary]
12,2,1			;x, y, z max for primary 

FI x74, x75, x76  	;x, y, z max for secondary
[secondary]
12,2,1 			;x, y, z max for secondary
;-------size of volumes found-------

;---test to see if same size
IF(x93.NE.x74) THEN
VM
echo "volumes not same size"
END
ENDIF

IF(x94.NE.x75) THEN
VM
echo "volumes not same size"
END
ENDIF

IF(x95.NE.x76) THEN
VM
echo "volumes not same size"
END
ENDIF

;---since x93-5, is equal to x74-6
x74=0
x75=0
x76=0

;---needed for large to be searched in regions
x87=INT(x95/x79)	  	;z slab height
x30=x87
x31=INT(x30/2)
IF((x31*2).EQ.x30) x87=x87-1	;ensure that slab height is odd

x24=INT(x93/x77)		;x region height
x30=x24
x31=INT(x30/2)
IF((x31*2).EQ.x30) x24=x24-1

x25=INT(x94/x78)		;y region width
x30=x25
x31=INT(x30/2)
IF((x31*2).EQ.x30) x25=x25-1

;---center jump value,used once per region, center
x16=INT((x24/2)+0.5)	;x center
x17=INT((x25/2)+0.5)	;y center
x18=INT((x87/2)+0.5)	;z center

;---error test to insure that the distance between 
;---centers is greater than width of serach cubes
;---error prints out 1.5* center jump value, this is intentional
x30=x16*1.5
IF(x80.GT.x30) THEN
 VM M
 echo "the size of the x search cube was larger than jump value"
 "\n                 "     {***x80} "                  "{***x30}
.
 EN
 ENDIF

x30=x17*1.5
IF(x11.GT.x30) THEN
 VM M
 echo "the size of the y search cube was larger than jump value"
 "\n                 "     {***x11} "                  "{***x30}
.
 EN
 ENDIF

x30=x18*1.5
IF(x12.GT.x30) THEN
 VM M
 echo "the size of the z search cube was larger than jump value"
"\n                 "     {***x12} "                   "{***x30}
.
 EN
 ENDIF
;---end error test

;---calculates pixels not included in search cube
x13 = x93 - (x80 * x77)	
x14 = x94 - (x11 * x78)
x15 = x95 - (x12 * x79)

;---space between looking in regions
IF(x79.GT.1) x88=INT((x95-(x87*x79))/(x79-1))	;z padding
IF(x78.GT.1) x27=INT((x94-(x25*x78))/(x78-1))	;y space
IF(x77.GT.1) x26=INT((x93-(x24*x77))/(x77-1))	;x space

;---print out how much volume is being serached for
VM M
echo "\nNot using " {****x13} "of " {****x93} "x for cubes"
 "\nNot using " {****x14} "of " {****x94} "y for cubes"
 "\nNot using " {****x15} "of " {****x95} "z for cubes"
 "\n In the X dim, your partition is "{***x24}" and your cube is "{***x80}
 "\n In the Y dim, your partition is "{***x25}" and your cube is "{***x11}
 "\n In the Z dim, your partition is "{***x87}" and your cube is "{***x12}
.

;---deal with created files and directories
;---input files ARE NOT DELT WITH HERE
FR G
[ft_slab_sq]output/FtSlabSq

FR G
[ft_slab]output/FtSlab

FR G
[slab]output/Slab

FR G
[doc_del]output/DocDel

FR G
[doc_out]VCOutput

FR G
[output]output

;---make output dir in UNIX
VM
mkdir -p [output] 
;---end dealing with files

;---call the workhorse program
@VCWH[x24,x25,x26,x27,x77,x78,x79,x80,x11,x12,x16,x17,x18,x87,x88,x91]

EN		;end program

