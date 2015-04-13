[x11,x22]
;
; Input [x11] : spatial frequency in Angstroms,
; Output [x22]: corresponding intensity (from 2nd column in scatter doc file)
;
; The filter is scaled to have 1 at 100A

FR
?X-ray scattering power spectrum?<scatter>

UD N,x53   ; number of point in scattering file
<scatter>

X52 = INT(x53/2)

DO LB1 X51=1,x53

UD IC,X52,X31,X22
<scatter>

X52=X52+1

IF(X52.gt.1500) GOTO LB2

UD IC,X52,X32,X23
<scatter>

IF(X11.LE.X31) THEN
   IF(X11.GE.X32)  THEN
      X77=(X31-X11)**2
      X78=(X32-X11)**2
      IF(X78.LT.X77) X22=X23
      GOTO LB2
   ENDIF
ENDIF

IF(X11.LT.X32) THEN
   X51=X52
   X52=INT((X53-X52)/2+X52)
ELSE
   X53=X52
   X52=INT((X52-X51)/2+X51)
ENDIF

LB1

LB2
X22=X22/37891.

RE
