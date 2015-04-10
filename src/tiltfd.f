
C++*********************************************************************
C
C TILTFD.FOR
C                                       
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C
C TILTFD(ARX,ARY,BRX,BRY,ATX,ATY,BTX,BTY,GAMR,GAMT,PHI,THA,GAMT2,IER)
C
C  COMPUTES THE TILT ANGLES CORRESPONDING TO A RECIPROCAL SPACE
C  PLANE INTERCEPTING THE SPIKES OF THE FOURIER TRANSFORM OF A 
C  PLANE REAL SPACE LATTICE.  TILTFD IS USED TO FIND THE IMAGE
C  PROJECTION DIRECTION.
C
C   INPUT:
C	ARX,ARY,BRX,BRY	: THE UNIT VECTORS OF THE AXIAL VIEW XTAL SPOTS
C	ATX,ATY,BTX,BTY : THE UNIT VECTORS OF THE XTAL PROJECTION WHO'S
C		PROJECTION ANGLES ARE TO BE FOUND
C   OUTPUT:
C	GAMR: ANGLE BETWEEN THE REFERENCE INPUT VECTORS, A TO B
C	GAMT: ANGLE BETWEEN THE TILT PROJECTION'S VECTORS
C	PHI: ASMUTHAL ANGLE BETWEEN THE AR VECTOR AND PROJECTION PLANE
C	THA: TILT ANGLE, BETWEEN AXIS AND NORMAL TO THE PROJECTION PLANE
C	GAMT2: COMPUTED ANGLE BETWEEN AT AND BT USING THA AND PHI
C		DETERMINED FROM VECTOR LENGTHS
C  ... NOTE THAT
C	THE COMPUTATION ASSUMES THE TWO SETS OF UNIT VECTORS ARE
C	TO THE SAME PAIR OF RECIPROCAL SPACE SPIKES. THAT IS,THE 
C	INDEXING INPUT TO (****) IS CONSISTENT. THE COMPUTATION OF
C	GAMT2 IS USE TO CHECK THIS CONSISTENCY.
C	  	TILTFD ALSO ASSUMES THE ANGLE BETWEEN THE REFERENCE
C	AXIAL PROJECTION VECTORS IS 90. DEGREES. THIS IS CHECKED WITH
C	GAMR.
C
C	IER = 0 NO ERROR
C	IER = 1 ERROR MESSAGE FROM ARCTRG CALL
C
C--*******************************************************************

	SUBROUTINE TILTFD(ARX,ARY,BRX,BRY,ATX,ATY,BTX,BTY,
     &                    GAMR,GAMT,PHI,THA,GAMT2,IER)

C       COMPUTE VECTOR LENGTHS AND THE ANGLES GAMT AND GAMR
 

C       I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
        SAVE

	REAL LRA,LRB,LTA,LTB
	COMMON/UNITS/LUN,NIN,NOUT

	DATA PI/3.14159/

	IER = 0
	LRA=SQRT(ARX*ARX+ARY*ARY)
	LRB=SQRT(BRX*BRX+BRY*BRY)
	LTA=SQRT(ATX*ATX+ATY*ATY)
	LTB=SQRT(BTX*BTX+BTY*BTY)
	COSGT=(ATX*BTX+ATY*BTY)/(LTA*LTB)
	COSGR=(ARX*BRX+ARY*BRY)/(LRA*LRB)
C       CALL ARCTRG(COSGR,DUM,GAMR,IER) !REPLACED BY ACOS BY NAIK 10/2/86
        GAMR = ACOS(COSGR)
C	IF(IER.NE.0)RETURN
C       CALL ARCTRG(COSGT,DUM,GAMT,IER) !      NAIK 10/2/86
        GAMT = ACOS(COSGT)
C	IF(IER.NE.0)RETURN

	WRITE(5,1)COSGR,GAMR,180*GAMR/PI,COSGT,GAMT,180*GAMT/PI
1	FORMAT(6F10.2)

	RA=(ATX*ATX+ATY*ATY)/(ARX*ARX+ARY*ARY)
	RB=(BTX*BTX+BTY*BTY)/(BRX*BRX+BRY*BRY)

C       COMPUTE PHI AND THA

	SNSPHI=(RA-1.)/(RA+RB-2.)
	D=1.
	IF(GAMR.GT.GAMT)D=-1.
	S=D*SQRT(SNSPHI)
C       CALL ARCTRG(S,PHI,DUM,IER)       ! NAIK 10/2/86
C	IF(IER.NE.0)RETURN
        PHI = ASIN(S)

	CSSTHA=(RB-1.)/(RB*RB+RB*RA-2.*RB-RA+1.)
	S=SQRT(CSSTHA)
C       CALL ARCTRG(S,DUM,THA,IER)      ! NAIK 10/2/86
C	IF(IER.NE.0)RETURN
        THA = ACOS(S)

	DN=0.5/SQRT(RA*RB)
	S=DN*SIN(2*PHI)*(1.-1./CSSTHA)
	IF(ABS(S).LT.1.)GOTO 110
	WRITE(NOUT,109)RA,RB,PHI,CSSTHA,S
109	FORMAT(' ** TILTFD ERROR, COS(GAMAT2)>1.'/
     &         ,'RA,RB,PHI,CSSTHA,S'/'   ',5(F10.3))
C110	CALL ARCTRG(S,DUM,GAMT2,IER)      ! NAIK 10/2/86
  110   GAMT2 = ACOS(S)
	RETURN
	END
