
C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C
C  LATICE(INTRY,IH,IK,SX,SY,W,CX,CY,AX,AY,BX,BY,ANGLE)
C
C  LATICE FINDS THE BEST FIT LATTICE UNIT VECTORS (AX,AY) (BX,BY) AND
C  LATTICE CENTER (CX,CY) FROM A SET OF THREE OR MORE POINTS
C  WITH LOCATION (SX,SY) AND INDEX (H,K)
C  W WEIGHTS THE RELATIVE IMPORTANCE OF POINTS
C
C  THREE OR MORE POINTS MUST NOT BE ON THE SAME LINE
C  TO USE LATICE:FIRST CALL LATICE(0,,,,,,,,,,,)
C			THIS CLEARS THE VECTORS FOR A NEW CALCULATION
C		:NEXT CALL LATICE(1,H,K,SX,SY,W,,,,,,)
C			THREE OR MORE TIMES
C			ENTER THE POSITIONS OF THREE OR MORE LATTICE
C			POINTS NOT ON A LINE
C		:FINALLY CALL LATICE(2,,,,,,CX,CY,AX,AY,BX,BY)
C			FOR THE RESULT UNIT VECTORS
C
C**********************************************************************

	SUBROUTINE LATICE(INTRY,IH,IK,SX,SY,W,CX,CY,AX,AY,BX,BY,ANGLE)

C       I THINK SAVE IS NEEDED FEB 99 al
        SAVE

	REAL H,K
	REAL XXM(3),YYM(3)

	COMMON/UNITS/ LUN,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
	COMMON/POINTR/IPOINT
	COMMON/LATTRX/XM(3,4),YM(3,4)
	DATA PI/3.14159/

	IF(INTRY.EQ.1)GO TO 100
	IF(INTRY.EQ.2)GO TO 200

C  INIT
	DO  J=1,3
           DO  I=1,4
             XM(J,I)=0.
             YM(J,I)=0.
	   ENDDO
	ENDDO
	IPOINT=0
	RETURN

C   SUM UP POINTS
100	IPOINT=IPOINT+1
C	WRITE(NOUT,711)IPOINT
C711	FORMAT( ' IPOINT =',1X,I5)
      H = IH
      K = IK
	XM(1,4)=XM(1,4)+W*SX
	XM(1,1)=XM(1,1)+W
	XM(1,2)=XM(1,2)+W*H
	XM(1,3)=XM(1,3)+W*K
	XM(2,4)=XM(2,4)+W*H*SX
	XM(2,1)=XM(1,2)
	XM(2,2)=XM(2,2)+W*H*H
	XM(2,3)=XM(2,3)+W*H*K
	XM(3,4)=XM(3,4)+W*K*SX
	XM(3,1)=XM(1,3)
	XM(3,2)=XM(2,3)
	XM(3,3)=XM(3,3)+W*K*K

	YM(1,4)=YM(1,4)+W*SY
	YM(2,4)=YM(2,4)+W*H*SY
	YM(3,4)=YM(3,4)+W*K*SY
	YM(1,1)=XM(1,1)
	YM(1,2)=XM(1,2)
	YM(1,3)=XM(1,3)
	YM(2,1)=XM(2,1)
	YM(2,2)=XM(2,2)
	YM(2,3)=XM(2,3)
	YM(3,1)=XM(3,1)
	YM(3,2)=XM(3,2)
	YM(3,3)=XM(3,3)
	XXM(1)=XM(1,4)
	XXM(2)=XM(2,4)
	XXM(3)=XM(3,4)
	YYM(1)=YM(1,4)
	YYM(2)=YM(2,4)
	YYM(3)=YM(3,4)

C	WRITE(NOUT,711)IPOINT
C	DO 117 I=1,3
C          WRITE(NOUT,712)I,(XM(I,J),J=1,4)
C117	CONTINUE
C712	FORMAT( ' XM(I,J),J=1,4 ; I= ',I2,2X,4(F9.3,1X))
	RETURN

C  COMPUTE BEST FIT
200	IF(IPOINT.GE.3)GO TO 210
	WRITE(NOUT,212)IPOINT
212	FORMAT(' ONLY',I3,'  POINTS ,THE PROGRAM NEEDS THREE OR MORE')
	RETURN

210	CONTINUE
	CALL SOLVE(XM,3,3)
	CALL SOLVE(YM,3,3)
	CX=XM(1,4)
	AX=XM(2,4)
	BX=XM(3,4)
	CY=YM(1,4)
	AY=YM(2,4)
	BY=YM(3,4)
	IF(AY.NE.0.)ANGLEA=ATAN2(AX,AY)*180./PI
	IF(AY.EQ.0.)ANGLEA=SIGN(90.,AX)
	IF(BY.NE.0.)ANGLEB=ATAN2(BX,BY)*180./PI
	IF(BY.EQ.0.)ANGLEB=SIGN(90.,BX)
	ANGLE=ANGLEA-ANGLEB
	XM(1,4)=XXM(1)
	XM(2,4)=XXM(2)
	XM(3,4)=XXM(3)
	YM(1,4)=YYM(1)
	YM(2,4)=YYM(2)
	YM(3,4)=YYM(3)

C      WRITE(NOUT,999)CX,CY,AX,AY,BX,BY
C999   FORMAT(1X,6F10.2)
C	WRITE(NOUT,711)IPOINT

	RETURN
	END
