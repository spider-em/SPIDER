C++*********************************************************************
C
C LATCEN.FOR
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
C   LATCEN:
C
C   LATCEN FINDS THE BEST FIT LATTICE UNIT VECTORS (AX,AY) (BX,BY) 
C   GIVEN THE LATTICE CENTER (CX,CY) AND A SET OF TWO OR MORE POINTS
C   WITH LOCATION (SX,SY) AND INDEX (H,K)
C   W WEIGHTS THE RELATIVE IMPORTANCE OF POINTS
C         THE TWO POINTS MUST NOT BE IN A LINE WITH THE CENTER
C
C   TO USE LATCEN:    FIRST CALL LATCEN(IERR,0,,,,,,,,,,,) 
C			[I.E., 13 ARGUMENTS]
C		THIS CLEARS THE VECTORS FOR A NEW CALCULATION
C		:NEXT CALL LATCEN(IERR,1,H,K,SX,SY,W,CX,CY,,,,)
C			TWO OR MORE TIMES
C			ENTER THE POSITIONS OF TWO OR MORE LATTICE
C			POINTS NOT ON A LINE WITH THE CENTER
C		:FINALLY CALL LATCEN(IERR,2,,,,,,,,AX,AY,BX,BY)
C			FOR THE RESULT UNIT VECTORS
C
C **********************************************************************

	SUBROUTINE LATCEN(IERR,INTRY,IH,IK,SX,SY,W,CX,CY,AX,AY,BX,BY)

 


	COMMON/UNITS/ LUN,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
        SAVE          IPOINT,XM,YM
	REAL          XM(2,3),YM(2,3),H,K

	IERR = 0
	IF (INTRY .EQ. 1) GO TO 100
	IF (INTRY .EQ. 2) GO TO 200

C       INIT
	DO  J=1,2
	   DO  I=1,3
	      XM(J,I)=0.
	      YM(J,I)=0.
	   ENDDO
	ENDDO
	IPOINT=0
	RETURN


C     SUM UP POINTS
100   IPOINT=IPOINT+1
      H = IH
      K = IK
	XM(1,3)=XM(1,3)+W*H*(SX-CX)
	XM(2,3)=XM(2,3)+W*K*(SX-CX)
	XM(1,1)=XM(1,1)+W*H*H
	XM(1,2)=XM(1,2)+W*H*K
	XM(2,1)=XM(1,2)
	XM(2,2)=XM(2,2)+W*K*K

	YM(1,3)=YM(1,3)+W*H*(SY-CY)
	YM(2,3)=YM(2,3)+W*K*(SY-CY)
	YM(1,1)=XM(1,1)
	YM(1,2)=XM(1,2)
	YM(2,1)=XM(2,1)
	YM(2,2)=XM(2,2)

	RETURN

C       COMPUTE FIT

200	IF (IPOINT .LT. 2) THEN
	   WRITE(NOUT,219)
	   IERR = 1
219	   FORMAT(' *** THE PROGRAM REQUIRES 2 OR MORE POINTS')
	   RETURN
        ENDIF

210	CALL SOLV2D(XM,TST)
	CALL SOLV2D(YM,TST)
	IF (TST .GE. 0.5) THEN
	   IERR = 1
	   WRITE(NOUT,239)
239	   FORMAT(' *** SINGULARITY IN EQUATIONS FOUND')
	   RETURN
        ENDIF

        AX=XM(1,3)
	BX=XM(2,3)
	AY=YM(1,3)
	BY=YM(2,3)

	END
