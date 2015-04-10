
C++*********************************************************************
C
C SOLV2D.F
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
C	NFAC 	NUMBER OF EIGENVECTORS REQUESTED
C       NUMIM	NUMBER OF IMAGES
C       NACT    NUMBER OF ACTIVE IMAGES
C	NPIX    NUMBER OF PIXELS PER IMAGE
C	LSAV	SEQUENTIAL IMAGE FILE (INPUT)
C	LIMA    IMAGE COORDINATE FILE (OUTPUT)
C	LPIX    PIXEL COORDINATE FILE (OUTPUT)
C	LEIG	EIGENSTUFF FILE (OUTPUT)
C	LTMP	SCRATCH FILE
C
C       ALL FILES MUST BE UNFORMATTED!
C	Q	ARRAY OF FOUR - BYTE WORDS WITH DIMENSION MEM
C	LUV	LOGICAL UTILITY VECTOR (1 = ACTIVE,  0 = INACTIVE)
C 
C **********************************************************************

	SUBROUTINE SOLV2D(A,TST)

	DIMENSION A(2,3)
	COMMON /UNITS/ LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

 

C     I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
      SAVE

	TST=0.
	D=A(1,1)*A(2,2)-A(1,2)*A(2,1)
	P=A(1,3)*A(2,2)-A(2,3)*A(1,2)
	Q=A(2,3)*A(1,1)-A(1,3)*A(2,1)
	IF(D.NE.0.)GO TO 3
	IF(P.EQ.0..AND.Q.EQ.0.)GO TO 7
	GO TO 6

3	IF (ABS(D).GE.(1.E-6)*ABS(P).AND.ABS(D).GE.(1.E-6)*ABS(Q))
     &      GO TO 5

6	CONTINUE
C       ********* DEBUG ***********
	WRITE(NDAT,95)D,P,Q,A(1,1),A(1,2),A(1,3),A(2,1),A(2,2),A(2,3)
95	FORMAT('  D,P,Q :',3(F14.6,3X)/
     1 '   A*,1          A*,2          A*,3'/
     2 '   ',3(F12.4,2X)/3(F12.4,2X))
C       ***********

	TST=1.
7	A(1,3)=0.
	A(2,3)=0.
	RETURN

5	A(1,3)=P/D
	A(2,3)=Q/D
	RETURN
	END
