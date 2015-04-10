C++*********************************************************************
C
C APPLYW.F
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
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C--*********************************************************************

	SUBROUTINE  APPLYW(CIRC,LCIRC,NUMR,WR,NRING,MAXRIN)

	INTEGER  NUMR(3,NRING),MAXRIN,NUMR3I,NUMR2I
	DIMENSION  CIRC(LCIRC),WR(NRING)

c$omp parallel do private(i,j,jc,NUMR3I,NUMR2I,W)
	DO I=1,NRING
	   NUMR3I=NUMR(3,I)
	   NUMR2I=NUMR(2,I)
	   W=WR(I)
	   CIRC(NUMR2I)=CIRC(NUMR2I)*W
	   IF (NUMR3I.EQ.MAXRIN)  THEN
              CIRC(NUMR2I+1)=CIRC(NUMR2I+1)*W
	   ELSE
              CIRC(NUMR2I+1)=CIRC(NUMR2I+1)*0.5*W
	   ENDIF

	   DO J=3,NUMR3I
              JC=J+NUMR2I-1
              CIRC(JC)=CIRC(JC)*W
           ENDDO
        ENDDO
	END
