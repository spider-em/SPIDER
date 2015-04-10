C++*********************************************************************
C
C CROSRNG_DS.F 
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
C CROSRNG_DS(CIRC1,CIRC2,LCIRC,NRING,MAXRIN,NUMR,QN,TOT,TT)
C
C  INPUT - FOURIER TRANSFORMS OF RINGS!!!
C  CIRC1 ALREADY MULTIPLIED BY WEIGHTS!!
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE CROSRNG_DS(CIRC1,CIRC2,LCIRC,NRING,MAXRIN,NUMR,
     &                        QN,TOT,TT)

	INTEGER           NUMR(3,NRING)
	DIMENSION         CIRC1(LCIRC),CIRC2(LCIRC)
	DOUBLE PRECISION  TT(*),QN

C       AUTOMATIC ARRAYS
	DIMENSION         T(MAXRIN+2)
	DOUBLE PRECISION  Q(MAXRIN+2)
	DOUBLE PRECISION  T7(-3:3)

        DATA  ICALL/0/

	IP = MAXRIN
#ifndef SP_LIBFFT
	IP=-LOG2(IP)
#endif


C       ZERO Q ARRAY
	Q = 0.0D0
     
	T(MAXRIN+1) = 0.0
	T(MAXRIN+2) = 0.0

 	DO I=1,NRING
	   NUMR3I = NUMR(3,I)
	   NUMR2I = NUMR(2,I)

	   T(1) = (CIRC1(NUMR2I)) * CIRC2(NUMR2I)

	   IF (NUMR3I .NE. MAXRIN)  THEN
	      T(NUMR3I+1) = (CIRC1(NUMR2I+1))*CIRC2(NUMR2I+1)
	      T(2)        = 0.0

	      DO J=3,NUMR3I,2
	         JC     = J+NUMR2I-1

	         T(J)   =  (CIRC1(JC))*CIRC2(JC) +
     &                     (CIRC1(JC+1))*CIRC2(JC+1)

	         T(J+1) = -(CIRC1(JC))*CIRC2(JC+1) +
     &                     (CIRC1(JC+1))*CIRC2(JC)
	      ENDDO

	      Q(1:NUMR3I+1) = Q(1:NUMR3I+1) + T(1:NUMR3I+1)

	   ELSE
	      T(2) = CIRC1(NUMR2I+1) * CIRC2(NUMR2I+1)

	      DO J=3,MAXRIN,2
	         JC     =  J+NUMR2I-1

	         T(J)   =  (CIRC1(JC))*CIRC2(JC) +
     &                     (CIRC1(JC+1))*CIRC2(JC+1)

	         T(J+1) = -(CIRC1(JC))*CIRC2(JC+1) +
     &                     (CIRC1(JC+1))*CIRC2(JC)
	      ENDDO
	      Q = Q + T
           ENDIF
	ENDDO
 
#ifdef SP_LIBFFT
	INV         = -1
	LDA         = 1
	Q(MAXRIN+1) = Q(2)
	Q(2)        = 0.0
	Q(MAXRIN+2) = 0.0
	CALL  ZDFFT1DU(INV,IP,Q,LDA,TT)

C       SKIP THE NORMALIZATION, DIVIDE THE MAXIMUM INSTEAD.
C	CALL DSCAL1D(IP,1.0D0/DBLE(FLOAT(IP)),Q,LDA)
#else
	CALL FFTR_D(Q,IP)
#endif

	QN = -1.0D20
	DO J=1,MAXRIN
	   IF (Q(J) .GE. QN)  THEN
	      QN   = Q(J)
	      JTOT = J
	   ENDIF
	ENDDO
#ifdef SP_LIBFFT
	QN=QN/MAXRIN
#endif

	DO K=-3,3
           J     = MOD(JTOT+K+MAXRIN-1,MAXRIN) + 1
	   T7(K) = Q(J)
	ENDDO

	CALL  PRB1D(T7,7,POS)

	TOT = FLOAT(JTOT) + POS

	END
