C++*********************************************************************
C
C APPLYWS.F
C               MODIFED FOR USING FFTW3           MAR 2003 ARDEAN LEITH
C               WEIGHTING BUG FOR POS. 2          JUN 2010 ARDEAN LEITH
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
C APPLYWS_NEW(CIRC,LCIRC,NUMR,WR,NRING)
C
C PURPOSE: APPLY WEIGTHS FROM ARRAY: WR TO: CIRC ARRAY VALUES
C	   WR()  = REAL(RING_RADIUS * PI / REAL(NVAL) *
C                  REAL(MAXRIN) / REAL(NVAL)
C
C NOTE:  FROM SPIDER ?--> 16.0  (IN 19MAY08) THIS WORKED WITH NON-
C        FFTW FOURIER FORMATTING WHICH HAD COMPACTED FOURIER WITH
C        SECOND LOCATION ON EACH RING REPRESENTING FINAL FFT TERM. 
C        THIS SECOND LOCATION HAD ONLY HALF WEIGHT.  
C        AFTER SWITCH TO FFTW THE SECOND LOCATION KEPT THE OLD
C        INCORRECT WEIGHT UP TILL:   2010/06/24 WHEN IT WAS FIXED
C        BY WEIGHTING THE SECOND TO LAST LOCATION OF EACH RING
C        THIS ACCOUNTS FOR CHANGES IN VALUES OF CCROT ON THOSE DATES.
C        IN JULY 2011 CHANGE IN APSH_** FOR CCROT SO THAT NEGATIVE 
C        VALUES (WHICH ARE NORMAL) AGAIN WERE REPORTED CORRECTLY.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE APPLYWS_NEW(CIRC,LCIRC,NUMR, WR,NRING)

      IMPLICIT NONE

      REAL    :: CIRC(LCIRC)
      INTEGER :: LCIRC
      INTEGER :: NUMR(3,NRING)
      INTEGER :: NRING
      REAL    :: WR(NRING)

      INTEGER :: IRING,IGO,NVAL,MAXRIN,J
      REAL    :: WT

      DO IRING=1,NRING
         IGO    = NUMR(2,IRING)
         NVAL   = NUMR(3,IRING)   ! INCLUDING FFT PAD
         MAXRIN = NUMR(3,NRING)   ! LENGTH OF LONGEST RING + FFT PAD

         WT     = WR(IRING)

C        APPLY WEIGHTS FOR CIRC LOCATIONS

         DO J=IGO,IGO+NVAL-1      ! LOOP OVER ALL FFT COEFFs
            CIRC(J) = CIRC(J) * WT
         ENDDO

         IF (NVAL .NE. MAXRIN) THEN
C           IF RING LENGTH < THAN THE MAX. RING LENGTH THEN WEIGHT
C           OF LAST REAL FFT COEF FOR THE RING IS HALF USUAL WEIGHT.
C           WHY? al
            CIRC(IGO+NVAL-2) = CIRC(IGO+NVAL-2) * 0.5 
         ENDIF
      ENDDO

      END


C       -------------------- APPLYWS --------------------------------
C       STILL USED IN: oracfmsk.f, oracfmskm.f:         

        SUBROUTINE APPLYWS(CIRC,LCIRC,NUMR,WR,NRING,MAXRIN)

	INTEGER :: NUMR(3,NRING), MAXRIN, NVAL,IGO
	REAL    :: CIRC(LCIRC), WR(NRING)

	DO I=1,NRING
	   IGO       = NUMR(2,I)
	   NVAL      = NUMR(3,I)

	   W         = WR(I)
	   CIRC(IGO) = CIRC(IGO)*W

	   IF (NVAL .EQ. MAXRIN)  THEN
	      CIRC(IGO+1) = CIRC(IGO+1) * W
	   ELSE
	      CIRC(IGO+1) = CIRC(IGO+1) * 0.5 * W
	   ENDIF

	   DO J=3,NVAL
	      JC       = J + IGO - 1
	      CIRC(JC) = CIRC(JC) * W
	   ENDDO
	ENDDO
	END

