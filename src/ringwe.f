
C ++********************************************************************
C                                                                      *
C RINGWE                                                                     *
C                                                                      *
C               MODIFED FOR USING FFTW3           MAR 2003 ARDEAN LEITH
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
C  RINGWE_NEW(WR,NUMR,NRING,MAXRIN)                                                                   *
C                                                                      *
C  PURPOSE:  FINDS WEIGHTS FOR RADIAL X-CORRELATION RINGS                                                            *
C                                                                      *
C  PARAMETERS:
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE RINGWE_NEW(WR,NUMR,NRING,MAXRIN)

	REAL, INTENT(OUT)   :: WR(NRING)
	INTEGER, INTENT(IN) :: NUMR(3,NRING), MAXRIN, NRING

	PI = 8.0*DATAN(1.0D0)

	DO I=1,NRING

#ifdef SP_LIBFFTW3
           LEN = NUMR(3,I) - 2
#else
           LEN = NUMR(3,I)
#endif
	   WR(I) = REAL(NUMR(1,I)) * PI / REAL(LEN) *
     &             REAL(MAXRIN) / REAL(LEN)
	ENDDO
	END

C       -------------------- RINGWE --------------------------------

	SUBROUTINE RINGWE(WR,NUMR,NRING,MAXRIN)

	REAL, INTENT(OUT)   :: WR(NRING)
	INTEGER, INTENT(IN) :: NUMR(3,NRING),MAXRIN

	PI = 8.0 * DATAN(1.0D0)

	DO I=1,NRING
           LEN   = NUMR(3,I)
	   WR(I) = REAL(NUMR(1,I)) * PI / REAL(LEN) *
     &             REAL(MAXRIN) / REAL(LEN)
	ENDDO
	END
