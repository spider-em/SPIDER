C++*********************************************************************
C
C  CCRD_2.F                      
C              PGI BUG                          FEB 10 2006 ArDean Leith
C              MOD PGI COMPILER BUG             FEB 19 2008 ArDean Leith
C              X RETURNS O                      APR 24 2009 ArDean Leith
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
C  CCRD_2(LUN2,X,BUF, LS,NSAM,NROW)
C
C  PURPOSE: CALCULATES CIRCULAR CROSCORRELATION,  X= F(X*CONJG(Y))
C
C  PARAMETERS:  X      FOURIER TRANSFORMS                 (SENT)
C                      F(X*CONJG(Y))                      (RET.)
C               BUF    I/O BUFFER                         (SENT)
C               LS     NSAM+2-MOD(NSAM,2)                 (SENT)
C               NSAM   READ DIMENSION                     (SENT)
C               NROW   NROW DIMENSION                     (SENT)   
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE CCRD_2(LUN2,X,BUF, LS,NSAM,NROW)

        COMPLEX          :: X(LS/2,NROW)
        COMPLEX          :: BUF(LS/2) 

        COMPLEX          :: CTEMP
        DOUBLE PRECISION :: PI2
 
        NNNN = LS / 2
        LRS  = 2 * NNNN

        PI2  = 8.0 * DATAN(1.0D0)
        ITMP = NSAM / 2
        SX   = PI2 * FLOAT(ITMP) / FLOAT(NSAM)
        ITMP = NROW / 2
        SY   = PI2 * FLOAT(ITMP) / FLOAT(NROW)

        DO J=1,NROW
           IY = J - 1
           IF (IY .GT. (NROW/2))  IY = IY - NROW
           ARGY = SY * IY
           CALL REDLIN(LUN2,BUF,LRS,J)

           DO I=1,NNNN
              ARG    = SX * (I-1) + ARGY
C PGI OP64 BUG?: O(I,J) = X(I,J)*CONJG(BUF(I))*CMPLX(COS(ARG),SIN(ARG))
              CTEMP  = CMPLX(COS(ARG),SIN(ARG))
              X(I,J) = X(I,J) * CONJG(BUF(I)) * CTEMP
           ENDDO
        ENDDO

        INS = -1
        CALL FMRS_2(X,NSAM,NROW,INS)

        END

C      --------------- CCRD_PH_2 -------------------------------------

        SUBROUTINE CCRD_PH_2(LUN2, X,Y, LS,NSAM,NROW)

C       REAL   X((NSAM+2-MOD(NSAM,2)),NROW,NSLICE)
C       ABOVE ON PGI COMPILER 7.1 FAILS TO COMPILE PROPERLY SOMETIMES

        IMPLICIT NONE

        REAL    :: X(LS,NROW)
        REAL    :: Y(LS)

        INTEGER :: LUN2,LS,NSAM,NROW
 
        REAL    :: SX,SY,ARGY,ARG,TMPR,TMPI,FNRM
        INTEGER :: ITMP,J,IY,I,INS
 
	REAL, PARAMETER  :: QUADPI = 3.1415926535897932384
	REAL, PARAMETER  :: PI2    = 2*QUADPI

        ITMP = NSAM / 2
        SX   = PI2 * FLOAT(ITMP) / FLOAT(NSAM)
        ITMP = NROW / 2
        SY   = PI2 * FLOAT(ITMP) / FLOAT(NROW)

C$omp   parallel do private(i,j,iy,argy,arg,tmpr,tmpi,fnrm)
        DO J=1,NROW
           IY = J - 1
           IF (IY .GT. (NROW/2)) IY = IY - NROW
           ARGY = SY * IY

           CALL REDLIN(LUN2,Y,LS,J)

           DO I=1,LS,2
              ARG      = SX * (I-1) / 2 + ARGY

    	      TMPR     = X(I,J)   * Y(I)  + X(I+1,J) * Y(I+1)
	      TMPI     = X(I+1,J) * Y(I)  - X(I,J)   * Y(I+1)

	      X(I,  J) = TMPR * COS(ARG) - TMPI * SIN(ARG)
	      X(I+1,J) = TMPI * COS(ARG) + TMPR * SIN(ARG)

	      FNRM     = 1 / SQRT(X(I,J)**2 + X(I+1,J)**2)

	      X(I,  J) = X(I,  J) * FNRM
	      X(I+1,J) = X(I+1,J) * FNRM
           ENDDO
        ENDDO

        INS = -1
        CALL FMRS_2(X,NSAM,NROW,INS)

        END
