C++*********************************************************************
C
C FOURING_Q.F
C             DE-OMP PARRALEL LOOP FOR FFTW3 USE    MAR 08 ARDEAN LEITH
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
C FOURING_Q(CIRC,LCIRC,NUMR,NRING,EO,MODE)
C
C PARAMETERS:
C             CIRC   - FT OF RINGS MULTIPLIED BY WEIGHTS   (SENT/RET)
C             EO  WEIGHTING?                                     RET. 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE FOURING_Q(CIRC,LCIRC,NUMR,NRING,EO,MODE)

        INTEGER           :: NUMR(3,NRING)
        DIMENSION         :: CIRC(LCIRC)
        CHARACTER*1       :: MODE
        DOUBLE PRECISION  :: E,EO,QT,PI

        PI = 4.0 * DATAN(1.0D0)
        IF (MODE .EQ. 'F')  PI = 2*PI
        E = 0.0

c$omp   parallel do private(i,j,nval,igo,qt),reduction(+:e)
        DO  I=1,NRING
           NVAL = NUMR(3,I) - 2
           IGO  = NUMR(2,I)
           QT   = REAL(NUMR(1,I)) * PI / REAL(NVAL)  ! PI * N / (N-2)

           DO J=IGO, IGO+NVAL-1
              E = E + QT * DBLE(CIRC(J)) * CIRC(J)
           ENDDO
        ENDDO

        EO = E

c       parallel do private(i,inv,nval,igo)(FAILS MAKING PLAN WITH FFTW3)
        DO I=1,NRING
           INV  = +1              ! FMRS USES INV AS ERROR RETURN
           NVAL = NUMR(3,I) - 2
           IGO  = NUMR(2,I)
           CALL FMRS_1(CIRC(IGO),NVAL,INV)
        ENDDO

        END

