C ++********************************************************************
C                                                                      *
C    FOURING.F                                                         *
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
C  PARAMETERS:
C             
C IMAGE_PROCESSING_ROUTINE                                             *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  FOURING(CIRC,LCIRC,NUMR,NRING,EO,MODE)

        INTEGER  NUMR(3,NRING)
        DIMENSION  CIRC(LCIRC)
        CHARACTER*1  MODE
        DOUBLE PRECISION  E,EO,QT,PI


        PI=4.0*DATAN(1.0D0)
        IF(MODE.EQ.'F')  PI=2*PI
        E=0.0
c$omp parallel do private(i,j,nsirt,qt,l),reduction(+:e)
        DO    I=1,NRING
           NSIRT=NUMR(3,I)
           QT=REAL(NUMR(1,I))*PI/REAL(NSIRT)
           DO    J=NUMR(2,I),NUMR(2,I)+NSIRT-1
              E=E+QT*DBLE(CIRC(J))*CIRC(J)
           ENDDO
        ENDDO
        EO=E
c$omp parallel do private(i,l)
        DO    I=1,NRING
           L=LOG2(NUMR(3,I))
           CALL  FFTR_Q(CIRC(NUMR(2,I)),L)
        ENDDO
        END
