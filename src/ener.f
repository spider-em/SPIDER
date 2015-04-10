
C ++********************************************************************
C                                                                      *
C ENER.F                                                               *
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
C
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C --********************************************************************

        DOUBLE PRECISION FUNCTION  ENER(CIRC,LCIRC,NRING,NUMR,MODE)
        DIMENSION  CIRC(LCIRC)
        INTEGER  NUMR(3,NRING)
        DOUBLE PRECISION EN,PI,TQ,ENERT
        CHARACTER*1  MODE
        PI=4.0D0*DATAN(1.0D0)
        IF (MODE .EQ. 'F')  PI=2*PI
        ENERT=0.0             
        DO I=1,NRING
           EN=0.0
           TQ=2.0*REAL(NUMR(1,I))*PI/REAL(NUMR(3,I))
           DO J=NUMR(2,I)+2,NUMR(2,I)+NUMR(3,I)-1
              EN=EN+TQ*CIRC(J)*CIRC(J)
           ENDDO
           EN=EN+TQ*(CIRC(NUMR(2,I))**2+CIRC(NUMR(2,I)+1)**2)/2.0
           ENERT=ENERT+EN/REAL(NUMR(3,I))
        ENDDO
        ENER=ENERT
        END

