
C++*******************************************************************
C
C DENOV.FOR
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
C  DENOV(LUN1,LUN2,NSAM,NROW
C
C         LUN1    LOGICAL UNIT NUMBER OF FILE
C         LUN2    LOGICAL UNIT NUMBER OF FILE
C         NSAM    NUMBER OF SAMPLES
C         NROW    NUMBER OF ROWS
C         FMAX    MAXIMUM
C         FMIN    MINIMUM
C
C
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C--*******************************************************************

        SUBROUTINE DENOV(LUN1,LUN2,NSAM,NROW,NSLICE)

        INCLUDE 'CMBLOCK.INC'

        DIMENSION DARRAY(NSAM)

        IF (IMAMI.NE.1) CALL NORM3(LUN1,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
        DIFF = FMAX - FMIN
        IF (DIFF .EQ. 0.0) THEN
           CALL ERRT(101,'ZERO DENSITY RANGE',NE)
           RETURN
        ENDIF

        WRITE(NOUT,20)FMIN,FMAX,FMIN
20      FORMAT(' FMIN = ',G10.3,'  FMAX = ',G10.3, '  OFFSET = ',G10.3)
        CALL RDPRM(SCAL,NOT_USED,'SCALING FACTOR')
        CALL RDPRM(OFFS,NOT_USED,'OFFSET')
        OFF = FMIN + OFFS

        DO I=1,NROW*NSLICE
           CALL REDLIN(LUN1,DARRAY,NSAM,I)
           DO  K=1,NSAM
              DARRAY(K)=AMOD((DARRAY(K)-OFF)*SCAL,DIFF)
           ENDDO
           CALL WRTLIN(LUN2,DARRAY,NSAM,I)
        ENDDO

        RETURN
        END

