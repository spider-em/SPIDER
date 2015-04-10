C++*********************************************************************
C
C  WRTVOL.F                   NEW MAR 2002 ARDEAN LEITH
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
C  WRTVOL(LUN,NSAM,NROW, ISLICEGO,ISLICEEND,VOL,IRTFLG)
C
C  PARAMETERS:  LUN                IO UNIT                        (SENT)
C               NSAM,NROW          SIZE                           (SENT)
C               ISLICEGO           STARTING SLICE                 (SENT)
C               ISLICEEND          ENDING SLICE                   (SENT)
C               VOL                VOLUME                         (SENT)
C               IRTFLG             ERROR FLAG                     (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE WRTVOL(LUN,NSAM,NROW, ISLICEGO,ISLICEEND,VOL,IRTFLG)

        IMPLICIT NONE

        INTEGER, INTENT(IN)  ::              LUN
        INTEGER, INTENT(IN)  ::              NSAM,NROW
        INTEGER, INTENT(IN)  ::              ISLICEGO,ISLICEEND
        REAL, DIMENSION(NSAM*NROW*(ISLICEEND-ISLICEGO+1)) :: VOL
        INTEGER, INTENT(OUT)  ::             IRTFLG

        INTEGER               :: ILOC,IREC 
        IRTFLG = 1

C       WRITE VOLUME CONTENTS
        ILOC = 1
        DO IREC = (ISLICEGO-1)*NROW+1, ISLICEEND*NROW
           CALL WRTLIN(LUN,VOL(ILOC),NSAM,IREC)
           ILOC = ILOC + NSAM
        END DO

        IRTFLG = 0

        RETURN
        END


