C++*********************************************************************
C
C  REDVOL.F       NEW AS GETVOLDAT               MAR 1999 ArDean Leith
C                 NO LONGER ALLOCATES            DEC 2000 ArDean Leith
C                 ADDED ISLICEGO                 MAY 2002 ArDean Leith
C                 REDVOL_SEL                     NOV 2011 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  REDVOL(LUN,NX,NY,ISLICEGO,ISLICEEND,VOL,IRTFLG)
C
C  PARAMETERS:  LUN                IO UNIT                        (SENT)
C               NX,NY              SIZE                           (SENT)
C               ISLICEGO           STARTING SLICE                 (SENT)
C               ISLICEEND          ENDING SLICE                   (SENT)
C               VOL                VOLUME                         (RET.)
C               IRTFLG             ERROR FLAG                     (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE REDVOL(LUN,NX,NY,ISLICEGO,ISLICEEND,VOL,IRTFLG)

        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: LUN
        INTEGER, INTENT(IN)  :: NX,NY
        INTEGER, INTENT(IN)  :: ISLICEGO,ISLICEEND
        REAL                 :: VOL(NX*NY*(ISLICEEND-ISLICEGO+1))
        INTEGER, INTENT(OUT) :: IRTFLG

        INTEGER              :: ILOC,IREC

        IRTFLG = 1

C       RECOVER VOLUME CONTENTS
        ILOC = 1
        DO IREC = (ISLICEGO-1)*NY+1, ISLICEEND*NY
           CALL REDLIN(LUN,VOL(ILOC),NX,IREC)
           ILOC = ILOC + NX
        ENDDO

        IRTFLG = 0

        END

C     ----------------------------- REDVOL_SEL  ------------------------

C       THE SAME AS REDVOL, CAN SPECIFY IF MPI_BCAST USED.

        SUBROUTINE REDVOL_SEL(LUN,NX,NY,ISLICEGO,ISLICEEND,
     &                        MPIBCAST, VOL,IRTFLG)

        IMPLICIT NONE

        INTEGER, INTENT(IN)  :: LUN
        INTEGER, INTENT(IN)  :: NX,NY
        INTEGER, INTENT(IN)  :: ISLICEGO,ISLICEEND
        REAL                 :: VOL(NX*NY*(ISLICEEND-ISLICEGO+1))
        LOGICAL, INTENT(IN)  :: MPIBCAST
        INTEGER, INTENT(OUT) :: IRTFLG

        INTEGER              :: ILOC,IREC

C       RECOVER VOLUME CONTENTS
        ILOC = 1
        DO IREC = (ISLICEGO-1)*NY+1, ISLICEEND*NY
           CALL REDLIN_SEL(LUN,NX,IREC,MPIBCAST, VOL(ILOC),IRTFLG)
           ILOC = ILOC + NX
        ENDDO

        END

