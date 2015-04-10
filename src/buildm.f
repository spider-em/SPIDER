C++*********************************************************************
C
C  BUILDM.F         MERGED WITH REANG                JUL 03 ARDEAN LEITH
C                   BYLIST ADDED                     SEP 03 ARDEAN LEITH
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
C  BUILDM(ILIST,DM,NANG,ANGBUF,FILLSS,SS,BYLIST,IRTFLG)
C
C  PURPOSE: BULID ROTATION MATRICES FROM THREE EULERIAN ANGLES
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE BUILDM(ILIST,DM,NANG,ANGBUF,FILLSS,SS,BYLIST,IRTFLG)

        INCLUDE 'CMBLOCK.INC' 

        DIMENSION  :: DM(9,NANG),ILIST(NANG),ANGBUF(4,NANG),SS(6,NANG)
        LOGICAL    :: FILLSS,BYLIST

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       READ ANGLES FROM THE DOCUMENT FILE.
C       ORDER IN THE DOCUMENT FILE IS PSI, THETA, PHI AND ANGLES 
C       ARE IN DEGREES! IN ANG ARRAY IT IS OTHER WAY AROUND
C       OUTPUT IS COMPACTED TO 1...NANG LINES (NOT BY SELECTOR)

        DO K=1,NANG

C          GET ANGLE SELECTOR INDEX FROM ILIST
           ITMP   = ILIST(K)

           ICOUNT = ANGBUF(1,ITMP)
           IF (ICOUNT .LE. 0) THEN
C             MISSING KEY
              CALL ERRT(102,'MISSING ANGLE FOR IMAGE',ITMP)
              IRTFLG = 1
              RETURN
           ENDIF

           KT = K
           IF (BYLIST) KT = ITMP

           CALL CANG(ANGBUF(4,ITMP),ANGBUF(3,ITMP),ANGBUF(2,ITMP),
     &               FILLSS,SS(1,KT),DM(1,KT))

           IF (VERBOSE) THEN
              IF (MYPID .LE. 0)WRITE(NOUT,333)K,(ANGBUF(J,ITMP),J=2,4)
333           FORMAT('  PROJECTION #',I7,
     &               '; PSI=',F6.1,' THETA=',F6.1,' PHI=',F6.1)
           ENDIF
         ENDDO

        IRTFLG = 0
        END

