C++*********************************************************************
C
C ROTAL3.F
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
C  ROTAL3(LUN1,LUN2,NSAM,NROW,NSLICE,MODE)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
 
        SUBROUTINE ROTAL3(LUN1,LUN2,NSAM,NROW,NSLICE,MODE)

        INCLUDE 'CMBLOCK.INC' 

        REAL, ALLOCATABLE, DIMENSION(:) ::  Q
        CHARACTER(LEN=3)                ::  MODE

        DIMENSION   P1(3),P2(3)
 
        MEMWANT = NSAM*NROW*NSLICE
        ALLOCATE(Q(MEMWANT),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(102,'ROTAL3 -- UNABLE TO ALLOCATE Q',MEMWANT) 
           GOTO 9999
        ENDIF
 
        CALL REDVOL(LUN1,NSAM,NROW,1,NSLICE,Q,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       GET  ANGLE
	CALL  RDPRM(ALPHA,NOT_USED,'ALPHA')

C       GET  ROTATION AXIS
        P1(3)   = HUGE(P1)
        CALL  RDPRM3S(P1(1),P1(2),P1(3),NOT_USED,
     &        'X, Y, & Z FOR FIRST POINT DEFINING ROTATION AXIS',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (P1(3) .EQ. HUGE(P1)) THEN
           CALL  RDPRM1S(P1(3),NOT_USED,
     &             'Z FOR FIRST POINT DEFINING ROTATION AXIS',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

        P2(3)   = HUGE(P2)
        CALL  RDPRM3S(P2(1),P2(2),P2(3),NOT_USED,
     &      'X, Y, & Z FOR SECOND POINT DEFINING ROTATION AXIS',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (P2(3) .EQ. HUGE(P1)) THEN
           CALL  RDPRM1S(P2(3),NOT_USED,
     &             'Z FOR SECOND POINT DEFINING ROTATION AXIS',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

        IF (MODE(3:3) .EQ. 'S') THEN
	   CALL ROTL3(LUN2,Q,NSAM,NROW,NSLICE,P1,P2,ALPHA)
        ELSE
	   CALL ROTL3Q(LUN2,Q,NSAM,NROW,NSLICE,P1,P2,ALPHA)
        ENDIF

9999    IF (ALLOCATED(Q)) DEALLOCATE(Q) 

	RETURN
        END
