
C++*********************************************************************
C
C   TILT.F
C                    USED OPFILE                  NOV 00 ARDEAN LEITH
C                    OPFILEC                      FEB  03 ARDEAN LEITH
C                    RDPRAF REMOVED               DEC  05 ARDEAN LEITH 
c
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
C   TILT(LUN)
C
C    PARAMETERS:     LUN     LOGICAL UNIT NUMBER OF FILE
C
C--*********************************************************************

	SUBROUTINE TILT(LUN)

        INCLUDE 'CMBLOCK.INC'  
        INCLUDE 'LABLOCK.INC'  
        INCLUDE 'CMLIMIT.INC' 

        INTEGER               :: NUMBER(200)
        CHARACTER(LEN=MAXNAM) :: FILNAM,FILPAT
        CHARACTER(LEN=1)      :: NULL
        REAL                  :: VALUES(3),FWA(4)
        REAL, PARAMETER       :: PI = 3.14159

        NULL = CHAR(0)

C       GET IMAGE FILE NAMES
4004    WRITE(NOUT,*) 'FOR IMAGE FILE NAMES:'
        NPROJ = 200
        CALL FILSEQP(FILPAT,NLET,NUMBER,200,NPROJ, 
     &     'FILE PREFIX OR TEMPLATE (EG. PIC****)',IRTFLG)
     
        CALL FILGET(FILPAT,FILNAM,NLETF,1,IRTFLG)
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'O',IFORM,NSAM,NROW,NDUM,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) IERR=4
	IF (IFORM .NE .-1) IERR=2
	IF (IERR  .NE. 0)  GOTO 999

C       FIRST FILE ASSUMED TO BE THE AXIAL PROJ. IT IS ONLY OPENED
C       TO KEEP DIMENSIONS FOR CHECKING AGAINST DIMS. OF SUBSEQUENT FILES.

	NSAM1=NSAM
	NROW1=NROW
	CLOSE(LUN)

	DO IPROJ=2,NPROJ

          CALL FILGET(FILPAT,FILNAM,NLETF,IPROJ,IRTFLG)

          MAXIM = 0
          CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'O',IFORM,NSAM,NROW,NDUM,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
          IF (IRTFLG .NE. 0)    IERR=4
          IF (IFORM .NE. -1)    IERR=2
          IF (NSAM  .NE. NSAM1) IERR=1
          IF (NROW  .NE. NROW1) IERR=1
          IF (IERR  .NE. 0) GOTO 999

           CALL RDPRA('AXIAL VIEW: UNIT VECT. AX,AY, BX,BY',
     &         4,0,.FALSE.,FWA,NVAL,IRTFLG)

           ARX = FWA(1)
           ARY = FWA(2)
           BRX = FWA(3)
           BRY = FWA(4)

           WRITE(NOUT,33)ARX,ARY,BRX,BRY
33         FORMAT(2X,4F8.2)

           CALL RDPRA('TILTD VIEW: UNIT VECT. AX,AY, BX,BY',
     &         4,0,.FALSE.,FWA,NVAL,IRTFLG)

           ATX = FWA(1)
           ATY = FWA(2)
           BTX = FWA(3)
           BTY = FWA(4)

          CALL TILTFD(ARX,ARY,BRX,BRY,ATX,ATY,BTX,BTY,GAMR,GAMT,
     &               PHI,THA,GAMT2,IER)
          IF (IER .NE. 0) THEN
            IERR = 35
            GOTO 999
          ENDIF

45        PHID=PHI*180./PI
          THAD=THA*180./PI
          GAMR=GAMR*180./PI
          GAMT=GAMT*180./PI
          GAMT2=GAMT2*180./PI

          WRITE(NOUT,50)IPROJ,PHID,THAD,GAMR,GAMT,GAMT2
50        FORMAT(' PROJ. NO.',I3,'PHI=',F8.2,'  THETA=',F8.2,/
     &      ' ANGLE BETW. UNIT VECTORS IN AXIAL PROJ.',F8.2,/
     &      ' ANGLE BETW. UNIT VECTORS IN TILTED PROJ.',F8.2,/
     &      ' ANGLE CHECK',F8.2)

          IRTFLG = -1      
          VALUES(1) = 1.0  
          VALUES(2) = PHI
          VALUES(3) = THA
          CALL SETLAB(LUN,NSAM,FDUM,14,3,VALUES,'F',IRTFLG)
          VALUES(1) = 1.0
          CALL SETLAB(LUN,NSAM,FDUM,21,1,VALUES,'U',IRTFLG)

          CLOSE(LUN)
	ENDDO

	RETURN


999	CALL ERRT(IERR,'TILT  ',NE)
	RETURN
	END

