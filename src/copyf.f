
C++*********************************************************************
C
C    COPYF.FOR -- CREATED JULY 17 1989 ardean leith                     
C                 USED OPAUXFILE                  FEB 99 ARDEAN LEITH
C                 OPFILEC                         FEB 03 ARDEAN LEITH
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
C    COPYF(LUN1,LUN2)
C
C    PURPOSE:   COPIES AN EDITABLE TXT FILE INTO A SPIDER IMAGE FILE
C
C    PARAMETERS:   LUN1      INPUT FILE UNIT NUMBER
C                  LUN2      OUTPUT FILE UNIT NUMBER
C--*********************************************************************

	SUBROUTINE COPYF(LUN1,LUN2)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUF(NBUFSIZ)

        CHARACTER(LEN=MAXNAM) :: FILNAM,FMT

        CHARACTER * 1         :: ANSW
        CHARACTER * 1         :: NULL = CHAR(0)
        LOGICAL               :: FREEFMT

C       OPEN INPUT FILE AS SEQUENTIAL ACCESS, FORMATTED
10      LENREC = 0
        CALL OPAUXFILE(.TRUE.,FILNAM,DATEXC,LUN1,LENREC,'O',
     &                       'EDITABLE IMAGE INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

12      CALL RDPRMC(ANSW,NC,.TRUE.,
     &     'ARE NX, NY, & NZ IN FIRST LINE OF FILE? (Y/N)',
     &     NULL,IRTFLG)
        IF (IRTFLG == -1)  GOTO 10

        IF (NC == 0 .OR. ANSW .NE. 'N') THEN
C          CAN GET NX OUT OF FILE
           READ(LUN1,*,IOSTAT=IOS) NX,NY,NZ
           IF (IOS .NE. 0) THEN
             CALL ERRT(101,'*** ERROR READING FILE',NDUM)
             GOTO 9999
           ENDIF

        ELSE
C          ASK USER FOR NX, ETC.
           CALL RDPRI3S(NX,NY,NZ,NOT_USED,
     &                'NX, NY & NZ',IRTFLG)
           IF (IRTFLG == -1) GOTO 12
        ENDIF

        IF (NZ .LE. 0) NZ = 1

16      CALL RDPRMC(FMT,NC,.TRUE.,
     &     'FORMAT DESCRIPTION (OR <CR> FOR FREE FORMAT)',
     &     NULL,IRTFLG)
        IF (IRTFLG == -1) GOTO 12

C       DEFAULT IS FREE FORMAT
        FREEFMT = (NC == 0)

20      IFORM = 1
        IF (NZ > 1) IFORM = 3
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'U',IFORM,NX,NY,NZ,
     &             MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG  == -1) GOTO 16
        IF (IRTFLG .NE. 0) GOTO 9999
        
	DO I=1,NY * NZ
C          IOSTAT NEEDED FOR PARTIAL LINES

           IF (FREEFMT) THEN
              READ(LUN1,*,        IOSTAT=IOS) (BUF(J),J=1,NX)
           ELSE
              READ(LUN1,FMT(1:NC),IOSTAT=IOS) (BUF(J),J=1,NX)
           ENDIF

           IF (IOS .NE. 0) THEN
              write(6,*) 'FMT:',fmt(1:10),':'
              write(6,*) 'Read error, iostat:',ios
           endif

           !write(6,*) (BUF(J),J=1,nx)

           CALL WRTLIN(LUN2,BUF,NX,I)
        ENDDO

9999	CLOSE(LUN1)
	CLOSE(LUN2)

	RETURN

        END


