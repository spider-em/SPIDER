
C++*********************************************************************
C
C UDOC.F                      LONG FILE NAMES     FEB   89 ArDean Leith
C                             CHAR VARIABLES      AUG   89 ArDean Leith
C                             DOC FILE LEFT OPEN  NOV   89 ArDean Leith
C                             USED LUNDOC         JUN   99 ArDean Leith
C                             TILLEND BUG         NOV   99 ArDean Leith
C                             ADDED NEEDREWIND    JUN   00 ArDean Leith
C                             OPENDOC PARAMETERS  DEC   00 ARDEAN LEITH
C                             ICOUNT > NLIST BUG  AUG   02 ARDEAN LEITH
C                             ERRT PARAM          OCT   10 ARDEAN LEITH
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
C  UDOC(CFUNC,NDOC)
C
C  PURPOSE:  SUPPORTS OPERATION TO RETRIEVE A LINE OF REGISTERS FROM 
C            DOC FILE. ALWAYS SOLICITS FILENAME,  OPENS FILE IF NAME
C            NOT SAME AS PREVIOUS DOC FILE USED BY THIS OPERATION.
C            REGISTERS ARE SPECIFED ON COMMAND LINE.
C            JUNE 2000 IF "UD S" NOW WILL REWIND AND TRY AGAIN IF IT
C            CAN NOT FIND KEY ON FIRST PASS THRU REMAINING PART OF FILE.
C             
C  TYPICAL USAGE: 
C            UD 11,X10,X11       (RETRIEVE)
C            UD X11,X10,X11      (RETRIEVE)
C            UD S,11,X10,X11     (SEQUENTIAL RETRIEVE)
C            UD E                (CLOSES LAST DOC. FILE OPENED WITH UD)
C            UD -5,X5,X10,X15    (RETRIEVE COMMENT KEY LINES (;KEY))
C
C  PARAMETERS:    CFUNC   OPERATION READ INTO DRIVER           (SENT)
C                 NDOC    LUN NUMBER OF FILE                   (SENT)
C
C  CALLS          UNSAV
C
C--*********************************************************************

	SUBROUTINE UDOC(CFUNC,NDOC)

C       SAVE IS NEEDED FEB 99 al
        SAVE

        INCLUDE 'CMLIMIT.INC' 
        INCLUDE 'CMBLOCK.INC'

        CHARACTER(LEN=*)            :: CFUNC
        INTEGER                     :: NDOC

        INTEGER, PARAMETER          :: MAXLIST=100
	REAL                        :: DLIST(MAXLIST)
                 
        CHARACTER(LEN=MAXNAM)       :: DOCNAM,DOC
        CHARACTER(LEN=MAXNAM), SAVE :: OLDNAM
        INTEGER, SAVE               :: LUNOLD

        CHARACTER(LEN=1)            :: NULL
	LOGICAL                     :: NEWFILE,COMOUT,TILLEND,NEEDREWIND

	DATA           LUNOLD/0/
        DATA           OLDNAM/'-'/

        NULL        = CHAR(0)

        IF (CFUNC(4:4) .EQ. 'E') THEN
C          ENDS USE OF THIS DOC FILE
           CLOSE(NDOC)
           OLDNAM = NULL
           LUNOLD = 0
           RETURN
        ENDIF

C       GET DOC FILE NAME
	CALL FILERD(DOC,NLET,NULL,'DOCUMENT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       DEFAULT TO NEED REWIND IF KEY NOT FOUND
        NEEDREWIND = .TRUE.

        IF (NDOC .NE. LUNOLD .OR. DOC .NE. OLDNAM) THEN
C         DIFFERENT DOC FILE OR DIFFERENT LUN FOR DOC FILE
          IF (LUNOLD .NE. 0) CLOSE(LUNOLD)
          LUNOLD = 0
          OLDNAM = NULL

C         OPEN THE DOC FILE
          DOCNAM = DOC
          CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOC,NIC,.FALSE.,' ',
     &                 .TRUE.,.TRUE.,.TRUE.,NEWFILE,IER)
          IF (IER .NE. 0) RETURN

C         ECHO FIRST COMMENT FROM DOC FILE UPON OPENING
          CALL LUNDOCSAYHDR(NDOC,NOUT,IRTFLG)

          OLDNAM     = DOC
          LUNOLD     = NDOC
C         NO NEED TO REWIND IF KEY NOT FOUND
          NEEDREWIND = .FALSE.
        ENDIF

        IGO     = 4
        TILLEND = .TRUE.
	IF (CFUNC(4:4) .EQ. 'S') THEN
C          USE SEQUENTIAL SEARCH MODE
           TILLEND = .FALSE.
           IGO     = 6
        ENDIF

C       REGISTER LINE ALREADY LOADED IN RDPR 
C       PARSE REGISTER LINE TO GET IKEY & NLIST 
        CALL REG_DOC_PARSE(CFUNC(IGO:),COMOUT,IKEY,NLIST,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       RETRIEVE THE LIST OF VALUES FOR THIS KEY, LOOK TILL EOF
C       IF TILLEND IS TRUE

10      CALL LUNDOCREDDAT(NDOC,IKEY,DLIST,NLIST,ICOUNT,
     &                    TILLEND,.FALSE.,IRTFLG)

        IF (IRTFLG .EQ. 0 .AND. ICOUNT .GT. 0) THEN
C         SUCCESSFUL RECOVERY OF KEY,  
          ICOUNT = MIN(ICOUNT,NLIST)
          CALL REG_SET_NSELA(ICOUNT,DLIST,.TRUE.,IRTFLG)

        ELSE IF (CFUNC(4:4) .EQ. 'S' .AND. NEEDREWIND) THEN
           REWIND(NDOC)
C          NO NEED TO REWIND IF KEY NOT FOUND ON SECOND PASS
           NEEDREWIND = .FALSE.
C          TRY SECOND PASS THRU FILE
           GOTO 10

        ELSE
           CALL ERRT(102,'KEY NOT FOUND',IKEY)
        ENDIF

        IF (TILLEND) THEN
           REWIND(NDOC)
        ENDIF

C       DO NOT CLOSE FILE UNTIL 'UD E' IS GIVEN!

	RETURN
	END


