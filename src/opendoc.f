
C++*********************************************************************
C
C OPENDOC.F           CHANGED PARAMETERS         DEC 2000 ARDEAN LEITH
C                     NO LONGER RETURN EXTENSION JAN 2001 ARDEAN LEITH
C                     LUNRET ADDED               JUL 2003 ARDEAN LEITH
C                     MPI                        SEP 2003 CHAO YANG
C                     CSTRING TOO SHORT          JUL 2006 ARDEAN LEITH
C                     'UD MAX' SUPPORT           MAY 2009 ARDEAN LEITH
C                     TIME ERROR IF NO MSG       JUN 2009 ARDEAN LEITH
C                     NO INCORE IF IRTFLG = -8   APR 2013 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOC,LUNRET,GETNAME,PROMPT,
C        ISOLDFILE,APPEND,MESSAGE,NEWFILE,IRTFLG)
C
C PURPOSE:     OPEN DOCUMENT FILE 
C
C OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOC,LUNRET,GETNAME,PROMPT,
C        ISOLDFILE,APPEND,MESSAGE,NEWFILE,IRTFLG)
C
C PARAMETERS:     DOCNAM        NAME OF DOC FILE                SENT/RET
C                 ADDEXT        LOGICAL FLAG TO ADD .EXT        SENT
C                 NLET          NO. OF CHAR. IN DOCNAM (NO EXT) RETURNED
C                 LUNDOC        IO UNIT FOR DOC FILE            SENT
C                 LUNRET        LUN OR INDEX FOR INCOREDOC FILE SENT/RET.
C                               (< 0 INDICATES INCORE FILE)
C                 GETNAME       FLAG TO REQUEST NAME            SENT
C                 PROMPT        PROMPT TO REQUEST NAME          SENT
C                 ISOLDFILE     FLAG THAT FILE IS EXISTING      SENT
C                 APPEND        FLAG TO OPEN FILE AS APPEND     SENT
C                 MESSAGE       FLAG TO WRITE COMMENT           SENT
C                 NEWFILE       FLAG THAT FILE WAS NEW          RETURNED
C	          IRTFLG        ERROR FLAG (0 IS NORMAL)        RETURNED
C                               (-9 ON INPUT IS DO NOT SAYOPEN) 
C                               (-8 ON INPUT IS DO NOT USE EXISTING
C                                INCORE FILE) 
C
C--*******************************************************************
		
        SUBROUTINE OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOC,LUNRET,GETNAME,
     &                  PROMPT,ISOLDFILE,APPEND,MESSAGE,NEWFILE,IRTFLG)

        IMPLICIT NONE

	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 
 
	CHARACTER (LEN=*)      :: DOCNAM
	LOGICAL                :: ADDEXT
	INTEGER                :: NLET,LUNDOC,LUNRET
	LOGICAL                :: GETNAME
	CHARACTER (LEN=*)      :: PROMPT
	LOGICAL                :: ISOLDFILE
	LOGICAL                :: APPEND
	LOGICAL                :: MESSAGE
        LOGICAL                :: NEWFILE
	INTEGER                :: IRTFLG

	CHARACTER (LEN=MAXNAM) :: DOCNAMPE
	CHARACTER (LEN=160)    :: CSTRING
        CHARACTER (LEN=12)     :: CDATT
        CHARACTER (LEN=8)      :: CTIMT
        CHARACTER (LEN=1)      :: NULL = CHAR(0)
	LOGICAL                :: EX
        LOGICAL                :: SAYOPEN,ISOPEN

        LOGICAL                :: EXTOK,GOTEXT,ICOK
        INTEGER                :: ICOMM,MYPID,MPIERR,LENP,NLETPE,MT,NE
        INTEGER                :: LUNOP,IOS
 
        INTEGER                :: lnblnkn

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYMPID

        SAYOPEN = (IRTFLG .NE. -9)
        ICOK    = (IRTFLG .NE. -8)

        IF (GETNAME) THEN
           CALL FILERD(DOCNAM,NLET,NULL,PROMPT,IRTFLG)
           IF (IRTFLG == -1) RETURN
        ELSE
           NLET = lnblnkn(DOCNAM)
        ENDIF
        IRTFLG = 2

        IF (ADDEXT) THEN
C          MERGE DOCNAM WITH DATEXC
           CALL FILNAMANDEXT(DOCNAM,DATEXC,DOCNAMPE,NLETPE,.TRUE.,
     &                       IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
           DOCNAMPE = DOCNAM
           NLETPE   = NLET
        ENDIF

C       SEE IF THIS FILE IS ALREADY OPEN IN-CORE
        CALL ISDOCINCORE(DOCNAM(1:NLET),LUNRET,MT,IRTFLG)
  
        IF (FCHAR(1:2) .EQ. 'UD' .AND. 
     &     (FCHAR(4:5) .NE. 'IC' .AND. 
     &      FCHAR(4:4) .NE. 'N'  .AND. 
     &      FCHAR(4:4) .NE. 'M')) THEN
C          (IF 'UD ??'  BUT NOT 'UD IC'..., DONT USE INCORE EVEN IF PRESENT)
           CONTINUE

        ELSEIF (LUNRET > 0 .AND. ICOK) THEN
C          IN-CORE FILE EXISTS, USE IT INSTEAD OF PHYSICAL FILE
           LUNRET = -LUNRET
           IRTFLG = 0
           RETURN
        ENDIF
        LUNRET = LUNDOC

C       SEE IF PHYSICAL FILE EXISTS
        IF (MYPID <= 0) THEN
           INQUIRE(FILE=DOCNAMPE(1:NLETPE),EXIST=EX,OPENED=ISOPEN,
     &             NUMBER=LUNOP,IOSTAT=IOS)
        ENDIF
#ifdef USE_MPI
c       write(6,*) ' opendoc; bcast ios: ',ios,mypid
        CALL BCAST_MPI('OPENDOC','IOS', IOS,1, 'I',ICOMM)
c       write(6,*) ' opendoc; bcast lunop: ',lunop,mypid
        CALL BCAST_MPI('OPENDOC','LUNOP', LUNOP,1, 'I',ICOMM)
c       write(6,*) ' opendoc; bcast isopen: ',isopen,mypid
        CALL BCAST_MPI('OPENDOC','ISOPEN', ISOPEN,1, 'L',ICOMM)
c       write(6,*) ' opendoc; bcast ex: ',ex,mypid
        CALL BCAST_MPI('OPENDOC','EX', EX,1, 'L',ICOMM)
#endif

        IF (IOS .NE. 0) THEN
	  WRITE(NOUT,*)' *** ERROR INQUIRING FILE: ',DOCNAMPE(1:NLETPE)
           CALL ERRT(100,' ',NE)
           IRTFLG = 1
           RETURN
        ENDIF

        NEWFILE = .NOT. EX

        IF (ISOLDFILE .AND. .NOT. EX) THEN
C          OLD DOC FILE SHOULD EXIST FOR READING, BUT DOES NOT!
	   WRITE(NOUT,*)' *** DOC FILE DOES NOT EXIST: ',
     &                  DOCNAMPE(1:NLETPE)
           CALL ERRT(100,' ',NE)
           IRTFLG = 1
           RETURN

        ELSEIF (ISOLDFILE .AND. EX) THEN
C          SEE IF FILE IS ALREADY OPEN
           IF (ISOPEN) THEN
              IF (MYPID <= 0) THEN
                 CLOSE(LUNOP)
	         WRITE(NOUT,*) ' FILE ALREADY OPEN, HAS BEEN CLOSED!'
              ENDIF
           ENDIF

C          OPEN EXISTING DOC FILE FOR READING
           IF (MYPID .LE. 0) THEN
              OPEN(UNIT=LUNDOC,FILE=DOCNAMPE(1:NLETPE),STATUS='OLD',
     &             IOSTAT=IOS)
           ENDIF
 
        ELSEIF (.NOT. ISOLDFILE .AND. EX) THEN
C          OPEN EXISTING EXISTING DOC FILE FOR WRITING
           IF (APPEND) THEN
              IF (MYPID .LE. 0) THEN
                 OPEN(UNIT=LUNDOC,FILE=DOCNAMPE(1:NLETPE),STATUS='OLD',
     &               ACCESS="SEQUENTIAL",POSITION="APPEND",IOSTAT=IOS)
              ENDIF
          ELSE
              IF (MYPID <= 0) THEN
                 OPEN(UNIT=LUNDOC,FILE=DOCNAMPE(1:NLETPE),STATUS='OLD',
     &               ACCESS='SEQUENTIAL',IOSTAT=IOS)
              ENDIF
           ENDIF

        ELSEIF (.NOT. ISOLDFILE) THEN
C          OPEN NEW DOC FILE FOR WRITING 
           IF (MYPID <= 0) THEN
              OPEN(UNIT=LUNDOC,FILE=DOCNAMPE(1:NLETPE),STATUS='UNKNOWN',
     &             IOSTAT=IOS)
           ENDIF
        ENDIF

        IF (IOS .NE. 0) THEN
	  WRITE(NOUT,*) ' *** ERROR OPENING DOC FILE: ',
     &                  DOCNAMPE(1:NLETPE)
          CALL ERRT(100,' ',NE)
          IRTFLG = 1
          RETURN
        ENDIF

        CALL DATE_2K(CDATT)
        CALL MYTIME(CTIMT)

        IF (NEWFILE .AND. MESSAGE) THEN
C          WRITE HEADER INTO FILE
           IF (MYPID .LE. 0) THEN
              WRITE(LUNDOC,90) PRJEXC(1:3),DATEXC(1:3),
     &                       CDATT(1:11),CTIMT,DOCNAMPE(1:NLETPE)
           ENDIF
 90        FORMAT(' ;' ,A,'/',A,3X,A,' AT ',A,3X,A)
        ENDIF

        IF (SAYOPEN .AND. MYPID .LE. 0) THEN
           IF (NEWFILE) THEN
              WRITE(NOUT,92) CDATT(1:11),CTIMT, DOCNAM(1:NLET)
 92           FORMAT('  ',A,' AT ',A,3X,' OPENED NEW DOC FILE: ',A)
              IF (USE_SPIRE) THEN
                 WRITE(CSTRING,92) CDATT(1:11),CTIMT, DOCNAM(1:NLET)
                 CALL SPIREOUT(CSTRING,IRTFLG)
              ENDIF
           ELSE
              WRITE(NOUT,*) ' OPENED EXISTING DOC FILE: ',DOCNAM(1:NLET)
           ENDIF
        ENDIF

        IRTFLG = 0

	END
