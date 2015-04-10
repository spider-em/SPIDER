
C++*********************************************************************
C
C SAVDOCQ.F           REWRITTEN                   SEPT 96  ArDean Leith
C                     USED LUNDOC                 JUNE 99  ArDean Leith
C                     OPENDOC PARAMETERS          DEC 2000 ARDEAN LEITH
C                     FLUSH                       SEP 2002 ArDean Leith
C                     INCORE OPENDOC              JUL 2003 ARDEAN LEITH
C                     LUNDOCWRTDAT FLUSHES        OCT 2003 ARDEAN LEITH
C                     REMOVED / FROM COMMENT      JUN 2008 ARDEAN LEITH
C                     USED: LUNDOCPUTCOM          NOV 2009 ARDEAN LEITH
C                     HEADER                      FEB 2010 ARDEAN LEITH
C                     SUPPORT FOR //COMMENT       AUG 2010 ARDEAN LEITH
C                     LONGER COMMENT KEYS         AUG 2010 ARDEAN LEITH
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
C  SAVDOCQ(DOCNAM,NLET,HEADER,IRTFLG)
C
C  PURPOSE:    SUBROUTINE TO SAVE PARTICULAR REGISTERS IN A 
C              FILE-BASED DOCUMENT FILE.  CALLED FROM COMMAND LINE.
C
C  OPERATION USAGE:
C             SD [KEY],[REG1],[REG2],[REG3]
C             SD E               (END DOC FILE USAGE, CLOSE FILE)
C             SD /NEW COMMENT FOR DOC. FILE 
C             SD //NEW COMMENT FOR DOC. FILE 
C                 
C  PARAMETERS:
C      DOCNAME   FILE NAME                            (SENT)
C      NLET      NUMBER OF CHAR IN DOCNAME            (SENT)
C      HEADER   PUT HEADER IN FILE ON OPEN            (SENT)
C      IRTFLG    ERROR FLAG                           (RETURNED)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE SAVDOCQ(DOCNAM,NLET,HEADER,IRTFLG)

        INCLUDE 'CMLIMIT.INC' 
        INCLUDE 'CMBLOCK.INC' 

        CHARACTER(LEN=*)      :: DOCNAM

C	MAX. NUMBER OF SIMULTANEOUS DOCUMENT FILES ALLOWED
	INTEGER, PARAMETER    :: MAXICDOCS = 20 
        CHARACTER(LEN=MAXNAM) :: DCFILE(MAXICDOCS)

C       MAXIMUM NUMBER OF REGISTERS SAVED
        INTEGER, PARAMETER    :: MAXLIST=50
	REAL                  :: DLIST(MAXLIST)
	INTEGER               :: ILIST(MAXLIST)
	LOGICAL               :: NEWFILE,COMOUT,HEADER,SDT
        LOGICAL               :: GETANS,UPPER,WANTSUB,SAYPRMT,SAYANS
        LOGICAL               :: ENDATSEMI,STRIP
        CHARACTER(LEN=160)    :: COMSTR
  
C       NEEDED FOR FUTURE CALLS
        SAVE           DCFILE, ICLAST

        DATA           DCFILE/MAXICDOCS*'*'/
        DATA           ICLAST/0/

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       SET ERROR RETURN
        IRTFLG = 1

        IF (FCHAR(4:4) .EQ. 'E') THEN
C          WANT TO CEASE USING A DOC FILE -----------------------------
C          DETERMINE WHICH OF THE OLD NAMES NEEDS TO BE CLOSED
           DO  IC=1,MAXICDOCS
             IF  (DOCNAM .EQ. DCFILE(IC)) THEN
C              CHANGE DCFILE SO THAT IT CANNOT BE MATCHED WITH ANY FILE NAME.
               DCFILE(IC)  = '*' 
               CLOSE(200+IC)
               IRTFLG     = 0
               RETURN
             ENDIF
           ENDDO

           IF (MYPID .LE. 0)WRITE(NOUT,*) ' DOCUMENT FILE NOT OPEN NOW '
           RETURN
        ENDIF


C       COMPARE NEW DOCUMENT FILE WITH OLD NAMES ---------------------

        IC     = 0

C       NAME IS MOST-LIKELY STILL THE SAME
        IF (ICLAST .GT. 0 .AND. ICLAST .LE. MAXICDOCS .AND.
     &      DOCNAM .EQ. DCFILE(ICLAST)) THEN
C          DOCNAM IS ALREADY IN-CORE
           IC    = ICLAST
           NDOC  = 200 + IC
           GOTO 14
        ENDIF

C       FILL FROM FRONT TO BACK
        DO ICT = 1,MAXICDOCS
           IF (DOCNAM .EQ. DCFILE(ICT)) THEN
C             DOCNAM IS ALREADY OPEN
              IC   = ICT
              NDOC = 200 + IC
              GOTO 14
           ENDIF

C          REMEMBER FIRST EMPTY LOCATION
           IF (IC .EQ. 0 .AND. DCFILE(ICT)  .EQ. '*') IC = ICT
        ENDDO

C       DOCUMENT NAME NOT FOUND IN LIST. 
        IF (IC .EQ. 0) THEN
C          LIST FULL, CLOSE SOMETHING & OPEN THIS
           IC = ICLAST + 1
           IF (IC .GT. MAXICDOCS) IC = 1
           DCFILE(IC) = '*'
           CLOSE(200 + IC)
        ENDIF

C       OPEN DOC FILE ON UNIT IC + 200 
        NDOC = IC + 200 
      
C       FILE MAY NOT HAVE BEEN ACCESSED YET TODAY, PUT NEW HEADER IN? 
        CALL OPENDOC(DOCNAM,.FALSE.,NLET,NDOC,NIC,.FALSE.,' ',
     &               .FALSE.,.TRUE.,HEADER,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NIC .LE. 0) THEN
           CALL ERRT(101,'USE <SD IC> FOR INCORE DOC FILE.',NE)
           RETURN
        ENDIF 

C       PUT NAME IN LIST OF OPENED FILES
        DCFILE(IC)= DOCNAM
        ICLAST = IC

C       PARSE REGISTER LINE, CHECK FOR IKEY & NLIST ------------------
14      CALL REG_DOC_PARSE(FCHAR(4:),COMOUT,IKEY,NLIST,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (COMOUT) THEN
C          JUST WANT TO PUT A TEXT COMMENT IN THE DOC FILE.

           WANTSUB = (FCHAR(5:5) .EQ. '/')
           IF (WANTSUB) THEN
C             SUBSTITUTE VARIABLES INTO COMMENT

C             DO NOT UPPERCASE THE COMMENT, STRIP AFTER ;
              UPPER     = .FALSE.
              SAYPRMT   = .FALSE.
              SAYANS    = .FALSE.
              ENDATSEMI = .TRUE.
              GETANS    = .FALSE.
              STRIP     = .TRUE.
              CALL RDPR(FCHAR(6:),NC,COMSTR,GETANS,
     &            UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              CALL LUNDOCPUTCOM(NDOC,COMSTR(1:NC),IRTFLG)
           ELSE
C             DO NOT SUBSTITUTE VARIABLES INTO COMMENT
              CALL LUNDOCPUTCOM(NDOC,FCHAR(5:),IRTFLG)
           ENDIF
           RETURN
        ENDIF


        IF (MAXLIST .LT. NLIST) THEN
           CALL ERRT(102,'MAX. NUMBER OF REGISTERS',MAXLIST)
           NLIST = MAXLIST
        ENDIF

C       RETRIEVE DATA FROM THE REGISTER(S) LISTED IN NSEL INTO: DLIST
        CALL REG_GET_NSELA(NLIST,DLIST,IRTFLG)

C       WRITE DATA FROM DLIST INTO DOC FILE FOR THIS KEY
        SDT =  (FCHAR(4:4) .EQ. 'T') 
        CALL LUNDOCWRTDAT(NDOC,IKEY,DLIST,NLIST,SDT,IRTFLG)

C       LEAVE FILE OPEN, FOR NEXT USE
        IRTFLG = 0

	RETURN
	END

