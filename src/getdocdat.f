C++*********************************************************************
C
C  GETDOCDAT.F   CREATED MAR 98 ArDean Leith
C                LUNDOCREDALL PARAMETERS          DEC   00 ARDEAN LEITH
C                LUNDOCREDSEQ RETURNS MAXY        APR 2003 ARDEAN LEITH
C                INCORE OPENDOC                   JUL 2003 ARDEAN LEITH
C                WANTSEQ BUG                      SEP 2003 ARDEAN LEITH
C                EMPTY DOC FILE BUG               DEC 2010 ARDEAN LEITH
C                (MYPID .LE. 0) CLOSE(LUNDOCFT)   MAR 2011 ARDEAN LEITH
C                (MYPID .LE. 0) CLOSE(LUNDOCFT)   MAR 2011 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C   GETDOCDAT(PROMPT,ASKNAME,DOCNAM,LUNDOC,GETSIZE,MAXXT,MAXYT,
C             RPOINTER,IRTFLG)
C
C   PARAMETERS:
C       PROMPT      PROMPT FOR DOCFILE NAME                       (SENT)
C       ASKNAME     LOGICAL VARIABLE TO ASK FOR DOCNAM            (SENT)
C       DOCNAM      DOC. FILE NAME                           (SENT/RET.)
C       LUNDOC      I/0 UNIT                                      (SENT)
C                   (< 0 IS FLAG FOR GETTING SEQUENTIAL DATA)
C       GETSIZE     LOGICAL VARIABLE TO ASK FOR ARRAY SIZE        (SENT)
C       MAXXT       COL. IN ARRAY RPOINTER                   (SENT/RET.)
C       MAXYT       ROWS IN ARRAY RPOINTER                   (SENT/RET.)
C                   IF GETSIZE IS TRUE, MAXXT & MAXYT SHOULD BE 
C                   SET TO ZERO ON ENTRY, OTHER-WISE THEY INDICATE
C                   MAXIMUM VALUES FOR RPOINTER ARRAY CREATION!!
C       RPOINTER    POINTER TO ARRAY (ALLOCATED HERE)             (RET.)
C       IRTFLG      ERROR FLAG (0 IS NORMAL)                      (RET.)
C
C  NOTES:   THIS PROGRAM ALWAYS ALLOCATES MEMORY OF ARRAY RPOINTER
C           EACH ROW IN RPOINTER HAS THE NUMBER OF THE REGISTERS IN THAT
C           KEY IN THE FIRST COLUMN (ZERO IF KEY NOT PRESENT) FOLLOWED 
C           BY THE REGISTER CONTENTS. (THIS MAY BE DIFFERENT IF GETTING
C           SEQUENTIAL DATA??? al)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

#ifdef NEVER 
        INTERFACE ! -------------UPDATE F90ALLOC.INC FOR PARAMETERS!!
        SUBROUTINE GETDOCDAT(PROMPT,ASKNAME,DOCNAM,LUNDOC,GETSIZE,
     &                       MAXXT, MAXYT,RPOINTER,IRTFLG)
        CHARACTER *(*), INTENT(IN) ::                 PROMPT
        LOGICAL, INTENT(IN) ::                        ASKNAME
        CHARACTER *(*), INTENT(INOUT) ::              DOCNAM
        INTEGER, INTENT(IN) ::                        LUNDOC
        LOGICAL, INTENT(IN) ::                        GETSIZE
        INTEGER, INTENT(INOUT) ::                     MAXXT
        INTEGER, INTENT(INOUT) ::                     MAXYT
        REAL, DIMENSION(:), POINTER  ::               RPOINTER
        INTEGER, INTENT(OUT) ::                       IRTFLG
        END SUBROUTINE GETDOCDAT
        END INTERFACE !--------------------
#endif

        SUBROUTINE GETDOCDAT(PROMPT,ASKNAME,DOCNAM,LUNDOC,GETSIZE,
     &                       MAXXT, MAXYT,RPOINTER,IRTFLG)

C       F90 SPECIFIC CODE
C       UPDATE F90ALLOC.INC FOR PARAMETERS!!

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC' 

        CHARACTER *(*), INTENT(IN)    :: PROMPT
        LOGICAL, INTENT(IN)           :: ASKNAME
        CHARACTER *(*), INTENT(INOUT) :: DOCNAM
        INTEGER, INTENT(IN)           :: LUNDOC
        LOGICAL, INTENT(IN)           :: GETSIZE
        INTEGER, INTENT(INOUT)        :: MAXXT
        INTEGER, INTENT(INOUT)        :: MAXYT
        REAL, POINTER                 :: RPOINTER(:,:)
        INTEGER, INTENT(OUT)          :: IRTFLG
        LOGICAL                       :: WANTSEQ,NEWFILE

        INTEGER,PARAMETER             :: NMAX = 10

        CHARACTER (LEN=180)           :: PROMPTT ! > THAN PROMPT!!

        LOGICAL,PARAMETER             :: ISOLDFILE = .TRUE.
        LOGICAL,PARAMETER             :: APPEND    = .FALSE.
        LOGICAL,PARAMETER             :: MESSAGE   = .FALSE.
        LOGICAL,PARAMETER             :: ASKNAME2  = .FALSE.
        LOGICAL                       :: ADDEXT    = .FALSE.

        INTEGER                       :: IWANTY,MAXX,MAXY
        INTEGER                       :: ICOMM,MYPID,MPIERR
        INTEGER                       :: LUNDOCFT,LENT,NLETD,NF,NLET
        INTEGER                       :: LUNDOCF,KEYUSED,NGOT

        INTEGER                       :: lnblnkn

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYMPID

        IRTFLG  = 1

C       FLAG TO RECOVER MULTIPLE VALUES HAVEING SAME KEY ('DOC RE')
        WANTSEQ  = (LUNDOC < 0)
        LUNDOCFT = ABS(LUNDOC)

        IF (ASKNAME) THEN
C          ~9 ON PROMPT ALLOWS NON-SPIDER EXTENSION TO BE INPUT
           LENT    = lnblnkn(PROMPT)
           PROMPTT = PROMPT(1:LENT) // '~9'

C          INPUT DOCUMENT FILENAME, APPEND DATEXC
           CALL FILERD(DOCNAM,NLETD,DATEXC,PROMPTT,NF)
           IF (NF .NE. 0) RETURN
        ENDIF

        !write(6,*) ' docnam: ',nletd,docnam(1:30)

C       INPUT DOCUMENT FILENAME,
        CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOCFT,LUNDOCF,ASKNAME2,
     &             ' ', ISOLDFILE,APPEND,MESSAGE,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       ECHO FIRST HEADER FROM FILE
        CALL LUNDOCSAYHDR(LUNDOCF,NOUT,IRTFLG)

        IF (GETSIZE) THEN
C          FIND MAXY & MAXY BY READING DOC FILE


           CALL LUNDOCINFO(LUNDOCF,MAXY,MAXX,KEYUSED,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CLOSE(LUNDOCFT)
              RETURN
           ENDIF
           IF (LUNDOCF > 0) REWIND(LUNDOCF)

           IF (WANTSEQ) THEN
C             WANT JUST DATA ARRAY
              IWANTY = MAXY
              MAXY   = KEYUSED
           ELSE
C             WANT KEYED ARRAY
              MAXX = MAXX + 1
           ENDIF
           IF (MAXXT > 0) MAXX = MIN(MAXX,MAXXT)
           IF (MAXYT > 0) MAXY = MIN(MAXY,MAXYT)
        ELSE
C          USE MAXX & MAXY SENT BY CALLING PROGRAM
           MAXX = MAXXT
           MAXY = MAXYT
        ENDIF

        IF (MAXX <= 0 .OR. MAXY <= 0) THEN
           WRITE(NOUT,*) ' *** WARNING EMPTY DOC FILE '
           CLOSE(LUNDOCFT)
           IRTFLG = 0
           MAXXT  = 0
           MAXYT  = 0
           RETURN
        ENDIF

        ALLOCATE(RPOINTER(MAXX,MAXY),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN  
           CALL  ERRT(46,'GETDOCDAT; RPOINTER',MAXX*MAXY)
           RETURN
        ENDIF
 
C       RETRIEVE DATA
        IF (WANTSEQ) THEN
C          TO RECOVER MULTIPLE VALUES HAVING SAME KEY
           CALL LUNDOCREDSEQ(LUNDOCF,RPOINTER,MAXX,MAXY,IWANTY,MAXYT,
     &                       IRTFLG)
        ELSE
           CALL LUNDOCREDALL(LUNDOCF,RPOINTER,MAXX,MAXY,.TRUE.,
     &                       NGOT,IRTFLG)
           MAXYT = MAXY
        ENDIF

        MAXXT  = MAXX
        IF (MYPID  <= 0)   CLOSE(LUNDOCFT)
        IF (IRTFLG .NE. 0) RETURN

        IRTFLG = 0

        END
