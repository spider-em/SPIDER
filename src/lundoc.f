
C++*********************************************************************
C                NEW                              APR 1999 ARDEAN LEITH
C                ADDED FILNAMSUB TO LUNDOCPARSE   JAN 2000 ARDEAN LEITH
C                NCOL (MAX NO OF COLUMNS) ADDED   JUN 2000 BILL BAXTER
C                REMOVED CONTINUATION LINES       JUL 2000 ARDEAN LEITH
C                LUNDOCREDALL PARAMETERS          DEC 2000 ARDEAN LEITH
C                CALL INDEXTOREG                  MAR 2001 ARDEAN LEITH
C                REMOVED FILNAMSUB IN LUNDOCPARSE APR 2001 ARDEAN LEITH
C                ADDED LUNDOCREDSLI               APR 2001 ARDEAN LEITH
C                INCREASED REGISTER NUMBERS       MAY 2001 ARDEAN LEITH
C                LUNDOCREDDAT GOBACK              SEP 2002 ARDEAN LEITH
C                LUNDOCREDSEQ RETURNS MAXY        APR 2003 ARDEAN LEITH
C                LUNDOCREDSEL NGOT BUG            MAY 2003 ARDEAN LEITH
C                INCORE SUPPORT                   JUL 2003 ARDEAN LEITH
C                K WRONG IN LUNDOCPARSE           AUG 2003 ARDEAN LEITH
C                NEW DOC FILE FORMAT              FEB 2004 ARDEAN LEITH
C                LUNDOCGETCOM DLIST OVERFLOW      APR 2004 ARDEAN LEITH
C                LUNDOCGETCOM ICOUNT BUG          APR 2004 ARDEAN LEITH
C                LUNDOCREDALLI ICOUNT BUG         JUN 2004 ARDEAN LEITH
C                LUNDOCREDALLI COL OVERFLOW BUG   JUL 2005 ARDEAN LEITH
C                LUNDOCSAYHDR MYPID BUG           SEP 2005 ARDEAN LEITH
C                [] DEFAULT FOR REGS.             NOV 2005 ARDEAN LEITH
C                LUNDOCREDNXT LOCDOC ARROW BUG    FEB 2007 ARDEAN LEITH
C                NEXTKEY() ADDED, MAXICDOC=12     FEB 2007 ARDEAN LEITH
C                TOLERATES BAD ; FROM JWEB        SEP 2007 ARDEAN LEITH
C                MPI BCAST CHANGES                OCT 2008 ARDEAN LEITH
C                NEXTKEY SET TO 1                 MAY 2009 ARDEAN LEITH
C                LUNDOCPUTCOM FORMAT              NOV 2009 ARDEAN LEITH
C                MPI_SET                          NOV 2009 ARDEAN LEITH
C                VAR. SUBSTITUED IN TEXT COMMENT  AUG 2010 ARDEAN LEITH
C                COMMENT KEY > 9999               AUG 2010 ARDEAN LEITH
C                LUNDOCGETCOM   END=998           OCT 2010 ARDEAN LEITH
C                ..GETCOM: READ (LUNDOC,81        OCT 2010 ARDEAN LEITH
C                ADDED MPI BARRIER                MAR 2011 ARDEAN LEITH
C                LUNDOCINFO ; / BUG               APR 2011 ARDEAN LEITH
C                LUNDOCREDLIN ; ERRT BUG          APR 2011 ARDEAN LEITH
C                FORMAT(  1PG13.6)                APR 2011 ARDEAN LEITH
C                LUNDOCWRTDATF FORMAT             APR 2012 ARDEAN LEITH
C                LABEL 92 UNDEFINED IF NOT MPI    MAR 2015 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C    LUNDOC
C
C    PURPOSE:    SUPPORT SUBROUTINES FOR DOCUMENT FILE HANDLING.
C     ------------------------- LUNDOCREDLIN -----------------------------
C     ------------------------- LUNDOCWRTDAT -----------------------------
C     ------------------------- LUNDOCWRTDATF ----------------------------
C     ------------------------- LUNDOCREDDAT -----------------------------
C     ------------------------- LUNDOCREDALL -----------------------------
C     ------------------------- LUNDOCREDALLI --------------------------
C     ------------------------- LUNDOCREDSLC -----------------------------
C     ------------------------- LUNDOCREDSEQ -----------------------------
C     ------------------------- LUNDOCREDSEL -----------------------------
C     ------------------------- LUNDOCREDANG -----------------------------
C     ------------------------- LUNDOCGETKEY -----------------------------
C     ------------------------- LUNDOCREDNXT -----------------------------
C     ------------------------- LUNDOCINFO -----------------------------
C     ------------------------- LUNDOCGETCOM -----------------------------
C     ------------------------- LUNDOCPUTCOM -----------------------------
C     ------------------------- LUNDOCSAYHDR -----------------------------
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

C*********************** DOCIC_INFO ******************************

        MODULE DOCIC_INFO

C          NEEDED FOR MAXNAM
           INCLUDE 'CMLIMIT.INC' 

C          ARRAY OF POINTERS TO DOCUMENT STORAGE ARRAYS
           INTEGER, PARAMETER :: MAXICDOCS = 12

           TYPE REAL_POINTER
              REAL, DIMENSION(:), POINTER :: IPT 
           END TYPE REAL_POINTER

           TYPE(REAL_POINTER), DIMENSION(MAXICDOCS) :: LOCDOC

           INTEGER, DIMENSION(MAXICDOCS) :: NUMCOLS
           INTEGER, DIMENSION(MAXICDOCS) :: NUMKEYS
           INTEGER, DIMENSION(MAXICDOCS) :: NLETOLDNAM
           INTEGER, DIMENSION(MAXICDOCS) :: NEXTKEY
           CHARACTER(LEN=MAXNAM)         :: OLDNAM(MAXICDOCS)

        END MODULE DOCIC_INFO


C     ------------------------- LUNDOCWRTDAT -----------------------------
C
C    LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,NLIST,IRTFLG)
C
C    PURPOSE:    WRITE DATA FOR A SPECIFED KEY OF A DOC. FILE INTO
C                INCORE OR FILE-BASED DOC. FILE
C
C    PARAMETERS:   LUNDOC  IO UNIT                                (SENT)
C                  IKEY    NUMBER OF KEY WANTED (<0 IS COMMENTED) (RET.)
C                  DLIST   ARRAY CONTAINING NUMBERS               (RET.)
C                  NLIST   NUMBER OF ELEMENTS IN ARRAY            (SENT)
C                  IRTFLG                                         (RET.)
C
C     NOTE:    CAN NOW WRITE OUT DOC FILES WITH 100 REGISTERS
C              BUT SOME OF THE READ PGMS CAN NOT HANDLE MORE THAN 9!!
C              (FOR USE WITH 'SD C')
C
C--*********************************************************************

      SUBROUTINE LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,NLIST,IRTFLG)

      USE DOCIC_INFO
      INCLUDE 'CMBLOCK.INC' 

      REAL, DIMENSION(*)          :: DLIST
      REAL, DIMENSION(:), POINTER :: IPQ

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

      IRTFLG = 0

C     RETURN IF NLIST == 0
      IF (NLIST <= 0) RETURN

C     IKEY IS THE KEY NUMBER.
C     NLIST IS THE NUMBER OF REGISTERS (VALUES) IN THE LINE.

      IRTFLG = 1

      IF (IKEY == 0) THEN
C        WANT TO CREATE ILLEGAL KEY NUMBER
         CALL ERRT(101,'ILLEGAL DOC. FILE KEY: 0',IDUM)
         RETURN

      ELSEIF (IKEY < 0) THEN

C        WANT TO CREATE COMMENT KEY ----------------------------------
         IF (LUNDOC <= 0) THEN
C           WANT TO CREATE COMMENT KEY IN INCORE DOC. FILE 
            CALL ERRT(101,'NO COMMENT KEYS IN IN-CORE FILES',IDUM)
            RETURN
         ENDIF

         IKEYT = -IKEY
         IF (IKEYT > 9999999) THEN
           WRITE(NOUT,90) IKEY
90         FORMAT(' *** COMMENT KEY:',I8,
     &            ' MUST BE IN RANGE -9999999...-1'/)
           CALL ERRT(101,'COMMENT KEY OUT OF RANGE -9999999 .. -1',IDUM)
           RETURN
         ELSEIF (NLIST > 9) THEN
           CALL ERRT(101,'COMMENT KEY CAN ONLY HAVE 9 REGISTERS',IDUM)
           RETURN
         ENDIF

         IF (MYPID <= 0) THEN
            WRITE(LUNDOC,91) IKEYT,NLIST,(DLIST(K),K=1,NLIST)
91          FORMAT(' ;',I8,1X,I1,1X,10000(1PG13.6,1X))
            CALL FLUSHFILE(LUNDOC)
         ENDIF

      ELSEIF (LUNDOC > 0) THEN
C        WANT TO WRITE REGULAR KEY TO DISK BASED DOC. FILE -----------

         IF (NLIST > 9 .OR. IKEY > 999999) THEN
C           MUST USE NEW FORMAT
            IF (MYPID <= 0)
     &         WRITE(LUNDOC,92) IKEY,NLIST,(DLIST(K),K=1,NLIST)
92             FORMAT(I10,' ',I4,10000(' ',1PG13.6))
         
         ELSEIF (IKEY <= 99999) THEN
C           WRITE LINE OF REGISTERS WITH KEY (I5)
            IF (MYPID <= 0) THEN
               WRITE(LUNDOC,93) IKEY,NLIST,(DLIST(K),K=1,NLIST)
93             FORMAT(I5,' ',I1,10000(' ',1PG13.6))
               IF (MYPID <= 0) CALL FLUSHFILE(LUNDOC)
            ENDIF

         ELSEIF (IKEY <= 999999) THEN
C           WRITE LINE OF REGISTERS WITH KEY (I6)
            IF (MYPID <= 0) THEN
               WRITE(LUNDOC,94) IKEY,NLIST,(DLIST(K),K=1,NLIST)
94             FORMAT(I6,' ',I1,10000(' ',1PG13.6))
               IF (MYPID <= 0) CALL FLUSHFILE(LUNDOC)
            ENDIF
         ENDIF

      ELSEIF (LUNDOC < 0) THEN
C        WANT TO WRITE REGULAR KEY TO INCORE DOC. FILE --------------

C        GET ARRAY SIZE FOR INCORE FILE (FIXED WHEN IT WAS CREATED)
         IC   = - LUNDOC
         IF (IC > MAXICDOCS) THEN
C           FILE LIST INDEX OUT OF RANGE
            CALL ERRT(102,'MAX. INCORE DOC. FILE. NUMBER',MAXICDOCS)
            RETURN
         ENDIF

         MAXX = NUMCOLS(IC)
         MAXY = NUMKEYS(IC)
         IF (NLIST > (MAXX - 1)) THEN
            WRITE(NOUT,95) NLIST
95          FORMAT('  *** NUMBER OF REGISTERS: ',I3)
            CALL ERRT(102,'REGISTER LIMIT FOR THIS INCORE FILE',MAXX-1)
            RETURN

         ELSEIF (IKEY > MAXY) THEN
            WRITE(NOUT,96) IKEY
96          FORMAT('  *** KEY:',I10)
            CALL ERRT(102,'KEY LIMIT FOR THIS INCORE FILE',MAXY)
            RETURN
         ENDIF

         IPQ       => LOCDOC(IC)%IPT
         ILOC      = (IKEY - 1) * MAXX + 1
         IPQ(ILOC) = NLIST

         DO IREG=1,NLIST
            IPQ(ILOC+IREG) = DLIST(IREG)
         ENDDO
      ENDIF

      IRTFLG = 0

      RETURN
      END




C     ------------------------- LUNDOCWRTDATF ----------------------------
C
C    LUNDOCWRTDATF(LUNDOC,IKEY,DLIST,NLIST,IRTFLG)
C
C    PURPOSE:    WRITE DATA FOR A SPECIFED KEY OF A DOC. FILE INTO
C                INCORE OR FILE-BASED DOC. FILE. USER SUPPLIES
C                LINE FORMAT.
C
C    PARAMETERS:   LUNDOC  IO UNIT                                (SENT)
C                  IKEY    NUMBER OF KEY WANTED (<0 IS COMMENTED) (RET.)
C                  DLIST   ARRAY CONTAINING NUMBERS               (RET.)
C                  NLIST   NUMBER OF ELEMENTS IN ARRAY            (SENT)
C                  FORMOUT LINE OUTPUT FORMAT                     (SENT)
C                  IRTFLG                                         (RET.)
C
C     NOTE:    CAN WRITE OUT DOC FILES WITH 100 REGISTERS
C              BUT SOME OF THE READ PGMS CAN NOT HANDLE MORE THAN 9!!
C
C--*********************************************************************

      SUBROUTINE LUNDOCWRTDATF(LUNDOC,IKEY,DLIST,NLIST,
     &                         FORMOUT,IRTFLG)

      USE DOCIC_INFO
      INCLUDE 'CMBLOCK.INC' 

      INTEGER                     :: LUNDOC,IKEY 
      REAL, DIMENSION(*)          :: DLIST
      INTEGER                     :: NLIST
      CHARACTER(LEN=*)            :: FORMOUT
      INTEGER                     :: IRTFLG
 
      REAL, DIMENSION(:), POINTER :: IPQ

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

      IRTFLG = 0

C     IKEY IS THE KEY NUMBER.
C     NLIST IS THE NUMBER OF REGISTERS (VALUES) IN THE LINE.

C     RETURN IF NLIST == 0
      IF (NLIST <= 0) RETURN

      IRTFLG = 1

      IF (IKEY == 0) THEN
C        WANT TO CREATE ILLEGAL KEY NUMBER
         CALL ERRT(101,'ILLEGAL DOC. FILE KEY: 0',IDUM)
         RETURN

      ELSEIF (IKEY < 0) THEN

C        WANT TO CREATE COMMENT KEY ----------------------------------
         IF (LUNDOC <= 0) THEN
C           WANT TO CREATE COMMENT KEY IN INCORE DOC. FILE 
            CALL ERRT(101,'NO COMMENT KEYS IN IN-CORE FILES',IDUM)
            RETURN
         ENDIF

         IKEYT = -IKEY
         IF (IKEYT > 9999999) THEN
           CALL ERRT(102,'COMMENT KEY OUT OF RANGE -9999999 .. -1',IKEY)
           RETURN
         ELSEIF (NLIST > 9) THEN
           CALL ERRT(101,'COMMENT KEY CAN ONLY HAVE 9 REGISTERS',IDUM)
           RETURN
         ENDIF

         IF (MYPID <= 0) THEN
            WRITE(LUNDOC,91) IKEYT,NLIST,(DLIST(K),K=1,NLIST)
91          FORMAT(' ;',I8,1X,I1,1X,10000(1PG13.6,1X))
            CALL FLUSHFILE(LUNDOC)
         ENDIF

      ELSEIF (LUNDOC > 0) THEN
C        WANT TO WRITE REGULAR KEY TO DISK BASED DOC. FILE -----------

         IF (MYPID <= 0) THEN
            WRITE(LUNDOC,FORMOUT) IKEY,NLIST,(DLIST(K),K=1,NLIST)
         
            IF (MYPID <= 0) CALL FLUSHFILE(LUNDOC)
         ENDIF

      ELSEIF (LUNDOC < 0) THEN
C        WANT TO WRITE REGULAR KEY TO INCORE DOC. FILE --------------

C        GET ARRAY SIZE FOR INCORE FILE (FIXED WHEN IT WAS CREATED)
         IC   = - LUNDOC
         IF (IC > MAXICDOCS) THEN
C           FILE LIST INDEX OUT OF RANGE
            CALL ERRT(102,'MAX. INCORE DOC. FILE. NUMBER',MAXICDOCS)
            RETURN
         ENDIF

         MAXX = NUMCOLS(IC)
         MAXY = NUMKEYS(IC)
         IF (NLIST > (MAXX - 1)) THEN
            WRITE(NOUT,95) NLIST
95          FORMAT('  *** NUMBER OF REGISTERS: ',I3)
            CALL ERRT(102,'REGISTER LIMIT FOR THIS INCORE FILE',MAXX-1)
            RETURN

         ELSEIF (IKEY > MAXY) THEN
            WRITE(NOUT,96) IKEY
96          FORMAT('  *** KEY:',I10)
            CALL ERRT(102,'KEY LIMIT FOR THIS INCORE FILE',MAXY)
            RETURN
         ENDIF

         IPQ       => LOCDOC(IC)%IPT
         ILOC      = (IKEY - 1) * MAXX + 1
         IPQ(ILOC) = NLIST

         DO IREG=1,NLIST
            IPQ(ILOC+IREG) = DLIST(IREG)
         ENDDO
      ENDIF

      IRTFLG = 0

      RETURN
      END




C     ------------------------- LUNDOCREDLIN -----------------------------
C
C    LUNDOCREDLIN(LUNDOC,WANTERRT,WARNIT,DBUF,MAXX,MAXY,WANTICCOL,
C                 IKEY,ICOUNT,IRTFLG
C
C    PURPOSE:    RECOVER LINE OF  DATA FROM A DOC FILE AND RETURN IT 
C                IN DBUF, DBUF KEEPS ICOUNT IN FIRST COLUMN.
C 
C    PARAMETERS:  LUNDOC     IO UNIT                             (SENT)
C                 WANTERRT   CALLERRT FLAG                       (SENT)
C                 WARNIT     WARNIT FLAG                    (SENT/RET.)
C                 DBUF       DATA ARRAY                     (SENT/RET.)
C                 MAXX       MAX X ARRAY DIMENSION               (SENT)
C                            (NUMBER OF REGISTERS + 1)!!
C                 MAXY       MAX Y ARRAY DIMENSION               (SENT)
C                            IF ZERO DOES NOT LOAD ARRAY JUST LINE
C                 WANTICCOLC WANT ICOUNT IN DBUF                 (SENT)
C                 IKEY       KEY  RETRIEVED                      (RET.)
C                 ICOUNT     COLS. RETRIEVED                     (RET.)
C                 IRTFLG     1=ERROR, 0 = NORMAL                 (RET.)
C
C    NOTE: KEY & ICOUNT MUST BE IN FIRST 80 COL.  SHOULD HANDLE 
C          OLD FORMAT DATA OK. NO LIMIT ON ICOUNT 
C 
C--*********************************************************************

      SUBROUTINE LUNDOCREDLIN(LUNDOC,WANTERRT,WARNIT,WANTICCOL,
     &                        DBUF,MAXX,MAXY,IKEY,ICOUNT,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 

      INTEGER, PARAMETER       :: LENRECLIN = 80  ! ONLY NEEDS START
      REAL,DIMENSION(*)        :: DBUF
      CHARACTER(LEN=LENRECLIN) :: RECLIN
      LOGICAL                  :: WARNIT,WANTERRT,NEWFORM,WANTICCOL

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

      IF (MYPID <= 0) THEN
         READ(LUNDOC,84,IOSTAT=IERR) RECLIN
84       FORMAT(A)
         IF (IERR < 0) RECLIN(1:2) = '!!'
      ENDIF
#ifdef USE_MPI
      CALL BCAST_MPI('LUNDOCREDLIN','RECLIN',RECLIN,LENRECLIN,'C',
     &               ICOMM,MPIERR)
#endif 

      IKEY   = 0
      ICOUNT = 0
      IRTFLG = 0

      IF (RECLIN(1:2) == '!!') THEN
C         END OF FILE
          IRTFLG = -1
          RECLIN = ''
          RETURN
       ENDIF

C      IGNORE COMMENT LINES & COMMENT KEY LINES. 
       IF (RECLIN(2:2) == ';') RETURN
       IF (RECLIN(1:2) == '; ') RETURN  !BAD JWEB COMMENT LINE
          
       NEWFORM = .TRUE.
       READ(RECLIN,*,IOSTAT=IERR) IKEY,ICOUNT

       IF (IERR > 0) THEN
C         ERROR ON READ, TRY OLD DOC. FILE FORMAT

          READ(RECLIN,83,IOSTAT=IERR) IKEY,ICOUNT
83        FORMAT(I6,I1,10000F12.6)
          NEWFORM = .FALSE.

          IF (IERR > 0) THEN
C            ERROR ON READ USING OLD FORMAT ALSO, RETURN
             IF (WANTERRT) THEN
                NLET = lnblnkn(RECLIN)
                IF (NLET <= 0) THEN
                   CALL ERRT(101,'EMPTY DOC FILE LINE',NE)
                ELSE
                   WRITE(NOUT,94) RECLIN(1:NLET)
94                 FORMAT('  *** UNABLE TO INTERPRET DOC FILE LINE ',
     &                    'STARTING WITH: ',A)
                   CALL ERRT(100,'LUNDOCREDLIN',NE)
                ENDIF
             ENDIF
             ICOUNT = 0
             IKEY   = 0
             IRTFLG = 1
             RETURN
          ENDIF
       ENDIF

       IF (ICOUNT <= 0) THEN
           IF (MYPID <= 0)
     &         WRITE(NOUT,*)' EMPTY DOCUMENT FILE LINE SKIPPED'
           RETURN

        ELSEIF (IKEY < 0) THEN
           IF (MYPID <= 0)
     &         WRITE(NOUT,*)' CONTINUATION LINE SKIPPED IN DOC FILE'
           RETURN

        ELSEIF (IKEY == 0) THEN
C          ZERO KEY ILLEGAL, PRINT WARNING MSG. BUT CONTINUE
           IF (MYPID <= 0)
     &        WRITE(NOUT,*)' ILLEGAL KEY NUMBER: 0  SKIPPED IN DOC FILE'
           RETURN

        ELSEIF (MAXY > 0 .AND. IKEY > MAXY) THEN
C          KEY THAT WILL NOT FIT IN DBUF SENDS ERROR MSG.
           IF (WARNIT) THEN 
               IF (MYPID <= 0)  WRITE(NOUT,93) MAXY
93             FORMAT('  KEYS GREATER THAN: ',I7,' NOT RETRIEVED')
               WARNIT = .FALSE.
           ENDIF
           IKEY   = 0
           ICOUNT = 0
           RETURN
        ENDIF

C        AVOID ARRAY OVERFLOW
         ICOUNT = MIN(ICOUNT,MAXX-1)

         IF (WANTICCOL) THEN
C           FIND DBUF LOCATION POINTER
            ILOC       = (IKEY -1) * MAXX + 1
            DBUF(ILOC) = ICOUNT
         ELSE
C           FIND DLIST LOCATION POINTER
            ILOC = 0
         ENDIF

         IF (NEWFORM .AND. MYPID <= 0) THEN
C            READ REGISTERS USING NEW DOC. FILE FORMAT
             BACKSPACE(LUNDOC)
             READ(LUNDOC,*,IOSTAT=IERR)IKEYT,ICOUNTT,
     &                                 (DBUF(ILOC+I),I=1,ICOUNT)
             IF (IERR .NE. 0) THEN
C               TRY READING AGAIN USING OLD FORMAT
                BACKSPACE(LUNDOC)
                READ(LUNDOC,83,IOSTAT=IERR) IKEYT,ICOUNTT,
     &                                      (DBUF(ILOC+I),I=1,ICOUNT)
             ENDIF
         ELSEIF (.NOT. NEWFORM  .AND. MYPID <= 0) THEN
C            READ REGISTERS USING OLD DOC. FILE FORMAT
             BACKSPACE(LUNDOC)
             READ(LUNDOC,83,IOSTAT=IERR) IKEYT,ICOUNTT,
     &                                   (DBUF(ILOC+I),I=1,ICOUNT)
         ENDIF
#ifdef USE_MPI
         CALL BCAST_MPI('LUNDOCREDLIN','DBUF',DBUF(ILOC+1),
     &                  ICOUNT,'R',ICOMM)
         CALL BCAST_MPI('LUNDOCREDLIN','IERR',IERR,1,'I',ICOMM)
#endif

         IF (IERR .NE. 0) THEN
C            ERROR ON REGISTER COLS READ
             IF (WANTERRT) THEN
                IF (MYPID <= 0) WRITE(NOUT,94) RECLIN
                CALL ERRT(100,'LUNDOCREDLIN',NE)
                RETURN
             ENDIF
             IKEY   = 0
             ICOUNT = 0
             IRTFLG = 1
             RETURN
         ENDIF

         END

C     ------------------------- LUNDOCREDDAT -----------------------------
C
C    LUNDOCREDDAT(LUNDOC,IKEY,DLIST,NMAX,ICOUNT,TILLEND,GOBACK,IRTFLG)
C
C    PURPOSE:    RECOVER A SPECIFED KEY FROM A DOC FILE.
C                WORKS ON INCORE AND FILE BASED DOC. FILES.
C
C    PARAMETERS:   LUNDOC  IO UNIT                               (SENT)
C                  IKEY    NUMBER OF KEY WANTED (<0 IS COMMENTED)(RET.)
C                  DLIST   ARRAY CONTAINING REGISTER NUMBERS     (RET.)
C                  NMAX    MAX DLIST ARRAY DIMENSION             (SENT)
C                  ICOUNT  NUMBER OF ELEMENTS IN ARRAY           (RET.)
C                  TILLEND KEEP READING TILL EOF EVEN IF FOUND   (SENT.)
C                  GOBACK  REWIND AND READ AGAIN IF NOT FOUND    (SENT)
C                  IRTFLG  1=ERROR                               (RET.)
C                          2=NOT FOUND                           (RET.)
C
C--*********************************************************************

      SUBROUTINE LUNDOCREDDAT(LUNDOC,IKEY,DLIST,NMAX,ICOUNT,
     &                        TILLEND,GOBACK,IRTFLG)

      USE DOCIC_INFO
      INCLUDE 'CMBLOCK.INC'

      REAL                        :: DLIST(NMAX)
      LOGICAL                     :: TILLEND,GOBACK,FIRST
      LOGICAL                     :: WANTERRT,WARNIT
      REAL, POINTER               :: IPQ(:)

      REAL                        :: DLISTT(NMAX)

      IRTFLG   = 1

      IF (IKEY < 0) THEN
C        WANT A COMMENTED KEY
         ICOUNT = NMAX
         CALL LUNDOCGETCOM(LUNDOC,IKEY,DLIST,ICOUNT,TILLEND,IRTFLG)
         RETURN

      ELSEIF (IKEY == 0) THEN
C        THIS IS ILLEGAL KEY, SOUND WARNING
         WRITE(NOUT,90) 
90       FORMAT('  WARNING; RETRIEVING ILLEGAL KEY NUMBER: 0')

      ELSEIF (LUNDOC < 0) THEN
C        WANT TO READ REGULAR KEY FROM INCORE DOC. FILE

C        GET ARRAY SIZE FOR INCORE FILE (FIXED WHEN IT WAS CREATED)
         IC   = - LUNDOC
         IF (IC > MAXICDOCS) THEN
C           FILE LIST INDEX OUT OF RANGE
            CALL ERRT(102,'MAX. INCORE DOC. FILE. NUMBER',MAXICDOCS)
            RETURN
         ENDIF

C        FIND HIGHEST KEY NUMBER IN THE INCORE DOC. FILE
         MAXY   = NUMKEYS(IC)
         MAXXIC = NUMCOLS(IC) 

         IF (IKEY > MAXY) THEN
            WRITE(NOUT,96) IKEY
96          FORMAT('  *** KEY:',I10)
            CALL ERRT(102,'KEY LIMIT FOR THIS INCORE FILE',MAXY)
            RETURN
         ENDIF

         IPQ     => LOCDOC(IC)%IPT
         ILOC    = (IKEY - 1) * MAXXIC + 1

C        MAKE SURE DLIST DOES NOT OVERFLOW
         ICOUNT  = IPQ(ILOC)
         ICOUNT  = MIN(NUMCOLS(IC)-1,ICOUNT,NMAX)

         IF (ICOUNT > 0) THEN
c           THIS KEY USED IN INCORE DOC. FILE

            DO IREG=1,ICOUNT
               DLIST(IREG) = IPQ(ILOC+IREG)
            ENDDO
       
            IRTFLG = 0
         ELSE
            IRTFLG = 2
         ENDIF
         RETURN
      ENDIF

      ICOUNT   = 0
      FIRST    = .TRUE.
      WARNIT   = .TRUE.
      WANTERRT = .FALSE.

      DO
C        READ NEXT LINE FROM DOC FILE
         CALL LUNDOCREDLIN(LUNDOC,WANTERRT,WARNIT,.FALSE.,
     &                     DLISTT,NMAX+1,0,IKEYT,ICOUNTT,IRTFLG)

         IF (IRTFLG < 0) THEN
C           END OF FILE
            IF (ICOUNT == 0 .AND. GOBACK .AND. FIRST) THEN
C              DID NOT FIND KEY, REWIND AND TRY AGAIN
               REWIND(LUNDOC)
               FIRST = .FALSE.
               CYCLE
            ENDIF

C           REPOSISTION LOCATION
            IF (TILLEND) REWIND(LUNDOC)

C           RETURN 2 ON EOF IF DID NOT FIND KEY
            IRTFLG = 2
            IF (ICOUNT > 0) IRTFLG = 0
            RETURN

         ELSEIF (IRTFLG > 0) THEN
C           ERROR RETRIEVING ARRAY
            RETURN
         ENDIF

         IF (IKEY == IKEYT) THEN
C           FOUND DESIRED KEY, MAY HAVE DUPLICATES OF KEY
            ICOUNT = ICOUNTT
            DLIST  = DLISTT

C           IF "TILLEND", KEEP READING TILL EOF
            IF (TILLEND) CYCLE
            RETURN
         ENDIF
      ENDDO 

      END

C     ------------------------- LUNDOCREDALL -----------------------------
C
C    LUNDOCREDALL(LUNDOC,DBUF,MAXX,MAXY,WANTERRT,NGOT,IRTFLG)
C
C    PURPOSE:    RECOVER ALL DATA FROM A DOC FILE AND RETURN IT IN DBUF
C                DBUF KEEPS KEY IN FIRST COLUMN.
C 
C    PARAMETERS:  LUNDOC     IO UNIT                             (SENT)
C                 DBUF       DATA ARRAY (MUST ALREADY EXIST)     (RET.)
C                 MAXX       MAX X ARRAY DIMENSION               (SENT)
C                            (NUMBER OF REGISTERS + 1)!!
C                 MAXY       MAX Y ARRAY DIMENSION               (SENT)
C                 WANTERRT   CALLERRT FLAG                       (SENT)
C                 NGOT       MAX. KEY # RETRIEVED                (RET.)
C                 IRTFLG     1=ERROR, 0 = NORMAL                 (RET.)
C
C--*********************************************************************

      SUBROUTINE LUNDOCREDALL(LUNDOC,DBUF,MAXX,MAXY,WANTERRT,
     &                        NGOT,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 

      REAL,DIMENSION(MAXX,MAXY) :: DBUF
      LOGICAL                   :: WARNIT,WANTERRT

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

      IF (LUNDOC < 0) THEN
C        RECOVER FROM INCORE DOC FILE
         IC = -LUNDOC
         CALL LUNDOCREDALLI(IC,DBUF,MAXX,MAXY,NGOTX,NGOT,IRTFLG)
         RETURN
      ENDIF

      IRTFLG = 1
      NGOT   = 0
      WARNIT = .TRUE.

C     CLEAR DBUF RETURNED ANSWER BUFFER KEY NUMBERS
      DBUF(1,:) = 0.0

C     WILL LOOP UNTIL IT RETURNS EOF OR UNTIL ERROR
      DO
         CALL LUNDOCREDLIN(LUNDOC,WANTERRT,WARNIT,.TRUE.,
     &                     DBUF,MAXX,MAXY,IKEY,ICOUNT,IRTFLG)

         IF (IRTFLG < 0) THEN
C           FINISHED READING FILE
            IRTFLG = 0
#ifdef USE_MPI
            !if (mypid <=0) write(6,*) 'End of doc file, mpi barrier'
            CALL MPI_BARRIER(ICOMM,MPIERR)
#endif
            RETURN

         ELSEIF (IRTFLG > 0) THEN
C           ERROR RETRIEVING ARRAY
            RETURN
         ENDIF

         IF (ICOUNT > 0 .AND. 
     &       (ICOUNT +2) <= MAXX .AND.
     &        IKEY <= MAXY) THEN
C            ZERO REMAINING VALUES IN REGISTER LIST FOR THIS LINE
             DBUF(ICOUNT+2:MAXX,IKEY) = 0.0
         ENDIF
         NGOT = MAX(NGOT,IKEY)
      ENDDO 

      END


C     ------------------------- LUNDOCREDALLI -----------------------------
C
C    LUNDOCREDALLI(IC,DBUF,MAXX,MAXY,IGOX,IENDX,NGOTX,NGOTY,IRTFLG)
C
C    PURPOSE:     RECOVER DATA FROM INCORE DOC FILE AND RETURN 
C                 IT IN DBUF. FIRST COLUMN IN DBUF IS KEY INUSE FLAG.
C                
C    PARAMETERS:  IC         IC UNIT                             (SENT)
C                 DBUF       DATA ARRAY (MUST ALREADY EXIST)     (RET.)
C                 MAXX       MAX X ARRAY DIMENSION               (SENT)
C                            (NUMBER OF REGISTERS + 1)!!
C                 MAXY       MAX Y ARRAY DIMENSION               (SENT)
C                 NGOTX      MAX. REG # RETRIEVED                (RET.)
C                 NGOTY      MAX. KEY # RETRIEVED                (RET.)
C                 IRTFLG     1=ERROR, 0 = NORMAL                 (RET.)
C
C--*********************************************************************

       SUBROUTINE LUNDOCREDALLI(IC,DBUF,MAXX,MAXY,NGOTX,NGOTY,IRTFLG)

       USE DOCIC_INFO

       INCLUDE 'CMBLOCK.INC' 

       REAL,DIMENSION(MAXX,MAXY)   :: DBUF
       REAL, DIMENSION(:), POINTER :: IPQ

       IRTFLG = 1

       IF (IC > MAXICDOCS) THEN
C         FILE LIST INDEX OUT OF RANGE
          CALL ERRT(102,'MAX. INCORE DOC. FILE. NUMBER',MAXICDOCS)
          RETURN
       ENDIF

C      FIND HIGHEST KEY NUMBER IN THE INCORE DOC. FILE
       MAXYIC = NUMKEYS(IC)
       MAXXIC = NUMCOLS(IC) 

       NGOTX = MIN(MAXX-1,MAXXIC-1)
       WRITE(NOUT,91) NGOTX 
91     FORMAT('  Number of incore registers recovered: ',I10)

       NGOTY  = MIN(MAXY,MAXYIC)
       WRITE(NOUT,90) NGOTY 
90     FORMAT('  Number of incore keys recovered: ',I10)

       IPQ       => LOCDOC(IC)%IPT
       ILOCD     = 0

       DO IKEY = 1,NGOTY
          DO ICOL = 1,NGOTX+1
             ILOCI           = (IKEY - 1) * MAXXIC + ICOL
             DBUF(ICOL,IKEY) = IPQ(ILOCI)
          ENDDO 
  
          IF ((NGOTX+1) < MAXX) THEN
C            FILL REST OF DBUF LINE WITH ZEROS
             DBUF(NGOTX+1:MAXX, IKEY) = 0.0
          ENDIF
       ENDDO

       IF (NGOTY < MAXY) THEN
C         FILL REST OF DBUF WITH ZEROS
          DBUF(1:MAXX, NGOTY+1:MAXY) = 0.0
       ENDIF

       IRTFLG = 0
       RETURN

       END


C     ------------------------- LUNDOCREDANG -----------------------------
C
C    LUNDOCREDANG(LUNDOC,DLIST,MAXY,NGOTY,MAXGOTY,IRTFLG)
C
C    PURPOSE:    RECOVER COLUMNS 1-3 OF DOC FILE AS REAL NUMBERS. IN
C                AN ARRAY. RETURNS NUMBER OF FILLED KEYS AND HIGHEST 
C                KEY THAT IS FILLED.
C                
C    PARAMETERS: LUNDOC     I/O UNIT                            (SENT)
C                DLIST      DATA ARRAY (MUST ALREADY EXIST)     (RET.)
C                MAXY       ARRAY DIMENSION                     (SENT)
C                NGOTY      NUMBER OF FILLED KEYS IN LIST       (RET.)
C                MAXGOTY    HIGHEST FILLED KEY IN LIST          (RET.)
C                IRTFLG     1=ERROR, 0 = NORMAL                 (RET.)
C
C--*********************************************************************

      SUBROUTINE LUNDOCREDANG(LUNDOC,DLIST,MAXY,NGOTY,MAXGOTY,IRTFLG)

      CALL LUNDOCREDSLC(LUNDOC,.FALSE.,IDUM,DLIST, 3,MAXY,
     &    .TRUE.,.FALSE.,1,3, 1,MAXY, NGOTY,MAXGOTY,IRTFLG)

      RETURN
      END



C     ------------------------- LUNDOCREDSEL -----------------------------
C
C    LUNDOCREDSEL(LUNDOC,ILIST,MAXY,NGOTY,MAXGOTY,IRTFLG)
C
C    PURPOSE:    RECOVER COLUMN 1 OF DOC FILE AS INTEGER NUMBERS IN
C                AN ARRAY. RETURNS NUMBER OF FILLED KEYS AND HIGHEST 
C                KEY THAT WAS FILLED.
C                
C    PARAMETERS: LUNDOC     I/O UNIT                            (SENT)
C                ILIST      DATA ARRAY (MUST ALREADY EXIST)     (RET.)
C                MAXY       ARRAY DIMENSION                     (SENT)
C                NGOTY      NUMBER OF FILLED KEYS IN LIST       (RET.)
C                MAXGOTY    HIGHEST FILLED KEY IN LIST          (RET.)
C                IRTFLG     1=ERROR, 0 = NORMAL                 (RET.)
C
C--*********************************************************************

      SUBROUTINE LUNDOCREDSEL(LUNDOC,ILIST,MAXY,NGOTY,MAXGOTY,IRTFLG)

      CALL LUNDOCREDSLC(LUNDOC,.TRUE.,ILIST,DUM,1,MAXY,
     &    .FALSE.,.FALSE. ,1,1, 1,MAXY,NGOTY,MAXGOTY,IRTFLG)

      RETURN
      END


C     ------------------------- LUNDOCREDSEQ -----------------------------
C
C    LUNDOCREDSEQ(LUNDOC,DLIST,MAXX,MAXY,IWANTY,NGOTY,IRTFLG)
C
C    PURPOSE:    RECOVER COLUMN 1..MAXX  DOC FILE AS REAL NUMBERS IN
C                AN ARRAY. RETURNS NUMBER OF FILLED KEYS.  ARRAY 
C                RECOVERS ALL KEYS EVEN IF DUPLICATED.
C                
C    PARAMETERS: LUNDOC     I/O UNIT                            (SENT)
C                DLIST      DATA ARRAY (MUST ALREADY EXIST)     (RET.)
C                MAXX,MAXY  ARRAY DIMENSION                     (SENT)
C                NGOTY      NUMBER OF FILLED KEYS IN DLIST      (RET.)
C                IRTFLG     1=ERROR, 0 = NORMAL                 (RET.)
C
C--*********************************************************************

      SUBROUTINE LUNDOCREDSEQ(LUNDOC,DLIST,MAXX,MAXY,IWANTY,
     &                        NGOTY,IRTFLG)

      CALL LUNDOCREDSLC(LUNDOC,.FALSE.,IDUM,DLIST,MAXX,MAXY,
     &    .FALSE.,.FALSE. ,1,MAXX, 1,IWANTY, NGOTY,MAXGOTY,IRTFLG)

      RETURN
      END




C     ------------------------- LUNDOCREDSLC -----------------------------
C
C    LUNDOCREDSLC(LUNDOC,USEINT,ILIST,DLIST,MAXX,MAXY,
C             KEEPKEYS, ERRSKIP,IGOX,IENDX, IGOY,IENDY, 
C             NGOTY,MAXGOTY,IRTFLG)
C
C    PURPOSE:    RECOVER SLICE OF SPECIFED KEYS AND REGISTERS FROM AN
C                INCORE DOC FILE AND RETURNS SLICE IN LIST OR DLIST. 
C                SELECTS WHETHER REGISTER CONTENTS ARE RETURNED AS 
C                INTEGERS OR REAL, AND WHETHER TO KEEP KEYS AS
C                INDICES OR JUST RECOVER SEQUENTIALLY IGNORING KEYS.
C                
C    PARAMETERS: LUNDOC     FILE I/O UNIT                       (SENT)
C                USEINT     SELECTS ILIST OR DLIST              (SENT)
C                ILIST      DATA ARRAY (MUST ALREADY EXIST)     (RET.)
C                DLIST      DATA ARRAY (MUST ALREADY EXIST)     (RET.)
C                MAXX       ARRAY DIMENSION                     (SENT)
C                MAXY       ARRAY DIMENSION                     (SENT)
C                KEEPKEYS   KEEPS KEYS AS INDICES               (SENT)
C                ERRSKIP    WARN IF KEY MISSING                 (SENT)
C                IGOX       STARTING REG.                       (SENT)
C                IENDX      ENDING REG.                         (SENT)
C                IGOY       STARTING KEY                        (SENT)
C                IENDY      ENDING KEY                          (SENT)
C                NGOTY      NUMBER OF FILLED KEYS IN LIST       (RET.)
C                MAXGOTY    HIGHEST FILLED KEY IN LIST          (RET.)
C                IRTFLG     1=ERROR, 0 = NORMAL                 (RET.)
C
C--*********************************************************************

      SUBROUTINE LUNDOCREDSLC(LUNDOC,USEINT,ILIST,DLIST,MAXX,MAXY,
     &    KEEPKEYS,ERRSKIP,IGOX,IENDX, IGOY,IENDY, NGOTY,MAXGOTY,IRTFLG)

      USE DOCIC_INFO

      INCLUDE 'CMBLOCK.INC' 

      INTEGER,DIMENSION(MAXX,MAXY) :: ILIST
      REAL,DIMENSION(MAXX,MAXY)    :: DLIST
      REAL, DIMENSION(MAXX)        :: PLIST
      LOGICAL                      :: USEINT,KEEPKEYS,ERRSKIP
      LOGICAL                      :: WANTERRT,WARNIT
      REAL, DIMENSION(:), POINTER  :: IPQ
 
      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

      IRTFLG    = 1
      NGOTX     = IENDX - IGOX + 1
      MAXGOTY   = 0
      NWANTY    = IENDY - IGOY + 1
      NGOTY     = 0

      IF (LUNDOC < 0) THEN
C        RECOVER SLICE FROM INCORE DOC FILE --------------------------
         IC = -LUNDOC
      
         IF (IC > MAXICDOCS) THEN
C           FILE LIST INDEX OUT OF RANGE
            CALL ERRT(102,'MAX. INCORE DOC. FILE. NUMBER',MAXICDOCS)
            RETURN
         ENDIF

C        FIND HIGHEST KEY NUMBER IN THE INCORE DOC. FILE
         MAXXIC = NUMCOLS(IC) 
         MAXYIC = NUMKEYS(IC)
         IPQ    => LOCDOC(IC)%IPT

         IENDYNOW = MIN(IENDY,MAXYIC)
         IF (IGOX > IENDX .OR. IENDX > (MAXXIC-1)) THEN
            CALL ERRT(102,'HIGHEST AVAILABLE INCORE REGISTER',MAXXIC-1)
            RETURN

         ELSEIF (IGOY > IENDYNOW) THEN
            CALL ERRT(102,'HIGHEST AVAILABLE INCORE KEY',MAXYIC)
            RETURN
         ENDIF

         DO IKEY = IGOY,IENDYNOW
            ILOCIC = (IKEY - 1) * MAXXIC + 1
            ICOUNT = IPQ(ILOCIC)

            IF (.NOT. KEEPKEYS) THEN
C              DO NOT SAVE MT LINE
               IF (ICOUNT <= 0) CYCLE

C              STORE CONSECUTIVELY
               NGOTY = NGOTY + 1
               IKEYT = NGOTY

               MAXGOTY = MAX(MAXGOTY,IKEY)
            ELSE
               NGOTY = NGOTY + 1
               IF (ERRSKIP .AND. ICOUNT <= 0) THEN
                  CALL ERRT(102,'MISSING KEY',IKEYT)
                  RETURN
               ENDIF
               IKEYT   = IKEY
               MAXGOTY = IKEY
            ENDIF

            IF (USEINT) THEN
               ILIST(IGOX:IENDX,IKEYT) = IPQ(ILOCIC+1:ILOCIC+NGOTX)
            ELSE
               DLIST(IGOX:IENDX,IKEYT) = IPQ(ILOCIC+1:ILOCIC+NGOTX)
            ENDIF
         ENDDO

      ELSE
C        RECOVER SLICE FROM FILE BASED DOC FILE -----------------------

         WARNIT   = .FALSE.
         WANTERRT = .TRUE.
         DO
C           READ NEXT LINE FROM DOC FILE
            CALL LUNDOCREDLIN(LUNDOC,WANTERRT,WARNIT,.FALSE.,
     &                     PLIST,MAXX+1,0,IKEY,ICOUNT,IRTFLG)

            IF (IRTFLG < 0) THEN
C              END OF FILE
               GOTO 999

            ELSEIF (IRTFLG > 0) THEN
C              ERROR RETRIEVING ARRAY
               RETURN
            ENDIF

            IF (ICOUNT > 0 .AND.
     &         (IKEY >= IGOY .AND. IKEY <= IENDY)) THEN
C              ONLY RETRIEVE REGISTERS FROM KEYS: IGOY...IENDY
               MAXGOTY = MAX(MAXGOTY,IKEY)
               NGOTY   = NGOTY + 1
               IF (.NOT. KEEPKEYS) IKEY = NGOTY

               IF (USEINT) THEN
                  ILIST(1:NGOTX,IKEY) = PLIST(IGOX:IENDX)
               ELSE
                  DLIST(1:NGOTX,IKEY) = PLIST(IGOX:IENDX)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

999   IF (MYPID <= 0) THEN
         WRITE(NOUT,91) NGOTY 
91       FORMAT('  NUMBER OF KEYS RECOVERED: ',I7)
      ENDIF

      IF (KEEPKEYS .AND. ERRSKIP .AND. NGOTY < NWANTY) THEN
         CALL ERRT(102,'MISSING KEYS',NWANTY - NGOTY)
      ENDIF

      IRTFLG = 0
      RETURN

      END


C     ------------------------- LUNDOCGETKEY -----------------------------
C
C    LUNDOCGETKEY(LUNDOC,DBUF,MAXX,MAXY,IKEY,PLIST,NLIST,IRTFLG)
C
C    PURPOSE:    RECOVERS SPECIFIC LINE (KEY) FROM DBUF. DBUF HAS 
C                BEEN RECOVERED FROM DOC FILE PRIOR TO CALLING THIS
C                SUBROUTINE.
C
C    PARAMETERS:  LUNDOC  IO UNIT                                (SENT)
C                 DBUF    DATA ARRAY                             (RET.)
C                 MAXX    MAX X ARRAY DIMENSION                  (SENT)
C                 MAXY    MAX Y ARRAY DIMENSION                  (SENT)
C                 PLIST   LIST OF VALUES                         (RET.)
C                 NLIST   NUMBER OF VALUES IN PLIST          (SENT/RET.)
C                 WANTERRT  CALLERRT FLAG                        (SENT)
C                 IRTFLG  1=ERROR, 0 = NORMAL                    (RET.)
C
C--*********************************************************************


      SUBROUTINE LUNDOCGETKEY(LUNDOC,DBUF,MAXX,MAXY,IKEY,PLIST,NLIST,
     &                        WANTERRT,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 

      DIMENSION         DBUF(*),PLIST(*)
      LOGICAL        :: WANTERRT

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

      IRTFLG = 1
      IF (IKEY > MAXY) THEN
C        KEY NOT FOUND IN DBUF
         IF (MYPID <= 0) WRITE(NOUT,90) IKEY
90       FORMAT('  *** KEY:',I7,'  NOT FOUND')
         IF (WANTERRT) CALL ERRT(100,'LUNDOCGETKEY',IDUM)
         RETURN

      ELSEIF (IKEY <= 0) THEN
         IF (MYPID <= 0) WRITE(NOUT,91) IKEY
91       FORMAT('  *** INVALID DOC FILE KEY REQUESTED: ',I7,/)
         IF (WANTERRT) CALL ERRT(100,'LUNDOCGETKEY',IDUM)
         RETURN
      ENDIF

C     FIND DBUF LOCATION POINTER
      ILOC = (IKEY - 1) * MAXX + 1

      NLISTGOT = DBUF(ILOC)
      IF (NLISTGOT <= 0) THEN
C        KEY NOT FOUND IN DBUF
         IF (WANTERRT) CALL ERRT(102,'EMPTY LINE FOR DOC FILE KEY',IKEY)
         RETURN
      ENDIF

      NLIST = MIN(NLIST,NLISTGOT)
      DO I = 1,NLIST
         PLIST(I) = DBUF(ILOC + I) 
      ENDDO
    
      IRTFLG = 0
      RETURN

      END

C     ------------------------- LUNDOCREDNXT -----------------------------
C
C    LUNDOCREDNXT(LUNDOC,IKEY,DLIST,NMAXDL,UNUSED,ICOUNT,IRTFLG)
C
C    PURPOSE:    SUBROUTINE TO RECOVER NEXT KEY LINE FROM A DOC FILE
C
C    PARAMETERS:  LUNDOC  IO UNIT                                (SENT)
C                 IKEY    NUMBER OF KEY                          (RET.)
C                 DLIST   ARRAY CONTAINING NUMBERS               (RET.)
C                 NMAXDL  MAX DLIST ARRAY DIMENSION              (SENT)
C                 UNUSED  UNUSED                                 (---)
C                 ICOUNT  NUMBER OF ELEMENTS RETURNED IN DLIST   (RET.)
C                 IRTFLG  1=ERROR                                (RET.)
C                         2=EOF                                  (RET.)
C
C   CALLER:       UNSDAL
C
C--*********************************************************************

      SUBROUTINE LUNDOCREDNXT(LUNDOC,IKEY,DLIST,NMAXDL,
     &                        UNUSED,ICOUNT,IRTFLG)

      USE DOCIC_INFO

      INCLUDE 'CMBLOCK.INC'

      REAL, DIMENSION(*)          :: DLIST
      REAL, DIMENSION(:), POINTER :: IPQ
      LOGICAL                     :: WARNIT,WANTERRT

      ICOUNT = 0
      IKEY   = 0
      IRTFLG = 1

      IF (LUNDOC < 0) THEN
C        INCORE DOC FILE
         NIC    =  -LUNDOC
         IPQ    => LOCDOC(NIC)%IPT 
         MAXY   =  NUMKEYS(NIC)
         MAXXIC =  NUMCOLS(NIC) 
         IKEY   =  NEXTKEY(NIC)
         !write(6,*) ' lundoc; looking for: ',ikey

         DO WHILE (IKEY > 0 .AND. IKEY <= MAXY) 
            ILOC   = (IKEY - 1) * MAXXIC + 1
C           FIND NUMBER OF REGS. FOR THIS KEY
            ICOUNT = IPQ(ILOC)

            IF (ICOUNT > 0) THEN
C              FOUND FILLED KEY, MAKE SURE DLIST DOES NOT OVERFLOW
               ICOUNT = MIN(ICOUNT,NMAXDL)
C              FILL DLIST
               DO I= 1,ICOUNT
                  DLIST(I) = IPQ(ILOC+I)
               ENDDO
               NEXTKEY(NIC) = IKEY + 1  ! NEXT KEY TO BE CHECKED
               IRTFLG       = 0
               RETURN
            ENDIF
            IKEY = IKEY + 1
         ENDDO

         NEXTKEY(NIC) = 1         ! FOR SETTING IKEY ON NEXT ACCESS
         !write(6,*) ' lundoc; set nextkey(',nic,'): 1'
         IRTFLG       = 2
         RETURN
      ENDIF

C     REGULAR DOC FILE
      WARNIT   = .FALSE.
      WANTERRT = .FALSE.
      DO
C        READ NEXT LINE FROM DOC FILE
         CALL LUNDOCREDLIN(LUNDOC,WANTERRT,WARNIT,.FALSE.,
     &                     DLIST,NMAXDL+1,0,IKEY,ICOUNT,IRTFLG)

         IF (IRTFLG < 0) THEN
C           END OF FILE, RETURN 2=NOT FOUND
            IRTFLG = 2
            RETURN

         ELSEIF (IRTFLG > 0) THEN
C           ERROR RETRIEVING PLIST FOR NEXT KEY
            RETURN
         ENDIF

         IF (ICOUNT > 0) THEN
C           HAVE RETRIEVED PLIST FOR NEXT KEY
            IRTFLG = 0
            RETURN
         ENDIF
      ENDDO

      END



C     ------------------------- LUNDOCINFO -----------------------------
C
C  LUNDOCINFO(NDOC,MAXKEYT,MAXREGT,KEYSINUSE,SAYIT,IRTFLG)
C
C  PURPOSE:  DETERMINES MAXKEY AND MAXREG INSIDE A DOCUMENT FILE. 
C
C  PARAMETERS:
C         NDOC          LOGICAL UNIT FOR DOC FILE                 (SENT)
C         MAXKEYT       NUMBER OF HIGHEST KEY                 (RETURNED)
C         MAXREGT       MAX. NUMBER OF REGISTERS PER LINE     (RETURNED)
C                       (THIS IS 1 LESS THAN COLS. NEEDED FOR
C                       RECOVERY WITH UNSDAL, ETC.)
C         KEYSINUSE     NUMBER OF USED KEYS                   (RETURNED)
C         SAYIT         FLAG TO ECHO NUMBERS                  (SENT)
C         IRTFLG        ERROR FLAG (O IS NORMAL RETURN)       (RETURNED)
C     
C--*********************************************************************

        SUBROUTINE LUNDOCINFO(NDOC,MAXKEYT,MAXREGT,KEYSINUSE,
     &                         SAYIT,IRTFLG)

        USE DOCIC_INFO

        INCLUDE 'CMBLOCK.INC' 

        INTEGER             :: NDOC,MAXKEYT,MAXREGT,KEYSINUSE,IRTFLG
        LOGICAL             :: SAYIT

        REAL, POINTER       :: IPQ(:)
        CHARACTER(LEN=80)   :: RECLIN  ! ONLY NEEDS START
        INTEGER             :: ICOMM,MYPID,MPIERR,NLET,IC,I

        INTEGER             :: lnblnk

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        IRTFLG = 1

        IF (NDOC < 0) THEN
C         INCORE DOC. FILE ----------------------------------------

C         GET ARRAY SIZE FOR INCORE FILE (FIXED WHEN IT WAS CREATED)
          IC   = - NDOC
          IF (IC > MAXICDOCS) THEN
C            FILE LIST INDEX OUT OF RANGE
             CALL ERRT(102,'MAX. INCORE DOC. FILE. NUMBER',MAXICDOCS)
             RETURN
          ENDIF

C         FIND NUMBER OF REGISTERS 
          MAXCOLS = NUMCOLS(IC) 
          MAXREGT = MAXCOLS - 1
          MAXKEYS = NUMKEYS(IC)

C         FIND HIGHEST KEY NUMBER IN USE
          IPQ  => LOCDOC(IC)%IPT
 
          DO I = 1,MAXKEYS
             ILOC   = (I - 1) * MAXCOLS + 1
             ICOUNT = IPQ(ILOC)
             IF (ICOUNT > 0) THEN
                 MAXKEYT   = I
                 MAXREGT   = ICOUNT
                 KEYSINUSE = I
             ENDIF
          ENDDO
          IRTFLG = 0
          RETURN
        ENDIF

C       REGULAR DOC. FILE ----------------------------------------

C       WANT TO FIND MAXIMUM KEY & REGISTER NUMBER IN USE 

        MAXKEYT   = 0
        MAXREGT   = 0
        KEYSINUSE = 0

        DO WHILE (.TRUE.)

#ifdef USE_MPI
           IF (MYPID <= 0) THEN
              READ(NDOC,'(A80)',END=100,IOSTAT=IER) RECLIN

              IF (RECLIN(2:2) == ';' .OR. 
     &            RECLIN(1:2) == '; ' ) THEN
C                COMMENT LINE
                 CYCLE
              ENDIF

              BACKSPACE(NDOC)   ! REREAD LINE
              READ(NDOC,*,END=100,IOSTAT=IER) RECLIN

              ENDFILE = .TRUE.    
              READ(NDOC,*,END=100,IOSTAT=IER) NKEY,NREGPLINE
              ENDFILE = .FALSE.
           ENDIF

            
 100       CALL BCAST_MPI('LUNDOCINFO','ENDFILE',ENDFILE,1,'L',ICOMM)
           CALL BCAST_MPI('LUNDOCINFO','IER',IER,1,'I',ICOMM)
           IF (ENDFILE) GOTO 799

           CALL BCAST_MPI('LUNDOCINFO','NKEY',NKEY,1,'I',ICOMM)
           CALL BCAST_MPI('LUNDOCINFO','NREGPLINE',NREGPLINE,1,
     &                    'I',ICOMM)
#else
           READ(NDOC,'(A80)',END=799,IOSTAT=IER) RECLIN
           IF (RECLIN(2:2) == ';' .OR. 
     &         RECLIN(1:2) == '; ') THEN
C             COMMENT LINE
              CYCLE
           ENDIF

           NLET = lnblnk(RECLIN)
           IF (NLET <= 0) THEN
C              BLANK LINE
               WRITE(NOUT,92) 
               CYCLE
           ENDIF

           BACKSPACE(NDOC)   ! REREAD LINE
           READ(NDOC,*,END=799,IOSTAT=IER) NKEY,NREGPLINE
           !write(6,*) 'ier,nkey,nregpline:',ier, nkey,nregpline,maxregt 
           
#endif

           IF (IER < 0) THEN
C             EOF ON READ
              EXIT

           ELSEIF (IER > 0) THEN
C             ERROR ON READ, TRY OLD DOC. FILE FORMAT
              BACKSPACE(NDOC)

              IF (MYPID <= 0) THEN
                 READ(NDOC,83,IOSTAT=IER) NKEY,NREGPLINE
              ENDIF
83            FORMAT(I6,I1)

#ifdef USE_MPI
              CALL BCAST_MPI('LUNDOCINFO','NKEY',NKEY,1,'I',ICOMM)
              CALL BCAST_MPI('LUNDOCINFO','NREGPLINE',NREGPLINE,1,
     &                       'I',ICOMM)
              CALL BCAST_MPI('LUNDOCINFO','IER',IER,1,'I',ICOMM)
#endif
           ENDIF

92        FORMAT('  WARNING; SKIPPING EMPTY LINE IN DOC FILE')
          IF (IER == 0) THEN
C             NOT A COMMENT LINE AND READS KEY & ICOUNT OK
              IF (NKEY < 0) THEN
                 WRITE(NOUT,90) 
90               FORMAT('  WARNING; CONTINUATION LINE IN DOC FILE')

              ELSEIF (NKEY == 0)THEN
                 WRITE(NOUT,91) 
91               FORMAT('  *** SKIPPING ILLEGAL KEY:0  IN DOC FILE')

              ELSEIF (NKEY > 0 .AND.  NREGPLINE == 0) THEN
                 WRITE(NOUT,92) 

              ELSE
C                REGULAR REGISTER LINE
                 IF (NKEY      > MAXKEYT) MAXKEYT = NKEY
                 IF (NREGPLINE > MAXREGT) MAXREGT = NREGPLINE
                 KEYSINUSE = KEYSINUSE + 1
              ENDIF

           ENDIF
        ENDDO

799     IF (SAYIT .AND. MYPID <= 0) then
              WRITE(NOUT,97) MAXREGT,KEYSINUSE,MAXKEYT
97         FORMAT('  Doc file has:',I4,'  registers and:',I6,
     &            ' keys,  Highest key in use:',I6,/)
        ENDIF

        IF (MYPID <= 0) REWIND(NDOC)

        IRTFLG = 0

        END



C     ------------------------- LUNDOCGETCOM -----------------------------
C
C    LUNDOCGETCOM(LUNDOC,IKEY,PLIST,NLIST,TILLEND,IRTFLG)
C
C    PURPOSE:    GET A SPECIFIED COMMENT KEY FROM FILE
C
C    PARAMETERS: NDOC    IO UNIT                                (SENT)
C                IKEY    COMMENT KEY WANTED (<0)                (SENT)
C                ILIST   ARRAY CONTAINING VALUES                (RET.)
C                NLIST   NUMBER OF ELEMENTS IN ARRAY       (SENT/RET.)
C                        MAX. MUST BE PROVIDED FROM CALLER
C                TILLEND FLAG TO KEEP READING TILL END          (SENT)
C                IRTFLG                                         (RET.)
C
C--*********************************************************************

        SUBROUTINE LUNDOCGETCOM(LUNDOC,IKEY,PLIST,NLIST,TILLEND,IRTFLG)

        INCLUDE 'CMBLOCK.INC' 

        REAL               :: PLIST(*)
        CHARACTER(LEN=180) :: RECLIN
        LOGICAL            :: TILLEND

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       SET ERROR RETURN FLAG
        IRTFLG = 1

        IKEYA  = ABS(IKEY)
        IF (NLIST <= 0) THEN
           CALL ERRT(101,'*** NUMBER OF REGISTERS NOT SPECIFIED',NE)
           RETURN

        ELSEIF (IKEY > 0) THEN
           CALL ERRT(101,'*** DID NOT REQUEST COMMENT KEY',NE)
           RETURN

        ELSEIF (IKEY > 9999999) THEN
           CALL ERRT(102,'ILLEGAL COMMENT KEY NUMBER',IKEY)
           RETURN
        ENDIF

        ICOUNT = 0

C -------------------------------------------------------

C       READ NEXT LINE FROM DOC FILE
#ifdef USE_MPI
10      CONTINUE   
        IF (MYPID == 0) THEN
           ERRFILE = .TRUE.
           ENDFILE = .TRUE.
           READ (LUNDOC,81,ERR=100,END=200) RECLIN
81         FORMAT(A120)
           ERRFILE = .FALSE.
           ENDFILE = .FALSE. 
        ENDIF

100     CALL BCAST_MPI('LUNDOCGETCOM','ERRFILE',ERRFILE,1,'L',ICOMM)
        IF (ERRFILE) GOTO 998

200     CALL BCAST_MPI('LUNDOCGETCOM','ENDFILE',ENDFILE,1,'L',ICOMM)
        IF (ENDFILE) GOTO 997

        CALL BCAST_MPI('LUNDOCGETCOM','RECLIN',RECLIN,180,'C',ICOMM)

#else
10      READ (LUNDOC,81,ERR=998,END=997) RECLIN
81      FORMAT(A)
#endif

        IF (RECLIN(2:2) == ';' .OR. RECLIN(1:2) == '; ' )THEN
C          COMMENT LINE FOUND, CHECK FOR COMMENTED KEY WITH ERR
           READ(RECLIN(3:),*,IOSTAT=IERR) KEYGOT,ICOUNTT
           IF (IERR .NE. 0) GOTO 10        ! NOT A COMMENT KEY LINE

C          MAKE SURE PLIST DOES NOT OVERFLOW
           ICOUNTT = MIN(ICOUNTT,NLIST)
           KEYGOTA = ABS(KEYGOT)

           IF (ICOUNTT > 0 .AND. KEYGOTA == IKEYA) THEN
C             READ THE COMMENT KEY DATA INTO PLIST
CCCC              READ(RECLIN(3:),*,ERR=50,END=998)IKEYDUM,ICOUNTDUM,
              READ(RECLIN(3:),*,ERR=50)IKEYDUM,ICOUNTDUM,
     &                                (PLIST(K),K=1,ICOUNTT)
              ICOUNT = ICOUNTT
           ENDIF

           IF (KEYGOTA == IKEYA .AND. .NOT. TILLEND) GOTO 997
        ENDIF

C       READ NEXT LINE OF DOC FILE
50      GOTO 10

C ---------------------------------------------------------------

997     IF (ICOUNT <= 0) THEN
C          END OF DOCUMENT FILE FOUND WITHOUT COMMENT KEY
           CALL ERRT(102,'COMMENT KEY NOT FOUND',IKEYA)
           RETURN
        
        ELSEIF (ICOUNT < NLIST) THEN
C          UNDERFLOW

           WRITE(NOUT,90) NLIST,ICOUNT
90         FORMAT('  *** WANTED: ',I3,' REGISTERS BUT ONLY GOT: ',I4/)
           CALL ERRT(100,'LUNDOCGETCOM',NE)
           NLIST  = ICOUNT
           RETURN

        ENDIF
        NLIST  = ICOUNT
        IRTFLG = 0

998     RETURN
        END


C     ------------------------- LUNDOCPUTCOM -----------------------------
C
C    LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)
C
C    PURPOSE:    PUT A TEXT COMMENT IN DOC. FILE (NOT A COMMENT KEY)
C
C    PARAMETERS: LUNDOC    DOC FILE I/O UNIT                      (SENT)
C                COMMENT   TEXT COMMENT                           (SENT)
C                IRTFLG    ERROR RETURN FLAG                      (RET.)
C
C--*********************************************************************

        SUBROUTINE LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)

        CHARACTER(LEN=*)  :: COMMENT

        INTEGER           :: LUNDOC,IRTFLG,ICOMM,MYPID,MPIERR,NC

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        IF (LUNDOC <= 0) RETURN

        IF (MYPID <= 0) THEN
           NC = lnblnkn(COMMENT)

           WRITE(LUNDOC,94) COMMENT(1:NC)
94         FORMAT(' ; ',A)
        ENDIF

        END


C     ------------------------- LUNDOCSAYHDR -----------------------------
C
C    LUNDOCSAYHDR(LUNDOC,LUNPUT,IRTFLG)
C
C    PURPOSE:    SUBROUTINE TO ECHO FIRST ENT LINE FROM A DOC FILE
C
C    PARAMETERS:   LUNDOC  DOC FILE IO UNIT                      (SENT)
C    PARAMETERS:   LUNPUT  OUTPUT IO UNIT                        (SENT)
C                  IRTFLG                                        (RET.)
C
C--*********************************************************************

      SUBROUTINE LUNDOCSAYHDR(LUNDOC,LUNPUT,IRTFLG)

        CHARACTER *120 RECLIN
#ifdef USE_MPI
        include 'mpif.h'
        LOGICAL ENDFILE, ERRFILE
        ICOMM   = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)
#else
        MYPID = -1
#endif
        IRTFLG = 1
        IF (LUNDOC <= 0) RETURN

C -------------------------------------------------------
        
C       READ NEXT LINE FROM DOC FILE
#ifdef USE_MPI
10      IF (MYPID == 0) THEN
           ERRFILE = .TRUE.
           ENDFILE = .TRUE.
           READ (LUNDOC,81,ERR=100,END=100) RECLIN
           ERRFILE = .FALSE.
           ENDFILE = .FALSE.
        ENDIF

100     CALL BCAST_MPI('LUNDOCSAYHDR','ENDFILE',ENDFILE,1,'L',ICOMM)
        CALL BCAST_MPI('LUNDOCSAYHDR','ERRFILE',ERRFILE,1,'L',ICOMM)
        IF (ENDFILE .OR. ERRFILE) GOTO 999

        CALL BCAST_MPI('LUNDOCSAYHDR','RECLIN',RECLIN,120,'C',ICOMM)
#else
10      READ (LUNDOC,81,ERR=999,END=999) RECLIN
#endif
81      FORMAT(A120)

        IF (RECLIN(2:2) == ';' .OR. RECLIN(1:2) == '; ')THEN
C          COMMENT LINE FOUND
           ILEN = LNBLNKN(RECLIN)
           IF (MYPID <= 0) WRITE(LUNPUT,*) RECLIN(1:ILEN)
           IRTFLG = 0
           GOTO 999
        ENDIF

C       READ NEXT LINE OF DOC FILE, UTIL WE GET COMMENT
        GOTO 10

C ---------------------------------------------------------------

999     IF (MYPID <= 0) REWIND(LUNDOC)
        RETURN
        END


