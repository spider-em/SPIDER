
C ++********************************************************************
C
C DOCS1 NEW                                      JUN  1999 ARDEAN LEITH
C       ADDED 'DOC RAN'                          AUG  1999 ARDEAN LEITH
C       USED LUNDOCWRTDAT                        AUG  1999 ARDEAN LEITH
C       ADDED 'DOC AND'                          SEPT 1999 ARDEAN LEITH
C       ADDED 'DOC SPLIT'                        OCT  1999 ARDEAN LEITH
C       OPENDOC PARAMETERS                       DEC  2000 ARDEAN LEITH
C       'DOC TOMINESET' ADDED                    JUN  2001 ARDEAN LEITH
C       CLOSED NDOCOUT IN ROUTINES               JUL  2001 ARDEAN LEITH
C       'DOC COM' ADDED                          DEC  2001 ARDEAN LEITH
C       'DOC COM' BUG                            MAY  2002 ARDEAN LEITH
C       'DOC COM' FILE NAME BUG                  SEP  2002 ARDEAN LEITH
C       'DOC COM' MAXY BUF                       JUN  2003 ARDEAN LEITH
C       INCORE OPENDOC                           JUL  2003 ARDEAN LEITH
C       MPI                                      OCT  2003 CHAO YANG
C       'DOC RAN' BUG                            JAN  2004 ARDEAN LEITH
C       'DOC OLD' ADDED                          FEB  2004 ARDEAN LEITH
C       'DOC AND' BUG                            FEB  2004 ARDEAN LEITH
C       'DOC MIR' KEYCOL BUG                     OCT  2004 ARDEAN LEITH
C       'DOC KEY' ADDED                          JUL  2005 ARDEAN LEITH
C       'DOC BOOT' ADDED                         JAN  2006 ARDEAN LEITH
C       'DOC SORT' REVERSE ORDER                 OCT  2010 ARDEAN LEITH
C       'DOC TOMINESET' REMOVED                  NOV  2010 ARDEAN LEITH
C       'DOC SORT A' APPEND OK                   NOV  2010 ARDEAN LEITH
C       'DOC STAT'                               NOV  2010 ARDEAN LEITH
C       'DOC ME' EMPTY DOC FILE BUG              DEC  2010 ARDEAN LEITH
C        other ops still need empty doc file traps!
C        'ME') KEYCOL = 1 BUG                    JAN  2010 ARDEAN LEITH 
C        MOVED 'SUB' and 'AND' TO DOCSUB         JAN  2010 ARDEAN LEITH 
C        MPI ERROR IN  DOCCREATE                 MAR  2011 ARDEAN LEITH 
C        UNIQUE IN 'DOC SORT'                    APR  2012 ARDEAN LEITH 
C        UNIQUE IN 'AT IT' BUG                   SEP  2013 ARDEAN LEITH 
C       'DOC SEP' ADDED                          MAR  2014 ARDEAN LEITH
C       'DOC PLOT' ADDED                         DEC  2014 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  DOCS1(MAXDIM)
C                                                                 
C  PURPOSE: MANIPULATES DOCUMENT FILES.
C   
C  NOTE: TSWITCH SAYS THIS IS A 2 LETTER OP
C                                                                    
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE DOCS1(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        INTEGER, INTENT(IN)          :: MAXDIM
 	CHARACTER(LEN=MAXNAM)        :: DOCNAM
 	CHARACTER(LEN=1)             :: NULL = CHAR(0)
        REAL, ALLOCATABLE            :: DLIST(:)

        INCLUDE 'F90ALLOC.INC'
        REAL,  POINTER               :: DOCBUF(:,:)

        INTEGER, PARAMETER           :: NDOCIN  = 70
        INTEGER, PARAMETER           :: LUNPLOT = 80

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        IF (FCHAR(4:5) == 'CR') THEN
C          CREATE NEW DOC FILE ------------------------- 'DOC CREATE'

           CALL DOCCREATE(MAXDIM)
           RETURN

        ELSEIF (FCHAR(4:5) == 'SH') THEN
C          SHUFFLE OLD DOC FILE ------------------------- DOC SHUFFLE'
C          (ALSO CALLED "SD SHUFFLE")
           CALL SHUFFLEDOC(MAXDIM)
           RETURN

        ELSEIF (FCHAR(4:5) == 'OL') THEN
C          COPY TO  OLD DOC FILE ------------------------ DOC COPY'
           CALL DOCDOWN()
           RETURN

        ELSEIF(FCHAR(4:5) == 'CO') THEN
C          COMBINE A SERIES OF INPUT DOC. FILES --------- DOC COMBINE'
           CALL DOCCOMBINE()
           RETURN

        ELSEIF(FCHAR(4:5) == 'SE') THEN
C          SEPARATE A SERIES OF INPUT DOC. FILES -------- DOC SEPARATE'
           CALL DOCSEPARATE()
           RETURN

        ELSEIF(FCHAR(4:5) == 'PL') THEN
C          GNUPLOT REGISTER CONTENTS  -------------------- 'DOC PLOT'

           CALL DPROFD_GPL(NDOCIN,LUNPLOT)
           RETURN
        ENDIF

C       OPEN EXISTING DOC FILE
C       MAXX IS 1 + NUM OF REGISTERS SINCE DOCBUF CONTAINS KEY ALSO
        MAXX    = 0
        MAXY    = 0
        NDOCINT = NDOCIN

C       DOC RENUMBER NEEDS SEQUENTIAL READ OF LINES NOT BY KEY
        IF (FCHAR(4:5) == 'RE') NDOCINT = -NDOCIN

        CALL GETDOCDAT('INPUT DOCUMENT',.TRUE.,DOCNAM,
     &                 NDOCINT,.TRUE.,MAXX, MAXY,DOCBUF,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        NLIST = MAX(MAXX,1)
        ALLOCATE(DLIST(NLIST), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'DOCS1; DLIST',NLIST)
           GOTO 9995
        ENDIF

        IF (FCHAR(4:5) == 'RA') THEN
C          CREATE RANDOM SELECTED DOC FILE ------------------ 'DOC RAN'
           CALL DOCRAN(MAXX, MAXY, DOCBUF(1,1), DLIST)

        ELSEIF (FCHAR(4:5) == 'BO') THEN
C          CREATE RANDOM SELECTED BOOTSTRAP DOC FILE -------- 'DOC BOOT'
           CALL DOCBOOT(MAXX, MAXY, DOCBUF(1,1), DLIST)

        ELSEIF (FCHAR(4:5) == 'RE') THEN
C          RENUMBER DOC FILE --------------------------------- 'DOC RE'
           CALL DOCRENUMBER(MAXX, MAXY, DOCBUF(1,1), DLIST)

        ELSEIF (FCHAR(4:5) == 'KE') THEN
C          REKEY DOC FILE ----------------------------------- 'DOC KEY'
           CALL DOCREKEY(MAXX, MAXY, DOCBUF(1,1), DLIST)

        ELSEIF (FCHAR(4:5) == 'AN' .OR.
     &          FCHAR(4:5) == 'SU') THEN
C          SUBTRACT DOC FILE ------------------------------- 'DOC SUB'
C          AND DOC FILE      ------------------------------- 'DOC AND'
           CALL DOCSUB(MAXX, MAXY, DOCBUF(1,1), DLIST,NLIST)

        ELSE
C          OTHER DOC FILE OPS--------------------------------- 'DOC ??'
C          SINCE DLIST MAY VARY IT IS NOT USED HERE
           CALL DOCSDO(MAXX, MAXY, DOCBUF(1,1))

        ENDIF

C       DEALLOCATE DOC. FILE MEMORY
9995    IF (ASSOCIATED(DOCBUF)) DEALLOCATE(DOCBUF)

C       DEALLOCATE DLIST MEMORY
        IF (ALLOCATED(DLIST))  DEALLOCATE(DLIST)

        RETURN
        END


C       --------------------- DOCSDO ----------------------------------

C       SORT THE INPUT DOC FILE-------------------------------- 'AT IT'
C       SORT THE INPUT DOC FILE----------------------------- 'DOC SORT'
C       MIRROR THE INPUT DOC FILE------------------------- 'DOC MIRROR'
C       MERGE THE TWO INPUT DOC FILES---------------------- 'DOC MERGE'
C       STATISTICS FROM DOC FILES -------------------------- 'DOC STAT'

	SUBROUTINE DOCSDO(MAXX, MAXY, DOCBUF)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

	CHARACTER(LEN=MAXNAM)        :: DOCNAM2,DOCNAM3
	CHARACTER(LEN=1)             :: NULL = CHAR(0)
	CHARACTER(LEN=MAXNAM)        :: ANSW
        LOGICAL                      :: NEWFILE,ERRI2,RENUMBER

        REAL                         :: DOCBUF(MAXX*MAXY)
        REAL,    ALLOCATABLE         :: SORTED(:),SORTED2(:)
        REAL,    ALLOCATABLE         :: DLIST(:)
        INTEGER, ALLOCATABLE         :: KEYLIST(:)

        INCLUDE 'F90ALLOC.INC'
        REAL, POINTER                :: DOCBUF2(:,:)
	CHARACTER(LEN=80)            :: PROMPT
        LOGICAL                      :: SENDIT,REP_KEY,UNIQUE
        LOGICAL                      :: REVERSE
        INTEGER                      :: NEEDINC
        LOGICAL                      :: APPENDOK,MESSAGE

        INTEGER, PARAMETER           :: NDOCIN2  = 71
        INTEGER, PARAMETER           :: NDOCOUT  = 72
        INTEGER, PARAMETER           :: NDOCOUT2 = 73

        IF (FCHAR(4:5) == 'SO') THEN
C                    12345678901234567890123456789012345678901234567890123456789
           PROMPT  ='COLUMN TO BE SORTED BY (0 IS KEY) (<0 TO REVERSE)' 
           NLETP    = 49
           RENUMBER = .TRUE.

        ELSEIF (FCHAR(4:5) == 'IT') THEN
           REVERSE = .FALSE.

        ELSEIF (FCHAR(4:5) == 'MI') THEN
           PROMPT   = 'COLUMN TO BE MIRRORED (0 IS KEY)' 
           NLETP    = 32
           RENUMBER = .FALSE.

        ELSEIF (FCHAR(4:5) == 'ME') THEN
           PROMPT   = 'COLUMN TO BE MERGED BY (0 IS KEY)'
           NLETP    = 35
           RENUMBER = .FALSE.

C          MERGE USES 2 INPUT DOC. FILES
C          MAXX2 IS 1 + NUM OF REGISTERS SINCE DOCBUF CONTAINS KEY ALSO
           MAXX2  = 0
           MAXY2  = 0
           CALL GETDOCDAT('SECOND INPUT DOCUMENT',.TRUE.,DOCNAM2,
     &                  NDOCIN2,.TRUE.,MAXX2, MAXY2,DOCBUF2,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

         ELSEIF (FCHAR(4:5) == 'ST') THEN
           PROMPT   = 'COLUMN TO BE ANALYZED (0 IS KEY)' 
           NLETP    = 40
        ENDIF

        IF (FCHAR(4:5) .NE. 'ST') THEN
C          OPEN OUTPUT DOCUMENT FILE
           CALL FILERD(DOCNAM3,NLET,NULL,'OUTPUT DOCUMENT',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9997

           APPENDOK = (FCHAR(4:9) == 'SORT A')
           MESSAGE  = .NOT. APPENDOK
           CALL OPENDOC(DOCNAM3,.TRUE.,NLET,NDOCOUT,NICDOCOUT,.FALSE.,
     &                  ' ',.FALSE.,APPENDOK,MESSAGE,NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9997
        ENDIF

        IF (FCHAR(4:5) == 'SO' .OR. 
     &      FCHAR(4:5) == 'MI' .OR.
     &      FCHAR(4:5) == 'ME') THEN
C          SORT THE INPUT DOC FILE------------------------- 'DOC SORT'
C          MIRROR THE INPUT DOC FILE----------------------- 'DOC MIRROR'
C          MERGE THE TWO INPUT DOC FILES------------------- 'DOC MERGE'

11         CALL RDPRI1S(KEYCOL,NOT_USED,PROMPT(1:NLETP),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9997

           REVERSE = (KEYCOL < 0)
           IF (FCHAR(4:5) .NE. 'ME') KEYCOL = ABS(KEYCOL)

           IF (FCHAR(4:5) == 'ME' .AND. MAXX > 0 )    THEN
              IF (ERRI2(KEYCOL,IDUM,1,-1,MAXX-1,0,0))   GOTO 11
           ELSEIF (MAXX > 0 .AND. MAXY > 0)             THEN
              IF (ERRI2(KEYCOL,IDUM,1, 0,MAXX-1,0,0))   GOTO 11
           ENDIF

           IF (FCHAR(4:5) == 'SO') THEN
              CALL RDPRMC(ANSW,NLET,.TRUE.,
     &       'COMPRESS & RENUMBER KEYS? (Y/N), REMOVE DUPLICATES (Y/N)',
     &       NULL,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9997
              CALL SSUPCAS(ANSW)
              RENUMBER = (ANSW(1:1) .NE. 'N')
              UNIQUE   = ( NLET > 1 .AND. INDEX(ANSW(2:),'Y') > 0)

           ELSEIF (RENUMBER) THEN
              CALL RDPRMC(ANSW,NLET,.TRUE.,
     &                 'COMPRESS & RENUMBER KEYS? (Y/N)',NULL,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9997
              CALL SSUPCAS(ANSW)
              RENUMBER = (ANSW .NE. 'N')
           ENDIF

        ELSEIF (FCHAR(4:5) == 'IT') THEN
           KEYCOL   = 1
           RENUMBER = .TRUE.
           UNIQUE   = .TRUE.

        ELSEIF (FCHAR(4:5) == 'SP') THEN
           CALL FILERD(DOCNAM2,NLET,NULL,'SECOND OUTPUT DOCUMENT',IRT)
           IF (IRT .NE. 0) RETURN

           CALL OPENDOC(DOCNAM2,.TRUE.,NLET,NDOCOUT2,NICDOCOUT2,.FALSE.,
     &                  ' ',.FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

         IF (FCHAR(4:5) == 'SO' .AND. KEYCOL == 0) THEN
C          NO NEED TO SORT LIST ----------------------------- 'DOC SORT'
           NVAL   = MAXX - 1
           NEWKEY = 0

           DO IROW = 1,MAXY
               ILOC   = (IROW - 1) * MAXX + 1
               ICOUNT = DOCBUF(ILOC)
               IF (ICOUNT > 0) THEN
                  IF (RENUMBER) THEN
C                    RENUMBER THE KEYS
                     NEWKEY = NEWKEY + 1
                  ELSE
                     NEWKEY = IROW
                  ENDIF
C                 PUSH VALUES INTO OUTPUT DOC. FILE
                  CALL LUNDOCWRTDAT(NICDOCOUT,NEWKEY,
     &                           DOCBUF(ILOC+1),NVAL,IRTFLG)
               ENDIF
           ENDDO
           GOTO 9990

        ELSEIF ((FCHAR(4:5) == 'ME') .AND. 
     &          ((KEYCOL .LT. 0) .OR.
     &           (MAXX  < 1 .OR. MAXY  < 1) .OR. 
     &           (MAXX2 < 1 .OR. MAXY2 < 1))) THEN
C          NO NEED TO SORT LISTS -------------------------- 'DOC MERGE'

           IF ((MAXX  < 1 .OR. MAXY  < 1) .AND. 
     &         (MAXX2 < 1 .OR. MAXY2 < 1)) THEN
                 ! OK   CALL ERRT(102,'EMPTY DOC. FILES',NDUM)
                 GOTO 9990
           ENDIF

           MAXXT = MAX(MAXX,MAXX2)

           ALLOCATE(DLIST(MAXXT+1),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'DLIST',MAXXT)
              GOTO 9990
           ENDIF

C          MERGING ALL KEYS
           KEYNEW = 0

C          MAKE SURE ALL REGISTERS ARE ZEROED IN OUTPUT
           NVAL = MAX(MAXX,MAXX2) - 1
           DO IREG = 1,NVAL
              DLIST(IREG) = 0.0
           ENDDO

           IF (MAXX > 0 .AND. MAXY > 0) THEN 
C             COPY VALUES FROM FIRST FILE
              DO IKEYT = 1,MAXY
                 ICOUNT = DOCBUF((IKEYT - 1) * MAXX + 1)
                 IF (ICOUNT > 0 ) THEN
C                   KEY EXISTS PUSH DLIST FROM FILE 1 INTO OUTPUT FILE
                    DO IREG = 1,MAXX
                       DLIST(IREG) = DOCBUF((IKEYT - 1) * MAXX + IREG)
                    ENDDO

C                   PUSH DLIST INTO DOC. FILE
                    KEYNEW = KEYNEW + 1
                    CALL LUNDOCWRTDAT(NICDOCOUT,KEYNEW,DLIST(2),
     &                                NVAL,IRTFLG)
                 ENDIF
              ENDDO
           ENDIF

           IF (MAXX2 > 0 .AND. MAXY2 > 0) THEN 
C             COPY VALUES FROM SECOND FILE
              DO IKEYT = 1,MAXY2
                 ICOUNT  = DOCBUF2(1,IKEYT)
                 IF (ICOUNT > 0) THEN
C                   KEY EXISTS, PUSH DLIST FROM FILE 2 INTO OUTPUT FILE
                    DO IREG = 1,MAXX2
                       DLIST(IREG) = DOCBUF2(IREG,IKEYT)
                    ENDDO

C                   PUSH DLIST INTO DOC. FILE
                    KEYNEW = KEYNEW + 1
                    CALL LUNDOCWRTDAT(NICDOCOUT,KEYNEW,DLIST(2),
     &                                NVAL,IRTFLG)
                 ENDIF
              ENDDO
           ENDIF
           GOTO 9990
     
        ELSEIF (FCHAR(4:5) == 'SO' .OR.
     &          FCHAR(4:5) == 'MI' .OR.
     &          FCHAR(4:5) == 'ME' .OR.
     &          FCHAR(4:5) == 'IT') THEN
C          NEED ONE OR MORE SORTED LISTS

C          SPLIT THE INPUT DOC FILE------------------------- 'DOC SPLIT'
C          SORT THE INPUT DOC FILE------------------------------ 'AT IT'
C          SORT THE INPUT DOC FILE--------------------------- 'DOC SORT'
C          MIRROR THE INPUT DOC FILE----------------------- 'DOC MIRROR'
C          MERGE THE TWO INPUT DOC FILES-------------------- 'DOC MERGE'

           ALLOCATE(SORTED(MAXY),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'SORTED IN DOCS1',MAXY)
              GOTO 9997
           ENDIF

C          SORTED RETURNS LIST OF KEYS SORTED BY USING VALUE IN KEYCOL
           CALL SORTIT(DOCBUF,MAXX,MAXY,KEYCOL,SORTED,
     &                 IKEYS,.TRUE.,UNIQUE,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9995

           IF (FCHAR(4:5) == 'ME') THEN
C             SORT THE SECOND LIST OF KEYS BUT RETURN VALUES IN SORTED2
              ALLOCATE(SORTED2(MAXY2),STAT=IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                 CALL ERRT(46,'SORTED IN DOCS1',MAXY2)
                 GOTO 9995
              ENDIF

C             RETURN SORTED LIST OF VALUES IN KEYCOL
              REP_KEY = (FCHAR(4:5) == 'ME')
              UNIQUE  = .FALSE. 
              CALL SORTIT(DOCBUF2(1,1),MAXX2,MAXY2,KEYCOL,SORTED2,
     &                    IKEYS2,REP_KEY,UNIQUE,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9995
          ENDIF
        ENDIF

        NLIST = MAXX
        IF (FCHAR(4:5) == 'MI') NLIST = MAX(2,KEYCOL+1)
        IF (FCHAR(4:5) == 'ME') NLIST = MAX(MAXX,MAXX2)
 
        ALLOCATE(DLIST(NLIST),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'DLIST IN DOCS1',NLIST)
           GOTO 9995
        ENDIF

        IF (FCHAR(4:5) == 'ST') THEN
C          DOC FILE COL STATISTICS--------------------------- 'DOC STAT'

           CALL RDPRI1S(KEYCOL,NOT_USED,PROMPT(1:NLETP),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9997
 
           IF (KEYCOL > MAXX) THEN
              CALL ERRT(102,'MAX. COL. IN DOC. FILE',MAXX)
              GOTO 9995
           ENDIF

           IREG = KEYCOL + 1    ! DOCBUF COLUMN
           NGOT = 0
           VMIN =  HUGE(VMIN)
           VMAX = -HUGE(VMIN)
           VSUM = 0.0
           VSQ  = 0.0

           DO KEY=1,MAXY
              NREGS = DOCBUF((KEY - 1) * MAXX + 1) 
              !write(6,*) nregs,DOCBUF((KEY - 1) * MAXX + IREG)
 
              IF (NREGS > 0 .AND. KEYCOL <= NREGS) THEN
                 VAL  = DOCBUF((KEY - 1) * MAXX + IREG)
                 VMIN = MIN(VMIN,VAL)
                 VMAX = MAX(VMAX,VAL)
                 VSUM = VSUM + VAL
                 VSQ  = VSQ  + VAL * VAL
                 NGOT = NGOT + 1
              ENDIF
           ENDDO

           IF (NGOT <= 0) THEN
              CALL ERRT(102,'NO VALUES IN COL.',KEYCOL)
              GOTO 9995
           ELSE
              VAVG = VSUM / FLOAT(NGOT)

              VSIG = 0.0
              VTOP  = VSQ - VSUM * VSUM / FLOAT(NGOT)
              IF (NGOT > 1) VSIG = SQRT(VTOP / FLOAT(NGOT - 1))    
           ENDIF

          WRITE(NOUT,90) KEYCOL,NGOT,VMIN,VMAX,VSUM,VAVG,VSIG
90        FORMAT('  COLUMN:', I2,
     &           '  VALUES:', I6,
     &           '  MIN:',    1PG10.3,
     &           '  MAX:',    1PG10.3,
     &           '  SUM:',    1PG12.5,
     &           '  AVG:',    1PG12.5,
     &           '  SIG:',    1PG12.5)

          CALL REG_SET_NSEL(1,5,FLOAT(NGOT),VMIN,VMAX,VSUM,VAVG,IRTFLG)
          CALL REG_SET_NSEL(6,6,VSIG,IRTFL)

        ELSEIF (FCHAR(4:5) == 'SO' .OR. FCHAR(4:5) == 'IT') THEN
C          SORT THE INPUT DOC FILE--------------------------- 'DOC SORT'
C          SORT THE INPUT DOC FILE--------------------------- 'AT IT'  '
           NEWKEY = 0

           I1    = 1
           I2    = IKEYS
           ISTEP = 1
           IF (REVERSE) THEN
C             DESCENDING ORDER
              I1    = IKEYS
              I2    = 1
              ISTEP = -1
           ENDIF

           DO IROW = I1,I2,ISTEP
               IKEY = SORTED(IROW)
               IF (RENUMBER) THEN
C                 RENUMBER THE KEY COLUMN
                  NEWKEY   = NEWKEY + 1
                  IKEYT    = NEWKEY
               ELSE
C                 KEEP ORIGINAL KEYS 
                  IKEYT    = IKEY
               ENDIF
               DO ICOL = 2,MAXX
                  DLIST(ICOL) = DOCBUF((IKEY - 1) * MAXX + ICOL)
               ENDDO

C              PUSH DLIST INTO DOC. FILE
               NVAL  = MAXX - 1
               CALL LUNDOCWRTDAT(NICDOCOUT,IKEYT,DLIST(2),NVAL,IRTFLG)

           ENDDO

        ELSEIF (FCHAR(4:5) == 'SP') THEN
C          SPLIT THE INPUT DOC FILE------------------------- 'DOC SPLIT'
           NVAL    = MAXX - 1
           NEWKEY1 = 0
           NEWKEY2 = 0
           DO IY = 1,MAXY,2
               KEY1 = DOCBUF((IY-1) * MAXX + 1)
               IF (KEY1 > 0) THEN
C                 PUSH KEY INTO FIRST FILE
                  DO ICOL = 2,MAXX
                     DLIST(ICOL) = DOCBUF((IY - 1) * MAXX + ICOL)
                  ENDDO
                  NEWKEY1 = NEWKEY1 + 1
C                 PUSH DLIST INTO FIRST DOC. FILE
                  CALL LUNDOCWRTDAT(NICDOCOUT,NEWKEY1,DLIST(2),
     &                              NVAL,IRTFLG)
               ENDIF

               IF (IY .LT. MAXY) THEN
C                 PUSH NEXT KEY INTO SECOND FILE
                  KEY2 = DOCBUF((IY) * MAXX + 1)
                  IF (KEY2 > 0) THEN
                     DO ICOL = 2,MAXX
                        DLIST(ICOL) = DOCBUF((IY) * MAXX + ICOL)
                     ENDDO
                     NEWKEY2 = NEWKEY2 + 1
                     CALL LUNDOCWRTDAT(NICDOCOUT2,NEWKEY2,DLIST(2),
     &                                 NVAL,IRTFLG)
                 ENDIF
              ENDIF
           ENDDO

        ELSEIF (FCHAR(4:5) == 'MI') THEN
C          'MIRRORING'  VALUES IN COLUMN: KEYCOL --------- 'DOC MIRROR'

           IF (KEYCOL == 0) THEN
C             MIRROR BY KEY EXISTANCE (NOT REGISTER CONTENTS)
              DLIST(2) = 1.0
              DO IKEY = 1,MAXY
                 IGOT = DOCBUF((IKEY - 1) * MAXX + 1)
                 IF (IGOT <= 0) THEN
C                   NONEXISTANT KEY, PUSH 1.0 INTO DOC. FILE
                    CALL LUNDOCWRTDAT(NICDOCOUT,IKEY,
     &                                DLIST(2),1,IRTFLG)
                 ENDIF
              ENDDO
              GOTO 9990
           ENDIF

           NEWKEY  = 0
           IKEY    = SORTED(1)
           LASTVAL = DOCBUF((IKEY - 1) * MAXX + KEYCOL + 1)

C          NO SENSE TO OTHER MISSING DOC COLUMNS ??
           DO I = 1,NLIST
              DLIST(I) = 0.0
           ENDDO

           DO IROW = 1,IKEYS
              IKEY = SORTED(IROW)
              IVAL = DOCBUF((IKEY - 1) * MAXX + KEYCOL + 1)
              IF (IVAL > (LASTVAL-1)) THEN
                 DO IT=LASTVAL+1,IVAL-1
C                   FILL IN MISSING VALUES FROM KEY COLUMN
                    NEWKEY   = NEWKEY + 1

C                   KEYS ARE RENUMBERED IF NOT FILLING FIRST COL.
                    DLIST(KEYCOL+1) = IT

C                   PUSH DLIST INTO DOC. FILE
                    NVAL  = NLIST - 1
                    IKEYT = NEWKEY
                    CALL LUNDOCWRTDAT(NICDOCOUT,IKEYT,DLIST(2),
     &                                NVAL,IRTFLG)
                 ENDDO
              ENDIF
              LASTVAL = IVAL
           ENDDO
  
        ELSEIF (FCHAR(4:5) == 'ME') THEN
C          MERGING  VALUES IN COLUMN: KEYCOL -------------- 'DOC MERGE'

           MAXXT = MAX(MAXX,MAXX2)

           IF (KEYCOL == 0) THEN
C             MERGING BY KEY
              MAXYT = MAX(MAXY,MAXY2)

              DO KEYT = 1,MAXYT
                 ICOUNT1 = 0
                 IF (KEYT <= MAXY) 
     &              ICOUNT1 = DOCBUF((KEYT - 1) * MAXX + 1)
                 ICOUNT2 = 0

                 IF (KEYT <= MAXY2) 
     &               ICOUNT2 = DOCBUF2(1,KEYT)

                 IF (KEYT <= MAXY .AND. ICOUNT1 > 0 .AND. 
     &              ICOUNT2 == 0) THEN
C                   KEY1 EXISTS AND KEY2 DOES NOT EXIST, PUSH DLIST 
C                   FROM FILE 1 INTO OUTPUT FILE
                    DO IREG = 1,MAXX
                       DLIST(IREG) = DOCBUF((KEYT - 1) * MAXX + IREG)
                    ENDDO

C                   PUSH DLIST INTO DOC. FILE
                    NVAL  = MAXX - 1
                    CALL LUNDOCWRTDAT(NICDOCOUT,KEYT,DLIST(2),
     &                                NVAL,IRTFLG)

                 ELSEIF (KEYT <= MAXY2 .AND. ICOUNT2 .NE. 0) THEN
C                   KEY2 EXISTS, PUSH DLIST FROM FILE 2 INTO OUTPUT FILE
                    DO IREG = 1,MAXX2
                       DLIST(IREG) = DOCBUF2(IREG,KEYT)
                    ENDDO

C                   PUSH DLIST INTO DOC. FILE
                    NVAL  = MAXX2 - 1
                    CALL LUNDOCWRTDAT(NICDOCOUT,KEYT,DLIST(2),
     &                                NVAL,IRTFLG)

                 ENDIF
              ENDDO

           ELSE                    
C             MERGING BY COLUMN OTHER THAN KEY
              KEYNEW = 0

C             POINT TO NEXT VALUE IN SORTED LIST FROM FILE 1
              IGO1     = 1
              KEY1     = SORTED(IGO1)
              VALNEXT1 = DOCBUF((KEY1 - 1) * MAXX + KEYCOL + 1)

C             POINT TO NEXT VALUE IN SORTED LIST FROM FILE 2
              IGO2     = 1
              KEY2     = SORTED2(IGO2)
              VALNEXT2 = DOCBUF2( KEYCOL + 1, KEY2)
             
C             FIND OUTPUT VALUES
              MAXKEYS = MAX(IKEYS,IKEYS2)

              DO WHILE (IGO1 <= IKEYS .OR. IGO2 <= IKEYS2)

C                FIND KEY FOR THIS SORTED VALUE FROM FIRST FILE
                 IF (IGO1 > IKEYS .AND. IGO2 <= IKEYS2) THEN
C                   FILE 1 FINISHED BUT STILL IN LIST FROM FILE 2 
C                   POINT TO NEXT VALUE IN SORTED LIST FROM FILE 2
                    KEY2 = SORTED2(IGO2)

C                   SAVE VALUES FROM FILE 2
                    NVAL   = MAXX2 - 1
                    KEYNEW = KEYNEW + 1
                    CALL LUNDOCWRTDAT(NICDOCOUT,KEYNEW,DOCBUF2(2,KEY2),
     &                                NVAL,IRTFLG)

C                   INCREMENT IGO2
                    IGO2 = IGO2 + 1

                 ELSEIF (IGO1 <= IKEYS .AND. IGO2 > IKEYS2) THEN
C                   FILE 2 FINISHED BUT STILL IN LIST FROM FILE 1 
                    KEY1     = SORTED(IGO1)
                    DO IREG = 1,MAXX
                       DLIST(IREG) = DOCBUF((KEY1 - 1) * MAXX + IREG)
                    ENDDO

C                   PUSH DLIST INTO DOC. FILE
                    NVAL  = MAXX - 1
                    KEYNEW = KEYNEW + 1
                    CALL LUNDOCWRTDAT(NICDOCOUT,KEYNEW,DLIST(2),
     &                                NVAL,IRTFLG)

C                   INCREMENT IGO1
                    IGO1 = IGO1 + 1

                 ELSEIF (IGO1 <= IKEYS) THEN
C                   STILL IN LIST FROM FILE 1 AND FILE 2

                    IF (VALNEXT1 .LT. VALNEXT2)THEN
C                      NOT IN FILE 2, SAVE VALUES FROM FILE 1
                       DO IREG = 1,MAXX
                          DLIST(IREG) = DOCBUF((KEY1 - 1) * MAXX + IREG)
                       ENDDO

C                      PUSH DLIST INTO DOC. FILE
                       NVAL   = MAXX - 1
                       KEYNEW = KEYNEW + 1
                       CALL LUNDOCWRTDAT(NICDOCOUT,KEYNEW,DLIST(2),
     &                                NVAL,IRTFLG)

C                      INCREMENT IGO1
                       IGO1 = IGO1 + 1
                       IF (IGO1 <= IKEYS) THEN
                          KEY1     = SORTED(IGO1)
                          VALNEXT1 = DOCBUF((KEY1 - 1)*MAXX+KEYCOL + 1)
                       ENDIF

                    ELSEIF (VALNEXT1 == VALNEXT2) THEN
C                      SAME KEYCOL VALUES IN BOTH, SAVE VALUES FROM  2

                       NVAL   = MAXX2 - 1
                       KEYNEW = KEYNEW + 1
                       CALL LUNDOCWRTDAT(NICDOCOUT,KEYNEW,
     &                                DOCBUF2(2,KEY2),NVAL,IRTFLG)

C                      INCREMENT IGO1
                       IGO1 = IGO1 + 1
                       IF (IGO1 <= IKEYS) THEN
                          KEY1     = SORTED(IGO1)
                          VALNEXT1 = DOCBUF((KEY1 - 1)*MAXX+KEYCOL + 1)
                       ENDIF
C                      INCREMENT IGO2
                       IGO2 = IGO2 + 1
                       IF (IGO2 <= IKEYS2) THEN
                          KEY2     = SORTED2(IGO2)
                          VALNEXT2 = DOCBUF2(KEYCOL + 1, KEY2)
                       ENDIF

                    ELSE
C                      VALNEXT1 IS > VALNEXT2,    
C                      UPDATE VALNEXT2, SAVING ANY PASSED VALUES FROM 2
                       NEEDINC = 1
                       DO IGO2T = IGO2,IKEYS2
                          KEY2     = SORTED2(IGO2T)
                          VALNEXT2 = DOCBUF2(KEYCOL +1,KEY2)

                          IF (VALNEXT1 .GE. VALNEXT2) THEN
C                            SAVE VALUES FROM FILE 2, KEEP GOING if >

                             NVAL   = MAXX - 1
                             KEYNEW = KEYNEW + 1
                             CALL LUNDOCWRTDAT(NICDOCOUT,KEYNEW,
     &                                DOCBUF2(2,KEY2),NVAL,IRTFLG)
                          ELSE
C                            POINT TO THIS VALUES IN SORTED LIST FROM #2
                             NEEDINC = 0
                             EXIT
                          ENDIF
                       ENDDO   !DO IGO2T = IGO2,IKEYS2

C                      INCREMENT IGO2
                       IGO2     = IGO2T + NEEDINC
                       IF (IGO2 <= IKEYS2) THEN
                          KEY2     = SORTED2(IGO2)
                          VALNEXT2 = DOCBUF2(KEYCOL +1,KEY2)
                       ENDIF

                    ENDIF   !   VALNEXT1 IS > VALNEXT2
                 ENDIF !   IF (IGO1 <= IKEYS)
              ENDDO
           ENDIF

        ENDIF

C       CLOSE THE OUTPUT DOC. FILE(S)
9990    CLOSE(NDOCOUT)
        CLOSE(NDOCOUT2)

C       DEALLOCATE ALLOCATABLE ARRAYS
9995    IF (ALLOCATED(DLIST))   DEALLOCATE(DLIST)
        IF (ALLOCATED(SORTED2)) DEALLOCATE(SORTED2)
        IF (ALLOCATED(SORTED))  DEALLOCATE(SORTED)
        CLOSE(NDOCIN2)

C       DEALLOCATE DOC. FILE MEMORY
9997    IF (FCHAR(4:5) == 'ME') THEN
C          USED TWO INPUT DOC FILES
           IF (ASSOCIATED(DOCBUF2)) DEALLOCATE(DOCBUF2)
        ENDIF

9999    RETURN
	END


C       ----------------------- SORTIT --------------------------------

        SUBROUTINE SORTIT(DOCBUF,MAXX,MAXY,KEYCOL,SORTED,
     &                    IKEYS,RET_KEY,UNIQUE,IRTFLG)

        REAL, DIMENSION(MAXX,MAXY), INTENT(IN)  :: DOCBUF
        REAL, DIMENSION(MAXY), INTENT(INOUT)    :: SORTED
        INTEGER,INTENT(IN)                      :: MAXX,MAXY,KEYCOL
        LOGICAL,INTENT(IN)                      :: RET_KEY,UNIQUE
        INTEGER,INTENT(OUT)                     :: IKEYS,IRTFLG

        REAL, DIMENSION(MAXY)                   :: RDUM,RKEYARAY

C       TRANSFER DATA TO SORT INPUT ARRAYS
        IKEYS    = 0
        KEYCOLP1 = KEYCOL + 1

        DO IROW = 1, MAXY
           IF (DOCBUF(1,IROW) > 0) THEN
C             KEY IS USED
              IKEYS = IKEYS + 1

              IF (RET_KEY) THEN
C                RETURN THE KEY NUMBER IN SORTED
                 SORTED(IKEYS) = IROW
              ELSE
C                RETURN THE VALUE IN SORTED
                 SORTED(IKEYS) = DOCBUF(KEYCOLP1,IROW)
              ENDIF

C             SORT BY THE VALUE IN COLUMN: KEYCOL
              RKEYARAY(IKEYS)  = DOCBUF(KEYCOLP1,IROW)
              RDUM(IKEYS)      = 0.0
           ENDIF
        ENDDO

C       SORT BY VALUE IN RKEYARAY, ONLY INTERESTED IN SORTED
        CALL SORT(RKEYARAY,RDUM,SORTED,IKEYS)

        IF (UNIQUE) THEN
C          ONLY WANT UNIQUE VALUES FROM COLUMN: KEYCOL
           IKEYSNEW = 0
C          INITIALIZE DLAST TO ENUSRE KEEPING FIRST VALUE IN KEYCOL
           IT       = SORTED(1)
           DLAST    = DOCBUF(KEYCOLP1,IT) + 1000

           DO I = 1,IKEYS
              IT  = SORTED(I)
              VAL = DOCBUF(KEYCOLP1,IT)
              IF (VAL .NE. DLAST) THEN
C                VAL IS NOT SAME AS PREVIOUS VALUE IN KEYCOL, KEEP IT
                 IKEYSNEW         = IKEYSNEW + 1
                 SORTED(IKEYSNEW) = IT
                 DLAST            = VAL
              ENDIF
           ENDDO
           IKEYS = IKEYSNEW
        ENDIF
        IRTFLG = 0

        RETURN
        END
 
C       ----------------------- DOCCREATE ----------------------------

	SUBROUTINE DOCCREATE(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        INTEGER,INTENT(IN)               :: MAXDIM
	CHARACTER(LEN=MAXNAM)            :: DOCNAM3
	CHARACTER(LEN=1)                 :: NULL
        LOGICAL                          :: ERRI2
        REAL,ALLOCATABLE                 :: DLIST(:)

        INTEGER, DIMENSION(1)            :: ILIST
        COMMON   ILIST

        INTEGER, PARAMETER                ::  NDOCOUT = 72

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

	NULL = CHAR(0)

C       CREATE OUTPUT DOC FILE ------------------------- 'DOC CREATE'

        CALL FILERD(DOCNAM3,NLET,NULL,'OUTPUT DOCUMENT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL OPENDOC(DOCNAM3,.TRUE.,NLET,NDOCOUT,NICDOCOUT,.FALSE.,' ',
     &                     .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

10      CALL RDPRI1S(KEYCOL,NOT_USED,
     &               'REGISTER TO BE FILLED (0 IS KEY)',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (ERRI2(KEYCOL,IDUM,1,0,6,0,0)) GOTO 10

        NUMB = MAXDIM
        CALL RDPRAI(ILIST,MAXDIM,NUMB,1,MAXDIM,'NUMBERS',
     &                 NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        NLIST = MAX(KEYCOL,1)
        ALLOCATE(DLIST(NLIST),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'DLIST IN DOCCREATE',NLIST)
           RETURN
        ENDIF

C       FILL UNUSED COLUMNS WITH ZERO'S
        DO I = 1, NLIST
           DLIST(I) = 0.0
        ENDDO

C       IF FILLING KEYS PUT A 1.0 IN FIRST REGISTER COL.
        IF (KEYCOL == 0) DLIST(1) = 1.0

        DO I=1,NUMB
           IF (KEYCOL <= 0) THEN
              IKEY          = ILIST(I)
           ELSE
              IKEY          = I
              DLIST(KEYCOL) = ILIST(I)
           ENDIF

C          PUSH DLIST INTO DOC. FILE
           CALL LUNDOCWRTDAT(NICDOCOUT,IKEY,DLIST,NLIST,IRTFLG)
        ENDDO

        IF (ALLOCATED(DLIST)) DEALLOCATE(DLIST)

#ifdef USE_MPI
        IF (MYPID == 0) THEN
           CALL FLUSHFILE(NICDOCOUT)
           CLOSE(NICDOCOUT)
        ENDIF
        CALL MPI_BARRIER(ICOMM,IERR)
#else
        IF (MYPID <= 0) CLOSE(NDOCOUT)
#endif

        RETURN
        END


C       ----------------------- DOCRAN --------------------------------

	SUBROUTINE DOCRAN(MAXX, MAXY, DOCBUF, DLIST)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

	CHARACTER(LEN=MAXNAM)             :: DOCNAM3
	CHARACTER(LEN=1)                  :: NULL

        REAL                              :: DOCBUF(MAXX*MAXY)
        REAL                              :: DLIST(MAXX)
        INTEGER, ALLOCATABLE              :: KEYLIST(:)
        LOGICAL                           :: NEWFILE

        INTEGER, PARAMETER                ::  NDOCOUT = 71

	NULL = CHAR(0)

C       RANDOMLY SAMPLE THE INPUT DOC FILE---------------- 'DOC RAN'

        CALL FILERD(DOCNAM3,NLET,NULL,'OUTPUT DOCUMENT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL OPENDOC(DOCNAM3,.TRUE.,NLET,NDOCOUT,NICDOCOUT,.FALSE.,' ',
     &                     .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL RDPRM1S(PERCENT,NOT_USED,'PERCENT WANTED',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (MAXX <= 0 .OR. MAXY <= 0) GOTO 9999

        ALLOCATE(KEYLIST(MAXY),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'KEYLIST IN DOCS1',MAXY)
           RETURN
        ENDIF

        IGOT = 0
        DO IKEY = 1,MAXY
C           KEYS MAY NOT BE CONSECUTIVE SO MUST MAKE A LIST
            ILOC = (IKEY - 1) * MAXX + 1
            IF (DOCBUF(ILOC) > 0) THEN
               IGOT          = IGOT + 1
               KEYLIST(IGOT) = IKEY
            ENDIF
        ENDDO

C       FIND NUMBER OF NEEDED KEYS
        NEEDED = PERCENT * IGOT * 0.01
        IRAN   = 0
        NLIST  = MAXX - 1
        DO
          IF (IRAN .GE. NEEDED) EXIT

          IRAN = IRAN + 1
C         CREATE RANDOM IVAL IN RANGE 0...IGOT-1
          CALL RANDOM_NUMBER(OUT)

          IVAL    = 1.5 + OUT * FLOAT(IGOT-1)
          IKEY    = KEYLIST(IVAL)

          DO ICOL = 2,MAXX
             DLIST(ICOL-1) = DOCBUF((IKEY - 1) * MAXX + ICOL)
          ENDDO

C         PUSH DLIST INTO DOC. FILE
          CALL LUNDOCWRTDAT(NICDOCOUT,IKEY,DLIST,NLIST,IRTFLG)

C         SELECT RANDOMLY WITHOUT DUPLICATION OF SELECTED VALUES
          KEYLIST(IVAL) = KEYLIST(IGOT)
          IGOT          = IGOT - 1
        ENDDO

C       DEALLOCATE ALLOCATABLE ARRAYS
        IF (ALLOCATED(KEYLIST)) DEALLOCATE(KEYLIST)
 
9999    CLOSE(NDOCOUT)

	END

C       ----------------------- DOCRENUMBER --------------------------

	SUBROUTINE DOCRENUMBER(MAXX, MAXY, DOCBUF, DLIST)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        INTEGER, INTENT(IN)               :: MAXX,MAXY
        REAL, DIMENSION(MAXX*MAXY)        :: DOCBUF
        REAL, DIMENSION(MAXX)             :: DLIST
	CHARACTER(LEN=MAXNAM)             :: DOCNAM3
	CHARACTER(LEN=1)                  :: NULL
        LOGICAL                           :: NEWFILE

        INTEGER, PARAMETER                ::  NDOCOUT = 72

	NULL = CHAR(0)

C       ------------- RENUMBER THE INPUT DOC FILE-------- 'DOC RENUMBER'

        CALL FILERD(DOCNAM3,NLET,NULL,'OUTPUT DOCUMENT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL OPENDOC(DOCNAM3,.TRUE.,NLET,NDOCOUT,NICDOCOUT,.FALSE.,' ',
     &                     .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       DOC REN OPERATION USES NON-STANDARD DOCBUF (WITHOUT KEYS)
C       HAVING LISTING OF ALL LINES IN THE DOC FILE IN ORDER
C       write(6,*) 'maxx,maxy:',maxx,maxy

        IF (MAXX <= 0 .OR. MAXY <= 0) GOTO 9999
           
        NLIST = MAXX

        DO IKEY = 1,MAXY
            DO ICOL = 1,MAXX
               DLIST(ICOL) = DOCBUF((IKEY - 1) * MAXX + ICOL)
            ENDDO

C           PUSH DLIST INTO DOC. FILE
            CALL LUNDOCWRTDAT(NICDOCOUT,IKEY,DLIST,NLIST,IRTFLG)
        ENDDO

9999    CLOSE(NDOCOUT)

	END

C       -------------RE KEY THE INPUT DOC FILE --------------- 'DOC KEY'

	SUBROUTINE DOCREKEY(MAXX, MAXY, DOCBUF, DLIST)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        INTEGER, INTENT(IN)   :: MAXX,MAXY
        REAL                  :: DOCBUF(MAXX,MAXY)
        REAL                  :: DLIST(MAXX)
	CHARACTER(LEN=MAXNAM) :: DOCNAM3

        LOGICAL               :: ADDEXT,GETNAME
        LOGICAL               :: ISOLD,APPEND,MESSAGE,NEWFILE
        INTEGER               :: IRTFLG, NLET, NICDOCOUT, NLIST 
        INTEGER               :: NEWKEY, IKEY, ICOUNT, ICOL, ICOUNTLAS        

        INTEGER, PARAMETER    :: NDOCOUT = 72

        ADDEXT  = .TRUE.
        GETNAME = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .FALSE.
        MESSAGE = .TRUE.
        IRTFLG  = -8         ! NO IC USE

        CALL OPENDOC(DOCNAM3,ADDEXT,NLET,NDOCOUT,NICDOCOUT,GETNAME,
     &           'RE-KEYED OUTPUT DOCUMENT',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
C       write(6,*) 'maxx,maxy:',maxx,maxy
        IF (MAXX <= 0 .OR. MAXY <= 0) GOTO 9999

C       COPY ONE OR MORE OTHER REGISTERS TO OUTPUT
        NLIST  = MAXX
        IF (FCHAR(8:8) == 'O') NLIST = 1

        NEWKEY = 0
        DO IKEY = 1,MAXY

C          DOCBUF HAS ICOUNT IN FIRST COL
           ICOUNT = DOCBUF(1,IKEY)

           IF (ICOUNT > 0) THEN
C             GOT VALID DOC FILE DATA LINE

C             PUT KEY IN FIRST COL OF OUTPUT DOC FILE
              DLIST(1) = IKEY

C             COPY NLIST INPUT DOC FILE COLUMNS
              DO ICOL = 2,NLIST
                  DLIST(ICOL) = DOCBUF(ICOL,IKEY)
              ENDDO

C             PUSH DLIST INTO OUTPUT DOC. FILE
              NEWKEY = NEWKEY + 1
              CALL LUNDOCWRTDAT(NICDOCOUT,NEWKEY,DLIST,NLIST,IRTFLG)

              ICOUNTLAS = ICOUNT
           ENDIF
        ENDDO
 
9999    CLOSE(NDOCOUT)

	END



C       ----------------------- DOCCOMBINE --------------------------

	SUBROUTINE DOCCOMBINE()

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        INCLUDE 'F90ALLOC.INC'
        REAL, DIMENSION(:,:), POINTER      :: DOCBUF
	CHARACTER(LEN=MAXNAM)              :: FILPAT,DOCNAM1,DOCNAM3
        LOGICAL                            :: NEWFILE

        INTEGER, PARAMETER                ::  NDOCIN  = 70
        INTEGER, PARAMETER                ::  NDOCIN2 = 71
        INTEGER, PARAMETER                ::  NDOCOUT = 72

C       COMBINE THE INPUT DOC FILES -------------------- 'DOC COMBINE'

C       SPACE FOR DOC FILE LIST FROM CMLIMIT
        NILMAX = NIMAX

C       ASK FOR DOC FILE LIST
        CALL FILELIST(.TRUE.,NDOCIN2,FILPAT,NLETP,INUMBR,NILMAX,NFILE,
     &      'TEMPLATE FOR DOC. FILE SERIES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NFILE > 0)  THEN
           WRITE(NOUT,2001) NFILE
2001       FORMAT('  Number of document files to be combined: ',I0)
        ELSE
           CALL ERRT(101,'No document files entered!',IER)
           GOTO 9999
        ENDIF

        CALL OPENDOC(DOCNAM3,.TRUE.,NLET,NDOCOUT,NICDOCOUT,.TRUE.,
     &       'OUTPUT DOCUMENT',.FALSE.,.TRUE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IKEYNOW = 0

        DO IFILE = 1,NFILE
C          DOC COM OPERATION USES STANDARD DOCBUF (WITH KEYS)

C          MAKE DOC FILE NAME
           CALL FILGET(FILPAT,DOCNAM1,NLETP,INUMBR(IFILE),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           CALL FILNAMANDEXT(DOCNAM1,DATEXC,DOCNAM1,NLET,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
 
           MAXX = 0
           MAXY = 0
           CALL GETDOCDAT(' ',.FALSE.,DOCNAM1,
     &                 NDOCIN,.TRUE.,MAXX, MAXY,DOCBUF,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           NLIST = MAXX - 1

           DO IKEY = 1,MAXY
              ICOUNT = DOCBUF(1,IKEY)
              IF (ICOUNT > 0) THEN
C                KEY EXISTS, PUSH LINE INTO COMBINED DOC. FILE
                 IKEYNOW = IKEYNOW + 1
                 CALL LUNDOCWRTDAT(NICDOCOUT,IKEYNOW,DOCBUF(2,IKEY),
     &                             NLIST,IRTFLG)
              ENDIF
           ENDDO
           CLOSE(NDOCIN)

C          DEALLOCATE DOC. FILE MEMORY
           IF (ASSOCIATED(DOCBUF)) DEALLOCATE(DOCBUF)
        ENDDO

9999    CLOSE(NDOCOUT)
        CLOSE(NDOCIN)

C       DEALLOCATE DOC. FILE MEMORY
        IF (ASSOCIATED(DOCBUF)) DEALLOCATE(DOCBUF)

        RETURN
	END



C       ----------------------- DOCDOWN ----------------------------

	SUBROUTINE DOCDOWN()

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	CHARACTER(LEN=MAXNAM)             :: DOCNAM
	CHARACTER(LEN=160)                :: RECLIN
        REAL, DIMENSION(9)                :: DLIST
        LOGICAL                           :: WARNIT,NEWFORM


	DATA NDOCINT,NDOCOUTT/70,72/

        WARNIT = .TRUE.

C       OPEN INPUT DOC. FILE
        CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOCINT,NDOCIN,.TRUE.,
     &           'INPUT DOCUMENT',.TRUE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       OPEN OUTPUT DOC. FILE
        CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOCOUTT,NDOCOUT,.TRUE.,
     &          'OUTPUT DOCUMENT',.FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        DO
          READ(NDOCIN,84,IOSTAT=IERR) RECLIN
84        FORMAT(A160)

          IF (IERR .LT. 0) THEN
C             END OF FILE
              GOTO 9999
           ENDIF

          IF (RECLIN(2:2) == ';') THEN
              WRITE(NDOCOUT,90,IOSTAT=IERR) RECLIN
90            FORMAT(A)
             
              CYCLE
          ENDIF

          NEWFORM = .TRUE.
          READ(RECLIN,*,IOSTAT=IERR) IKEY,ICOUNT

          IF (IERR > 0) THEN
C            ERROR ON READ, TRY OLD DOC. FILE FORMAT
             READ(RECLIN,83,IOSTAT=IERR) IKEY,ICOUNT
83           FORMAT(I6,I1,10000F12.6)

             NEWFORM = .FALSE.

             IF (IERR > 0) THEN
C               ERROR ON READ USING OLD FORMAT ALSO, RETURN
                WRITE(NOUT,91) RECLIN
91              FORMAT(' *** UNABLE TO INTERPRET DOC FILE LINE: ',A)
                CALLERRT(100,'DOCCOPY',NE)
                GOTO 9999
             ENDIF
          ENDIF

          IF (ICOUNT <= 0) THEN
              WRITE(NOUT,*) ' EMPTY DOCUMENT FILE LINE SKIPPED'
              CYCLE

          ELSEIF (IKEY .LT. 0) THEN
              WRITE(NOUT,*) ' CONTINUATION LINE SKIPPED IN DOC FILE'
              CYCLE

          ELSEIF (IKEY == 0) THEN
C            KEY THAT WILL NOT FIT IN DBUF SENDS ERROR MSG.
             WRITE(NOUT,*)' SKIPPED ILLEGAL KEY NUMBER: 0 IN DOC FILE'
             CYCLE

          ELSEIF (IKEY > 999999) THEN
C            KEY THAT WILL NOT FIT IN OLD DOC FILE SENDS ERROR MSG.
             IF (WARNIT) THEN
                WRITE(NOUT,93) IKEY
93              FORMAT('  ** KEY: ',I9,'  NOT RETRIEVED')
                WARNIT = .FALSE.
             ENDIF

          ELSEIF (ICOUNT > 9) THEN
C            KEY THAT WILL NOT FIT IN OLD DOC FILE SENDS ERROR MSG.
             IF (WARNIT) THEN
                WRITE(NOUT,*) ' ** REGISTERS > 9 NOT RETRIEVED'
                WARNIT = .FALSE.
             ENDIF
             ICOUNT = 9
          ENDIF

          BACKSPACE(NDOCINT)
          IF (NEWFORM) THEN
C            TRY NEW DOC. FILE FORMAT
             READ(NDOCIN,*,IOSTAT=IERR)IKEYT,ICOUNTT,
     &                                 (DLIST(I),I=1,ICOUNT)
C            IF ERROR ON READ, TRY OLD DOC. FILE FORMAT
             IF (IERR .NE. 0) THEN
C               TRY READING AGAIN USING OLD FORMAT
                BACKSPACE(NDOCINT)
                NEWFORM = .FALSE.
             ENDIF
          ENDIF
          IF (.NOT. NEWFORM) THEN
C            TRY OLD DOC. FILE FORMAT
             READ(NDOCIN,83,IOSTAT=IERR) IKEYT,ICOUNTT,
     &                                   (DLIST(I),I=1,ICOUNT)
          ENDIF

          IF (IERR == 0) THEN
             IF (IKEY <= 99999) THEN
                WRITE(NDOCOUT,95)IKEY,ICOUNT,(DLIST(I),I=1,ICOUNT)
95              FORMAT(I5,' ',I1,9G12.3)
             ELSE
                WRITE(NDOCOUT,96)IKEY,ICOUNT,(DLIST(I),I=1,ICOUNT)
96              FORMAT(I6,I1,9G12.3)
             ENDIF
          ENDIF
       ENDDO

9999   CLOSE(NDOCIN)
       CLOSE(NDOCOUT)

       END


C       ----------------------- DOCBOOT--------------------------------

	SUBROUTINE DOCBOOT(MAXX, MAXY, DOCBUF, DLIST)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

	CHARACTER(LEN=MAXNAM)              :: DOCNAM3
	CHARACTER(LEN=1)                   :: NULL

        REAL                               :: DOCBUF(MAXX*MAXY)
        REAL                               :: DLIST(MAXX)
        INTEGER, ALLOCATABLE               :: KEYLIST(:)
        INTEGER, ALLOCATABLE               :: KEYLISTOUT(:)
        LOGICAL                            :: NEWFILE

        INTEGER, PARAMETER                 ::  NDOCOUT = 71

	NULL = CHAR(0)

C       RANDOMLY SAMPLE THE INPUT DOC FILE---------------- 'DOC BOOT'

        CALL FILERD(DOCNAM3,NLET,NULL,'OUTPUT DOCUMENT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL OPENDOC(DOCNAM3,.TRUE.,NLET,NDOCOUT,NICDOCOUT,.FALSE.,' ',
     &                     .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN


c        CALL RDPRI1S(IWANT,NOT_USED,
c     &               'NUMBER OF SELECTIONS WANTED',IRTFLG)
c        IF (IRTFLG .NE. 0) RETURN
 
        IF (MAXX <= 0 .OR. MAXY <= 0) GOTO 9999

        ALLOCATE(KEYLIST(MAXY),KEYLISTOUT(MAXY),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'DOCS1; KEYLIST',2*MAXY)
           RETURN
        ENDIF

        IGOT = 0
        DO IKEY = 1,MAXY
C           KEYS MAY NOT BE CONSECUTIVE SO MUST MAKE A LIST
            ILOC = (IKEY - 1) * MAXX + 1
            IF (DOCBUF(ILOC) > 0) THEN
C              KEY IS IN USE
               IGOT          = IGOT + 1
               KEYLIST(IGOT) = IKEY
            ENDIF
        ENDDO
        IWANT = IGOT

C       SELECT IWANT ENTRIES RANDOMLY WITH POSSIBLE DUPLICATION 
C       OF THE SELECTED ENTRIES
        DO IRAN = 1,IWANT

C         CREATE RANDOM IVAL IN RANGE 1...IGOT
          CALL RANDOM_NUMBER(OUT)

          IVAL             = MIN(IGOT,MAX(1,INT(OUT*IGOT+0.5)))
          KEYLISTOUT(IRAN) = KEYLIST(IVAL)
        ENDDO

C       MAKE LIST OF THE VALUES IN COL 1 FOR THE SELECTED ENTRIES
        DO I = 1,IWANT
           IKEY       = KEYLISTOUT(I)
           KEYLIST(I) = DOCBUF((IKEY - 1) * MAXX + 2)
        ENDDO

C       SORT KEYS ORDERING KEYLIST BY VALUES IN FIRST COLUMN
        CALL SORTINT(KEYLIST, KEYLISTOUT, IWANT)

        NLIST  = MAXX - 1
        DO I = 1,IWANT

          IKEY = KEYLISTOUT(I)
          DO ICOL = 2,MAXX
             DLIST(ICOL-1) = DOCBUF((IKEY - 1) * MAXX + ICOL)
          ENDDO

C         PUSH DLIST INTO DOC. FILE WITH RENUMBERING OF KEYS
          CALL LUNDOCWRTDAT(NICDOCOUT,I,DLIST,NLIST,IRTFLG)
        ENDDO

C       DEALLOCATE ALLOCATABLE ARRAYS
        IF (ALLOCATED(KEYLIST))    DEALLOCATE(KEYLIST)
        IF (ALLOCATED(KEYLISTOUT)) DEALLOCATE(KEYLISTOUT)
 
9999    CLOSE(NDOCOUT)

	END

