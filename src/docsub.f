
C ++********************************************************************
C
C DOCSUB  REMOVED FROM DOCS1                     JAN 2011 ARDEAN LEITH
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
C  DOCSUB
C                                                                 
C  PURPOSE: 'DOC SUB'
C       SUBTRACT THE 2ND INPUT DOC FILE----------------- 'DOC SUBTRACT'
C       AND CONTENTS OF TWO DOC FILES ----------------------- 'DOC AND'
C   
C  NOTE: TSWITCH SAYS 'DOC' IS A 2 LETTER OP
C                                                                    
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE DOCSUB(MAXX, MAXY, DOCBUF,DLIST,NLIST)

C       NLIST IS 1 + NUM OF REGISTERS SINCE DOCBUF CONTAINS KEY ALSO

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        REAL                         :: DOCBUF(MAXX*MAXY)
        REAL                         :: DLIST(NLIST)

	CHARACTER(LEN=MAXNAM)        :: DOCNAM2,DOCNAM3
	CHARACTER(LEN=1)             :: NULL = CHAR(0)
	CHARACTER(LEN=1)             :: ANSW
        LOGICAL                      :: NEWFILE,ERRI2,RENUMBER

        REAL,    ALLOCATABLE         :: SORTED(:),SORTED2(:)
        INTEGER, ALLOCATABLE         :: KEYLIST(:)

        INCLUDE 'F90ALLOC.INC'
        REAL, POINTER                :: DOCBUF2(:,:)
	CHARACTER(LEN=80)            :: PROMPT
        LOGICAL                      :: SENDIT,REP_KEY,UNIQUE
        LOGICAL                      :: REVERSE
        LOGICAL                      :: APPENDOK,MESSAGE

        INTEGER, PARAMETER           :: NDOCIN2  = 71
        INTEGER, PARAMETER           :: NDOCOUT  = 72
        INTEGER, PARAMETER           :: NDOCOUT2 = 73


        IF (FCHAR(4:5) .EQ. 'SU') THEN
C                      123456789 123456789 123456789 1234567890
           PROMPT   = 'COLUMN TO BE SUBTRACTED BY (0 IS KEY)' 
           NLETP    = 37
           RENUMBER = .FALSE.

C          SUBTRACT USES 2 INPUT DOC. FILES
C          MAXX2 IS 1 + NUM OF REGISTERS SINCE DOCBUF CONTAINS KEY ALSO
           MAXX2  = 0
           MAXY2  = 0
           CALL GETDOCDAT('SECOND INPUT DOCUMENT',.TRUE.,DOCNAM2,
     &                  NDOCIN2,.TRUE.,MAXX2, MAXY2,DOCBUF2,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

        ELSEIF (FCHAR(4:5) .EQ. 'AN') THEN
C                      123456789 123456789 123456789 1234567890
           PROMPT   = 'COLUMN TO BE CHECKED (0 IS KEY)' 
           NLETP    = 31
           RENUMBER = .FALSE.

C          AND USES 2 INPUT DOC. FILES
C          MAXX2 IS 1 + NUM OF REGISTERS SINCE DOCBUF CONTAINS KEY ALSO
           MAXX2  = 0
           MAXY2  = 0
           CALL GETDOCDAT('SECOND INPUT DOCUMENT',.TRUE.,DOCNAM2,
     &                  NDOCIN2,.TRUE.,MAXX2, MAXY2,DOCBUF2,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
           CALL ERRT(101,"UNKNOWN 'DOC' OPERATION",IDUM)
           RETURN 
        ENDIF

C       OPEN OUTPUT DOCUMENT FILE
        CALL FILERD(DOCNAM3,NLET,NULL,'OUTPUT DOCUMENT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9990

        APPENDOK = .FALSE.
        MESSAGE  = .TRUE.
        CALL OPENDOC(DOCNAM3,.TRUE.,NLET,NDOCOUT,NICDOCOUT,.FALSE.,
     &               ' ',.FALSE.,APPENDOK,MESSAGE,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9990


        IF (FCHAR(4:5) .EQ. 'SU' .OR.
     &      FCHAR(4:5) .EQ. 'AN' ) THEN
C          SUBTRACT THE 2ND INPUT DOC FILE--------------- 'DOC SUBTRACT'
C          AND THE INPUT DOC FILES --------------------------- 'DOC AND'

           CALL RDPRI1S(KEYCOL,NOT_USED,PROMPT(1:NLETP),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9990

           REVERSE = (KEYCOL < 0)
           KEYCOL  = ABS(KEYCOL)

           IF (MAXX > 0 .AND. MAXY > 0) THEN
              IF (ERRI2(KEYCOL,IDUM,1, 0,MAXX-1,0,0)) GOTO 9990
           ENDIF
        ENDIF

        !write(6,*) ' maxx:',maxx,maxy,' maxx2:',maxx2,maxy2
     
        IF ((FCHAR(4:5) .EQ. 'AN') .AND. 
     &      (MAXX  < 1 .OR. MAXY  < 1) .AND. 
     &      (MAXX2 < 1 .OR. MAXY2 < 1)) THEN
C          BOTH FILES EMPTY, THUS EMPTY OUTPUT ----------- 'AN' BOTH MT
           GOTO 9990

        ELSEIF ((FCHAR(4:5) .EQ. 'SU') .AND.  
     &          (MAXX  < 1 .OR. MAXY  < 1)) THEN 
C          SUBTRACT BUT FIRST FILE IS EMPTY, ------------ 'SU' 1'S MT
C          EMPTY OUTPUT
           GOTO 9990

        ELSEIF ((FCHAR(4:5) .EQ. 'SU') .AND.
     &          (MAXX2  < 1 .OR. MAXY2  < 1)) THEN 
C          SUBTRACT,  SECOND FILE IS EMPTY.  COPY VALUES FROM FIRST FILE

           ALLOCATE(SORTED(MAXY),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'DOCS1; SORTED',MAXY)
              GOTO 9990
           ENDIF

C          SORTED RETURNS LIST OF KEYS SORTED BY USING VALUE IN KEYCOL
           UNIQUE = .FALSE.
           CALL SORTIT(DOCBUF,MAXX,MAXY,KEYCOL,SORTED,
     &                 IKEYS,.TRUE.,UNIQUE,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9990

           KEYNEW = 0
           NVAL   = MAXX - 1
           DO IKEYT = 1,MAXY
              ICOUNT = DOCBUF((IKEYT - 1) * MAXX + 1)
              IF (ICOUNT > 0 ) THEN
C                KEY EXISTS PUSH DLIST FROM FILE 1 INTO OUTPUT FILE
                 DO IREG = 1,MAXX
                    DLIST(IREG) = DOCBUF((IKEYT - 1) * MAXX + IREG)
                 ENDDO

C                PUSH DLIST INTO DOC. FILE
                 KEYNEW = KEYNEW + 1
                 CALL LUNDOCWRTDAT(NICDOCOUT,KEYNEW,DLIST(2),
     &                             NVAL,IRTFLG)
              ENDIF
           ENDDO
           GOTO 9990

        ELSEIF (FCHAR(4:5) .EQ. 'SU' .OR.
     &          FCHAR(4:5) .EQ. 'AN' .AND. KEYCOL .GT. 0) THEN
C          NEED ONE OR MORE SORTED LISTS

C          SUBTRACT THE 2ND INPUT DOC FILE--------------- 'DOC SUBTRACT'
C          AND THE INPUT DOC FILES --------------------------- 'DOC AND'

           ALLOCATE(SORTED(MAXY),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'DOCS1; SORTED',MAXY)
              GOTO 9990
           ENDIF

           UNIQUE = .FALSE.
C          SORTED RETURNS LIST OF KEYS SORTED BY USING VALUE IN KEYCOL
           CALL SORTIT(DOCBUF,MAXX,MAXY,KEYCOL,SORTED,
     &                 IKEYS,.TRUE.,UNIQUE,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9990

           IF (FCHAR(4:5) .EQ. 'SU' .OR. FCHAR(4:5) .EQ. 'ME') THEN
C             SORT THE SECOND LIST OF KEYS BUT RETURN VALUES IN SORTED2
              ALLOCATE(SORTED2(MAXY2),STAT=IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                 CALL ERRT(46,'DOCS1; SORTED',MAXY2)
                 GOTO 9990
              ENDIF

C             RETURN SORTED LIST OF VALUES IN KEYCOL
              REP_KEY = .FALSE.
              UNIQUE  = .FALSE. 
              CALL SORTIT(DOCBUF2(1,1),MAXX2,MAXY2,KEYCOL,SORTED2,
     &                    IKEYS2,REP_KEY,UNIQUE,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9990

          ELSEIF (FCHAR(4:5) .EQ. 'AN') THEN
C             SORT THE SECOND LIST OF KEYS BUT RETURN VALUES IN SORTED2
              ALLOCATE(SORTED2(MAXY2),STAT=IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                 CALL ERRT(46,'DOCS1; SORTED',MAXY2)
                 GOTO 9990
              ENDIF

C             RETURN SORTED LIST OF VALUES IN KEYCOL
              REP_KEY = .FALSE.
              UNIQUE  = .FALSE. 
              CALL SORTIT(DOCBUF2(1,1),MAXX2,MAXY2,KEYCOL,SORTED2,
     &                    IKEYS2,REP_KEY,UNIQUE,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9990
          ENDIF
        ENDIF

        IF (FCHAR(4:5) .EQ. 'SU') THEN
C          VALUES IN FIRST FILE KEYCOL BUT NOT IN 2ND: ------- 'DOC SUB.'

           IF (KEYCOL .EQ. 0) THEN
C             EXCLUSIVE SUBTRACTING BY KEY
              DO IKEYT = 1,MAXY
                 ICOUNT1 = DOCBUF((IKEYT - 1) * MAXX + 1)
                 IF (ICOUNT1 .GT. 0 .AND. IKEYT .LE. MAXY2) THEN
C                   KEY EXISTS IN FIRST FILE, CHECK EXISTANCE IN 2ND
                     
                    ICOUNT2 = DOCBUF2(1,IKEYT)
                    IF (ICOUNT2 .EQ. 0) THEN
C                      NEED TO KEEP THIS KEY
                       DO IREG = 1,MAXX
                          DLIST(IREG) = DOCBUF((IKEYT - 1) *MAXX +IREG)
                       ENDDO

C                      PUSH DLIST INTO DOC. FILE
                       NVAL  = MAXX - 1
                       CALL LUNDOCWRTDAT(NICDOCOUT,IKEYT,DLIST(2),
     &                                NVAL,IRTFLG)
                    ENDIF

                 ELSEIF (ICOUNT1 .GT. 0) THEN
C                   KEY EXISTS IN FIRST FILE, 2ND FILE FINISHED
                     
                    DO IREG = 1,MAXX
                       DLIST(IREG) = DOCBUF((IKEYT - 1) * MAXX + IREG)
                    ENDDO

C                   PUSH DLIST INTO DOC. FILE
                    NVAL  = MAXX - 1
                    CALL LUNDOCWRTDAT(NICDOCOUT,IKEYT,DLIST(2),
     &                                NVAL,IRTFLG)
                 ENDIF
              ENDDO

           ELSE                    
C             SUBTRACTING BY COLUMN OTHER THAN KEY
C             POINT TO NEXT VALUE IN SORTED LIST FROM FILE 2
              IGO2     = 1
              VALNEXT2 = SORTED2(IGO2)
             
C             FIND OUTPUT VALUES
              DO IK = 1,IKEYS
C                FIND VALUE IN KEYCOL OF FIRST DOC. FILE

C                FIND KEY FOR THIS SORTED VALUE FROM FIRST FILE
                 KEY1 = SORTED(IK)
C                FIND SORTED VALUE FROM FIRST FILE
                 VAL1 = DOCBUF((KEY1 - 1) * MAXX + KEYCOL + 1)

                 SENDIT = .TRUE.
                 IF (IGO2 .GT. IKEYS2) THEN
C                   VALUE DOES NOT EXIST IN FILE 2 (ALL FILE 2'S DONE)
                    SENDIT = .TRUE.

                 ELSE
C                   MUST CHECK TO SEE IF VALUE EXISTS IN FILE 2
                    IF (VAL1 .LT. VALNEXT2) THEN
C                      VALUE 1 IS BELOW NEXT VALUE 2, SAVE THIS LINE
                       SENDIT = .TRUE.

                    ELSEIF (VAL1 .EQ. VALNEXT2 .AND.
     &                      IGO2 .LE. IKEYS2) THEN
C                      VALUE 1 = NEXT VALUE 2, MUST INCREASE IGO2
                       SENDIT   = .FALSE.
                       IGO2     = IGO2 + 1
                       IF (IGO2 .LE. IKEYS2) VALNEXT2 = SORTED2(IGO2)

                    ELSEIF (VAL1 .GT. VALNEXT2 .AND.
     &                      IGO2 .LT. IKEYS2) THEN
C                      VALUE 1 > NEXT VALUE 2, MUST INCREASE VALNEXT2
                       SENDIT = .TRUE.
                       DO IGO2T = IGO2+1,IKEYS2
                          IF (VAL1 .EQ. SORTED2(IGO2T)) THEN
C                            KEEP GOING TILL VALUE 2 > 
                             SENDIT = .FALSE.
                          ELSEIF (VAL1 .LT. SORTED2(IGO2T)) THEN
C                            POINT TO THIS VALUES IN SORTED LIST FROM #2
                             IGO2     = IGO2T
                             VALNEXT2 = SORTED2(IGO2)
                             EXIT
                         ENDIF
                       ENDDO
                    ENDIF                    
                 ENDIF

                 IF (SENDIT) THEN
C                   VALUE FROM KEYCOL NOT IN 2'ND DOC. FILE
C                   PUT LINE FROM FIRST DOC. FILE INTO OUTPUT FILE
                    DO IREG = 1,MAXX
                       DLIST(IREG) = DOCBUF((KEY1 - 1) * MAXX + IREG)
                    ENDDO

C                   PUSH DLIST INTO DOC. FILE
                    NVAL  = MAXX - 1
                    IKEYT = KEY1
                    CALL LUNDOCWRTDAT(NICDOCOUT,IKEYT,DLIST(2),
     &                                NVAL,IRTFLG)
                 ENDIF
              ENDDO
           ENDIF

        ELSEIF (FCHAR(4:5) .EQ. 'AN') THEN
C          VALUES IN FIRST FILE KEYCOL AND ALSO IN 2ND: ------ 'DOC AND'

           IF (KEYCOL .EQ. 0) THEN
C             'ANDING' BY KEY
              DO KEYT = 1,MAXY
                 ICOUNT1 = DOCBUF((KEYT - 1) * MAXX + 1)
                 IF (ICOUNT1 .NE. 0 .AND. KEYT .LE. MAXY2) THEN
C                   KEY EXISTS IN FIRST FILE, CHECK EXISTANCE IN 2ND
                     
                    ICOUNT2 = DOCBUF2(1,KEYT)
                    IF (ICOUNT2 .GT. 0) THEN
C                      NEED TO KEEP THIS KEY
                       DO IREG = 1,MAXX
                          DLIST(IREG) = DOCBUF((KEYT - 1) * MAXX + IREG)
                       ENDDO

C                      PUSH DLIST INTO DOC. FILE
                       NVAL  = MAXX - 1
                       CALL LUNDOCWRTDAT(NICDOCOUT,KEYT,DLIST(2),
     &                                NVAL,IRTFLG)
                    ENDIF

                 ELSEIF (KEY1 .NE. 0) THEN
C                   KEY EXISTS IN FIRST FILE, 2ND FILE FINISHED, QUIT
                    EXIT
                 ENDIF
              ENDDO

           ELSE                    
C             'ANDING' BY COLUMN OTHER THAN KEY
C             POINT TO NEXT VALUE IN SORTED LIST FROM FILE 2
              IGO2     = 1
              VALNEXT2 = SORTED2(IGO2)
              VALSENT  = MIN(SORTED(1),SORTED2(1)) - 1.0
              IKEYNOW  = 0

C             FIND OUTPUT VALUES
              DO IK = 1,IKEYS
C                FIND VALUE IN KEYCOL OF FIRST DOC. FILE

C                FIND KEY FOR THIS SORTED VALUE FROM FIRST FILE
                 KEY1 = SORTED(IK)

C                FIND SORTED VALUE FROM FIRST FILE
                 VAL1 = DOCBUF((KEY1 - 1) * MAXX + KEYCOL + 1)

                 SENDIT = .FALSE.
                 IF (IGO2 .GT. IKEYS2) THEN
C                   VALUE DOES NOT EXIST IN FILE 2 (ALL FILE 2'S DONE)
                    EXIT

                 ELSE
C                   MUST CHECK TO SEE IF VALUE EXISTS IN FILE 2
                    IF (VAL1 .LT. VALNEXT2) THEN
C                      VALUE 1 IS BELOW NEXT VALUE 2, DO NOT SAVE
                       SENDIT = .FALSE.

                    ELSEIF (VAL1 .EQ. VALNEXT2 .AND.
     &                      IGO2 .LE. IKEYS2) THEN
C                      VALUE 1 = NEXT VALUE 2, MUST INCREASE IGO2
                       SENDIT   = VALSENT .NE. VAL1
                       IGO2     = IGO2 + 1
                       IF (IGO2 .LE. IKEYS2) VALNEXT2 = SORTED2(IGO2)

                    ELSEIF (VAL1 .GT. VALNEXT2 .AND.
     &                      IGO2 .LT. IKEYS2) THEN
C                      VALUE 1 > NEXT VALUE 2, MUST INCREASE VALNEXT2
                       SENDIT = .FALSE.
                       DO IGO2T = IGO2+1,IKEYS2
                          IF (VAL1 .EQ. SORTED2(IGO2T)) THEN
C                            KEEP GOING TILL VALUE 2 > 
                             SENDIT = .TRUE.
                          ELSEIF (VAL1 .LE. SORTED2(IGO2T)) THEN
C                            POINT TO THIS VALUES IN SORTED LIST FROM #2
                             IGO2     = IGO2T
                             VALNEXT2 = SORTED2(IGO2)
                             EXIT
                         ENDIF
                       ENDDO
                    ENDIF                    
                 ENDIF

                 IF (SENDIT) THEN
C                   VALUE FROM KEYCOL ALSO IN 2'ND DOC. FILE
C                   PUT LINE FROM FIRST DOC. FILE INTO OUTPUT FILE
                    DO IREG = 1,MAXX
                       DLIST(IREG) = DOCBUF((KEY1 - 1) * MAXX + IREG)
                    ENDDO

C                   PUSH DLIST INTO DOC. FILE
                    NVAL    = MAXX - 1
                    IKEYNOW = IKEYNOW + 1
                    IF (KEYCOL .GT. 0) IKEYNOW = KEY1
                    CALL LUNDOCWRTDAT(NICDOCOUT,IKEYNOW,DLIST(2),
     &                                NVAL,IRTFLG)
                    VALSENT = DLIST(KEYCOL+1)
                 ENDIF
              ENDDO
           ENDIF
        ENDIF

C       CLOSE THE OUTPUT DOC. FILE(S)
9990    CLOSE(NDOCOUT)
        CLOSE(NDOCOUT2)

C       DEALLOCATE ALLOCATABLE ARRAYS
        IF (ALLOCATED(SORTED2))  DEALLOCATE(SORTED2)
        IF (ALLOCATED(SORTED))   DEALLOCATE(SORTED)
        IF (ASSOCIATED(DOCBUF2)) DEALLOCATE(DOCBUF2)
        CLOSE(NDOCIN2)

9999    RETURN
	END


