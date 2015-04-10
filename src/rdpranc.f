
C++*********************************************************************
C
C  RDPRANC.F  -- CREATED JULY 91 FROM RDPRAI               ARDEAN LEITH
C                SPLIT OUT FROM RDPRAI & RDPRAF  AUG. 1999 ARDEAN LEITH
C                FIXED OVERFLOW BUG AUG 99                 ARDEAN LEITH
C                FIXED FILLRANGE 1-1 BUG JAN 00            ARDEAN LEITH
C                USED GETNEXTOKEN FOR IBM JAN 00           ARDEAN LEITH
C                FIXED OVERFLOW BUG DEC 00                 ARDEAN LEITH
C                INTEGER RDPRANC  MAY 02                   ARDEAN LEITH
C                [VAR-VAR] BUG MAY 07                      ARDEAN LEITH
C
C **********************************************************************
C * AUTHOR:   ArDean Leith                                             *
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
C    RDPRANC(ILIST,NMAX,NUMBER,ILOW,IHI,PROMPT,LIMITED,IRTFLG)
C
C    PURPOSE:  READ AN ALPHANUMERIC STRING, CHECK FOR ANY SPECIAL OPERATION,
C              RETURNS ARRAY OF INTEGERS, THE NUMBER OF VALUES, OR A
C              FLAG  TO INDICATE THAT ONE SHOULD RETURN TO PREVIOUS QUESTION.
C
C              ALLOWABLE STRINGS           NUMBERS ENTERED IN ARRAY
C                 1.0,2,83.,66E1           1, 2, 8, 660
C                 1.0 2 83. 66E1           1, 2, 8, 660
C                 (CAN PUT IN A SERIES OF NUMBERS WITH -)
C                 1.0-5                    1, 2, 3, 4, 5
C                 X11,X12....              CONTENTS OF X11 & X12...
C                 X11-X12                  CONTENTS OF X11....X12
C
C  PARAMETERS : ILIST     ARRAY FOR ANSWERS                      (RET.)
C               NMAX      MAX LENGTH OF ARRAY                    (SENT)
C               NUMBER    ON ENTRY IS MAX NUMBER OF ANSWERS (SENT/RET.) 
C                         TO BE RETURNED!!
C                         (<0 ON ENTRY IS FLAG TO ACCEPT NULL 
C                         RESPONSE)
C                         ON RETURN IS NUMBER OF ANSWERS ACTUALLY 
C                         RETURNED
C               ILOW      BOTTEM OF RANGE (CAN'T BE BELOW THIS)   (SENT)
C               IHI       TOP OF RANGE (CAN'T BE ABOVE THIS)      (SENT)
C                         ILOW & IHI CURRENTLY UNUSED!!!!!!!!!!!!!!
C               PROMPT    SOLICITATION MESSAGE                    (SENT)
C               LIMITED   REQUIRES NMAX VALUES / CALL, SO THAT    (SENT)
C                           IT CAN BE USED IN BATCH DO-LOOP
C               IRTFLG    RETURN FLAG (0 IS NORMAL, -1 IS GOTO    (SENT)
C                              PREVIOUS QUESTION)
C      
C    CALLED BY:   
C
C    NOTE: CAN PERFORM REGISTER SUBSTITUTION NOW al
C          THIS DIFFERS FROM: RDPRA IN THAT IT DOES NOT EVALULATE
C          EXPRESSIONS AND INTERPRETS 1-10 AS A LIST: 1,2,3,..10
C          DO NOT CONFUSE THE TWO SUBROUTINES!!
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE RDPRANC(ILIST,NMAX,NUMBER,ILOW,IHI,PROMPT,
     &                     LIMITED,IRTFLG)        

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER            :: IBUF(NIMAX)
        INTEGER            :: ILIST(NMAX)

        CHARACTER *(*)     :: PROMPT

C       MAXANS IS LENGTH OF ANS
        PARAMETER            (MAXANS = 600)
        CHARACTER * 600    :: ANS

        CHARACTER * 80     :: ANST
	CHARACTER * 1      :: NULL,A
        LOGICAL            :: LOOPNOPAREN,LIMITED
        LOGICAL            :: RANGE,RFLAG

        MAXB = NIMAX
        NULL = CHAR(0)

C       DODO PROOF THE INPUT
        NVAL = ABS(NUMBER)
        IF (NVAL == 0)    NVAL = 1
        IF (NVAL > NMAX) NVAL = NMAX

        NUMBER    = 0
        RANGE  = .FALSE.

C       GET INPUT LINE(S)
20      IF (PROMPT(:1) == '~') THEN
C          USE PROMPT FOR INPUT LINE
           ANS   = PROMPT(2:)
           NCHAR = LEN(PROMPT) - 1
        ELSE
           IRTFLG = -999   ! DO NOT UPPERCASE IT
           CALL RDPRMC(ANS,NCHAR,.TRUE.,PROMPT,NULL,IRTFLG)
           IF (NCHAR == 0 .OR. IRTFLG .NE. 0) RETURN
        ENDIF

21     IF (ANS(NCHAR:NCHAR) == ',') THEN
C          INPUT CONTINUATION LINE 
           IRTFLG = -999   ! DO NOT UPPERCASE IT
           CALL RDPRMC(ANST,NCHAR2,.TRUE.,
     &          'NEXT LINE OF INPUT',NULL,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           NCHAR = NCHAR + NCHAR2
           IF ((NCHAR + NCHAR2) > MAXANS) THEN
C             ANS HAS ALREADY OVERFLOWED!!!
              CALL ERRT(102,'ANS OVERFLOW',MAXANS)
              GOTO 20
           ENDIF
           ANS(NCHAR+1:) = ANST(1:NCHAR2)
           NCHAR = NCHAR + NCHAR2

C          SEE IF THERE ARE MORE LINES OF INPUT
           GOTO 21
        ENDIF

        IRTFLG      = 1

        LOOPNOPAREN = (ANS(1:1) .NE. '(' .AND. 
     &                 NLOOP > 1 .AND.
     &                 ILOOP > 1 .AND.
     &                 LEGACYPAR)
C       SEE IF REGISTERS IN INPUT BEFORE ANY SEMICOLON
        ISEMICOLON  = INDEX(ANS,';')
        IF (ISEMICOLON == 0) ISEMICOLON = NCHAR
        LOOPNOPAREN = LOOPNOPAREN .AND. 
     &               (INDEX(ANS(1:ISEMICOLON),'[') == 0)

        IF (LOOPNOPAREN .AND. .NOT. LIMITED) THEN
C          WITHIN DO-LOOP BUT NO '()' AND INPUT NUMBER IS NOT LIMITED
C          CAN NOT HAVE: X,Y,Z,A.... INPUT WITHIN A DO LOOP
           WRITE(NOUT,96)
96         FORMAT(' *** ERROR, WITHIN A DO LOOP THIS QUESTION NEEDS ',
     &            '() AROUND ANSWER.')
           CALL ERRT(100,'RDPRANC',NDUM)
           RETURN

        ELSEIF (NCHAR == 1 .AND. ANS(1:1) == '*') THEN
           IRTFLG = -1
           RETURN
        ENDIF

        IEND     = 0
        IGORANGE = -1
        RANGE    = .FALSE.

        DO           ! ----------------- LOOP --------------------------
           IFIRST = IEND + 1

C          GET TOKEN (CHAR. STRING DELIMITED BY A ", ( ) ] -")
           CALL GETNEXTTOKEN_N(ANS,NCHAR,IFIRST,IGO,IEND)

C          SEE IF ALL TOKENS FROM STRING HAVE BEEN EVALUATED
           IF (IGO <= 0) EXIT

C           write(6,*) '  RDPRANC TOKEN: ',ANS(IGO:IEND)
C           write(6,*) '  IGO...IEND: ',IGO,IEND

C          TOKEN RETURNED, SET START POSITION OF NEXT TOKEN IN ANS
           IF (ANS(IGO:IEND) == '-') THEN
C             '-' IS START OF A RANGE DELIMITER
              RANGE = .TRUE.
              CYCLE
           ENDIF
 
C          EVALUATE NUMBER OR REGISTER
           IF (ANS(IGO:IGO) == '[')THEN
C             TOKEN STARTS WITH '[', SO IS A REGISTER REFERENCE
C             GET REGISTER VALUE & NUMBER FROM TOKEN
              CALL REG_GET_VAR(0,ANS(IGO:IEND),.FALSE.,FNUM,IREG,
     &                         LASTVAR,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ELSE
C             EXTRACT REAL NUMBER    
              READ(ANS(IGO:IEND),*,IOSTAT=IERR) FNUM
              IF (IERR .NE. 0) THEN
                 WRITE(NOUT,*) ' *** SYNTAX ERROR: ',ANS(1:IEND)
                 CALL ERRT(100,'RDPRANC',NDUM)
                 IF (PROMPT(1:1) == '~') RETURN
                 GOTO 20
              ENDIF
           ENDIF

           IENDRANGE = FNUM

           IF (RANGE .AND. IENDRANGE < IGORANGE) THEN
              WRITE(NOUT,97)IGORANGE,IENDRANGE
   97         FORMAT(' *** ERROR: INVALID RANGE: ',I10,'...',I10)
              CALL ERRT(100,'RDPRANC',NDUM)
              GOTO 20

           ELSEIF (RANGE .AND. IENDRANGE == IGORANGE) THEN
C             RANGE IS SAME AS IGORANGE, WHICH IS ALREADY HANDLED
              RANGE = .FALSE.
              CYCLE

           ELSE 
              IF (RANGE) THEN
C                FILL WITH A RANGE OF INTEGER VALUES, STARTED FROM IGORANGE 
                 IGORANGE = IGORANGE + 1
              ELSE
                 IGORANGE = IENDRANGE 
              ENDIF
                   
C             FILL IBUF OR ILIST
              DO IVAL = IGORANGE,IENDRANGE
                 NUMBER = NUMBER + 1
                 IF (LOOPNOPAREN) THEN
C                   FILL IBUF AS A TEMPORARY ARRAY
                    IF (NUMBER > MAXB) THEN
C                      IBUF WILL OVERFLOW ! 
                       CALL ERRT(102,'IBUF OVERFLOW',MAXB)
                       RETURN
                    ENDIF
                    IBUF(NUMBER)  = IVAL
                 ELSE
C                   FILL ILIST
                    IF (NUMBER > NMAX) THEN
C                      ILIST WILL OVERFLOW ! 
                       CALL ERRT(102,'LIST OVERFLOW',NMAX)
                       RETURN
                    ENDIF
                    ILIST(NUMBER) = IVAL
                 ENDIF
              ENDDO
              IGORANGE = IENDRANGE 
           ENDIF

C          GET ANY REMAINING NUMBERS FROM THE INPUT STRING
           IF (.NOT. LOOPNOPAREN .AND. NUMBER >= NVAL) EXIT
           RANGE = .FALSE.
        ENDDO

C------------------- END OF ANS PARSING LOOP -----------------------

        IF (NUMBER == 0) THEN
C          GOT NULL INPUT, ACCEPT IT, DO NOT ALTER ILIST
           IRTFLG = 0
           RETURN
        ENDIF

        IF (LOOPNOPAREN) THEN
C          INPUT HAS NO () AROUND IT, AND IS IN A DO-LOOP
C          USES DIFFERENT SET OF INPUTS FOR EACH INDEX OF CURRENT LOOP
 
C          ILOC IS POINTER TO CURRENT LOCATION IN IBUF
           ILOC = (ILOOP - 1) * NVAL + 1

C          NGOT IS NUMBER OF VALUES LEFT IN IBUF
           NGOT = NUMBER - ILOC + 1

           IF (NGOT < NVAL)  THEN
C             NEEDS MORE INPUT TO GET "NVAL" INPUTS, READ ANOTHER LINE.
              WRITE(NOUT,95) NVAL,ILOOP,NGOT
   95         FORMAT(' *** ERROR: NEEDS AT LEAST: ',I4,
     &          '  INPUT VALUES FOR LOOP: ',I4,' BUT ONLY GOT: ',I4/)
              CALL ERRT(100,'RDPRANC',NDUM)
              GOTO 20
           ENDIF

           DO I = 1,NVAL
              ILIST(I) = IBUF(I+ILOC-1)
           ENDDO
           NUMBER = NVAL
        ENDIF
C----------------------------------------------
        IRTFLG = 0

c         write (NOUT,*) ' number of files: ',number
c         DO I = 1,number
c             write (NOUT,*) number, ': ',ILIST(I)
c          ENDDO
c        write (NOUT,*) ' -----------'


        RETURN

        END
       
