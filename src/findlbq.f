
C++*********************************************************************
C
C  FINDLBQ.F -- ADAPTED FOR CHAR.                  AUG  89 ArDean Leith
C               MERGED WITH SEARCHQ STUFF          SEPT 97 ArDean Leith
C               ADDED IFLEVEL                      DEC. 97 ArDean Leith
C               INCORE PROCS                       JAN  01 ArDean Leith
C               LUNDONOW                           FEB  01 ArDean Leith
C               LNBLNKN                            MAY  04 ArDean Leith
C               DOC INSIDE LOOP BUG                JUL  07 ArDean Leith
C               NINSAVEOF                          NOV  09 ArDean Leith
C               PASS LOOP END WHILE HUNTING LB     DEC  09 Ardean Leith
C               PASS LB6 WHILE HUNTING LB7 BUG     JAN  10 Ardean Leith
C               IF CYCLE WITHIN DO BUG             JAN  10 Ardean Leith
C               ! COMMENT DELIMITER                DEC  11 Ardean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C  FINDLBQ(WANTLABEL,IDOTOP,NLOOPT,IDOSTK,NEWLOOP, IFLEVEL,IRTFLG)
C
C  PURPOSE: SEARCHES FOR 'LB*', OR 'ENDDO' IN SPIDER INCOMING
C           OPERATIONS. IF SUCCESSFUL RETURNS CURRENT LINE IN: FCHAR
C
C  PARAMETERS:
C        WANTLABEL      LABEL WE ARE SEARCHING FOR                (SENT)
C                          CAN BE LB#, ELSE, ENDIF, ENDDO
C        IDOTOP         DO-LOOP NESTING LEVEL                     (SENT)
C        NLOOPT         INSIDE DO-LOOP NOW IF > 0                 (SENT)
C	 IDOSTK         DO-LOOP STACK                             (SENT)
C	 NEWLOOP	IF PASS CURRENT A DO-LOOP END -            (RET)
C                          MUST POP DO-LOOP STACK IN CALLER         
C        IFLEVEL        IF CLAUSE NESTING LEVEL            (SENT & RET.) 
C        IRTFLG		ERROR NUMBER ZERO IS NORMAL               (RET.)
C
C--*******************************************************************

	SUBROUTINE FINDLBQ(WANTLABEL,IDOTOP,NLOOPT,IDOSTK,NEWLOOP,
     &                     IFLEVEL,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        COMMON /LUNDOECHO/ LUNDONOW,NDOLINE,NINSAVEOF

        CHARACTER (LEN=*)      :: WANTLABEL
        CHARACTER (LEN=MAXNAM) :: ANSW
        LOGICAL                :: KEEPGO
        INTEGER                :: IDOSTK(7,*)
        INTEGER                :: DOLABEL
        INTEGER                :: NEWLOOP 

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

	NEWLOOP    = 0
        KEEPGO     = .TRUE.
        NEEDENDIF  = 0
        NEEDENDO   = 1

C       KEEP READING AND DISCARDING LINES FROM INPUT AS NEEDED
        DO WHILE (KEEPGO) 
           ANSW(4:4) = ' '

           IF (COPT .EQ. 'B' .AND. NINSAVEOF .EQ. 0) THEN
C             INCREMENT PROCEDURE LINE COUNTER
              IBCNT = IBCNT + 1

C             READ FROM CURRENT STORED PROCEDURE LINE= IBCNT
              CALL PROC_GETPLINE(IBCNT,0,ANSW,NCHAR,IOS)
           ELSE
C             READ FROM FILE OPENED ON NIN (FOR INTERACTIVE LOOP)
              READ(NIN,FMT='(A)',IOSTAT=IOS) ANSW
              NCHAR = LNBLNKN(ANSW)
              IF (LUNDONOW .GT. 0) THEN
C               MUST COPY INPUT LINE TO CURRENT INTERACTIVE DO-LOOP FILE
                WRITE(LUNDONOW,*) ANSW(1:NCHAR)
              ENDIF
           ENDIF

           IF (IOS .NE. 0) THEN
              WRITE(NOUT,91) WANTLABEL
 91           FORMAT(' *** JUMP DESTINATION NEVER FOUND: ',A)
              CALL ERRT(100,'FINDLBQ',NE)
              IRTFLG = 1
              RETURN
           ENDIF

C          REMOVE BLANKS FROM INPUT STRING
           CALL SHRINKQ(ANSW,MAXNAM,ANSW,NCHAR)

C          IGNORE ANY COMMENTS AT END OF INPUT STRING
           ISEMICOL = SCAN(ANSW(1:NCHAR),';!')
           IF (ISEMICOL .GT. 0) THEN
              ANSW(ISEMICOL:NCHAR) = CHAR(0)
              NCHAR = ISEMICOL - 1
           ENDIF
           IF (NCHAR .LE. 0) CYCLE

C          CONVERT INPUT STRING TO ALL UPPER CASE
           CALL SSUPCAS(ANSW(1:NCHAR))

C          SEE IF INPUT LINE CONTAINS 'THEN'
           IIFTHEN = INDEX(ANSW(1:NCHAR),'THEN')

C          FIND LAST NON-BLANK IN WANTLABEL            
           LENLB = LNBLNKN(WANTLABEL)

C          FIND LAST NON-BLANK IN ANSW            
           LENANSW = LNBLNKN(ANSW)
           LENT    = INDEX(ANSW,CHAR(0)) - 1
           IF (LENT .GT. 0) LENANSW = MIN(LENANSW,LENT)

        !write(6,*)' want: ',wantlabel(:lenlb),' answ: ',answ(1:lenansw)

C          -------------------------------------- SEARCHING FOR "ENDDO"
C          DO NOT CARE ABOUT IF'S, OR LABELS IN HERE

          IF (WANTLABEL(:LENLB) .EQ. 'ENDDO' .AND.
     &         ANSW(1:2)         .EQ. 'DO'   .AND.
     &         ANSW(3:3)         .NE. 'C'    .AND.
     &         ANSW(3:4)         .NE. 'LB') THEN

C              HUNTING FOR 'ENDDO' AND PASSING NESTED 'DO...ENDDO'
               NEEDENDO = NEEDENDO + 1
               !write(6,*) ' needendop1:',needendo

           ELSEIF (WANTLABEL(:LENLB) .EQ. 'ENDDO' .AND.
     &             ANSW(1:2)      .EQ. 'IF' .AND.
     &             IIFTHEN .GT. 5 ) THEN
C              HUNTING FOR 'ENDDO' AND FOUND A NEW 'IF...THEN' CLAUSE
               NEEDENDIF = NEEDENDIF + 1

           ELSEIF (WANTLABEL(:LENLB) .EQ. 'ENDDO' .AND.
     &             ANSW(1:5)      .EQ. 'ENDIF') THEN
C              HUNTING FOR 'ENDDO' AND FOUND AN 'ENDIF' 
               NEEDENDIF = NEEDENDIF - 1

           ELSEIF (WANTLABEL(:LENLB) .EQ. 'ENDDO' .AND.
     &             ANSW(1:5)         .EQ. 'ENDDO') THEN

C              HUNTING FOR 'ENDDO' AND FOUND IT
               NEEDENDO = NEEDENDO - 1
               !write(6,*) ' needendom1:',needendo

C              IF FOUND CORRESPONDING ENDDO, HALT INPUT
               IF (NEEDENDO .LE. 0) THEN
                  KEEPGO  = .FALSE.
                  IFLEVEL = IFLEVEL + NEEDENDIF
                  IF (NEEDENDIF .GT. 0) THEN
                     WRITE(NOUT,*)' WARNING: JUMP INTO "IF" CLAUSE'
                  ENDIF
              ENDIF

C          ----------------------------------------- SEARCHING FOR "LB"
C          MUST CARE ABOUT IF'S, DO's, 'LB**'s, AND ENDDO's IN HERE

           ELSEIF (ANSW(1:2) .EQ. 'LB' .AND.
     &         ANSW(:LENANSW) .EQ. WANTLABEL(:LENLB)) THEN
C              HUNTING FOR SPECIFIED 'LB??' AND FOUND IT, HALT INPUT
               KEEPGO  = .FALSE.

               IFLEVEL = IFLEVEL + NEEDENDIF
               IF (NEEDENDIF .GT. 0) THEN
                  WRITE(NOUT,*)' WARNING: JUMP INTO "IF" CLAUSE'
               ENDIF
               !write(6,*) ' found:',ANSW(:LENANSW),':',WANTLABEL(:LENLB),':'

           ELSEIF (WANTLABEL(1:2) .EQ. 'LB' .AND.
     &             ANSW(1:2)      .EQ. 'IF' .AND.
     &             IIFTHEN .GT. 5 ) THEN
C              HUNTING FOR 'LB..' AND FOUND A NEW 'IF...THEN' CLAUSE
               NEEDENDIF = NEEDENDIF + 1

           ELSEIF (WANTLABEL(1:2) .EQ. 'LB' .AND.
     &             ANSW(1:5)      .EQ. 'ENDIF') THEN
C              HUNTING FOR 'LB..' AND FOUND A 'ENDIF' 
               NEEDENDIF = NEEDENDIF - 1

           ELSEIF (WANTLABEL(1:2) .EQ. 'LB' .AND.
     &             ANSW(1:2)      .EQ. 'DO' .AND.
     &             ANSW(3:3)      .NE. 'C'  .AND.
     &             ANSW(3:4)      .NE. 'LB') THEN
C              HUNTING FOR 'LB..' AND FOUND A NEW UNLABELED 'DO' CLAUSE

               NEEDENDO = NEEDENDO + 1
               !write(6,*) ' needendop1:',needendo

           ELSEIF (WANTLABEL(1:2) .EQ. 'LB' .AND.
     &             ANSW(1:5)      .EQ. 'ENDDO') THEN
C              HUNTING FOR 'LB..' AND PASSING A "ENDDO" CLAUSE

               NEEDENDO = NEEDENDO - 1
               !write(6,*) ' needendom1:',needendo

           ENDIF

C          -------------------------------- PASS END OF ACTIVE DO-LOOP?

C          SET NEWLOOP IF WE PASS BY CURRENT DO-LOOP LABEL
           !write(6,*) ' keepgo,nloopt,newloop:',keepgo,nloopt,newloop

	   IF (KEEPGO .AND. NLOOPT .GT. 0) THEN 
              IF (ANSW(1:2) .EQ. 'LB') THEN
C                PASSED BY LINE LABEL, MAY BE END OF A LOOP
                 CALL GETLBNO(ANSW,LBNO,IRTFLG)
                 DO I = IDOTOP,1,-1
                    IF (LBNO .EQ. IDOSTK(6,I)) THEN
C                      FOUND END OF ACTIVE DO-LOOP
                       NEWLOOP = NEWLOOP - 1
                       !write(6,*) ' LBNO,IDOSTK(6,I):',LBNO,IDOSTK(6,I)
                    ENDIF
                 ENDDO
                 !write(6,990)newloop,idotop,(idostk(6,i),i=1,idotop)
990              format(' newloop:',i4,' idotop:',i4,' stk:',10i5)

              ELSEIF (ANSW(1:5) .EQ. 'ENDDO') THEN
C                PASSED BY ENDDO, END OF A LOOP

                 IF (NEEDENDO .LE. 0) THEN
C                   MAY BE A CURRENT ACTIVE DO-LOOP
                    IF (IDOSTK(6,IDOTOP).LT.0) THEN
C                       USING LB LESS DO LOOP SYNTAX
                        NEWLOOP = NEWLOOP - 1
                    !write(6,990)newloop,idotop,(idostk(6,i),i=1,idotop)
                    ENDIF
                  ENDIF
              ENDIF
           ENDIF
        ENDDO

C       FOUND DESIRED LABEL OR ENDIF CAN RETURN NOW
        IRTFLG = 0

C       SIMULATE ECHO OF OPERATION TO RESULTS FILE
        FCHAR = ANSW

        IF (MYPID .LE. 0) WRITE(NDAT,90) ANSW(1:NCHAR)
 90     FORMAT(' .OPERATION:',5X,A)

        RETURN
	END




C++*********************************************************************
C
C  FINDENDIF -- ADAPTED FOR CHAR.                  AUG  89 ArDean Leith
C               MERGED WITH SEARCHQ STUFF          SEPT 97 ArDean Leith
C               ADDED IFLEVEL                      DEC. 97 ArDean Leith
C               INCORE PROCS                       JAN  01 ArDean Leith
C               LUNDONOW                           FEB  01 ArDean Leith
C               LNBLNKN                            MAY  04 ArDean Leith
C               DOC INSIDE LOOP BUG                JUL  07 ArDean Leith
C               NINSAVEOF                          NOV  09 ArDean Leith
C               PASS LOOP END WHILE HUNTING LB     DEC  09 Ardean Leith
C **********************************************************************
C **********************************************************************
C
C  FINDENDIF(WANTLABEL,IFLEVEL,IRTFLG)
C
C  PURPOSE: SEARCHES FOR LB*, ELSE, ENDIF, ENDDO IN SPIDER INCOMING
C           OPERATIONS. IF SUCCESSFUL RETURNS CURRENT LINE IN: FCHAR
C
C  PARAMETERS:
C        WANTLABEL      LABEL WE ARE SEARCHING FOR                (SENT)
C                          CAN BE 'ELSE' OR 'ENDIF'
C        IFLEVEL        IF CLAUSE NESTING LEVEL            (SENT & RET.) 
C        IRTFLG		ERROR RETURN, ZERO IS NORMAL              (RET.)
C
C--*******************************************************************

	SUBROUTINE FINDENDIF(WANTLABEL,IFLEVEL,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        COMMON /LUNDOECHO/ LUNDONOW,NDOLINE,NINSAVEOF

        CHARACTER (LEN=*)      :: WANTLABEL
        CHARACTER (LEN=MAXNAM) :: ANSW
        LOGICAL                :: KEEPGO

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        KEEPGO     = .TRUE.
        NEEDENDIF  = 0

C       KEEP READING AND DISCARDING LINES FROM INPUT AS NEEDED
        DO WHILE (KEEPGO) 
           ANSW(4:4) = ' '

           IF (COPT .EQ. 'B' .AND. NINSAVEOF .EQ. 0) THEN
C             INCREMENT PROCEDURE LINE COUNTER
              IBCNT = IBCNT + 1

C             READ FROM CURRENT STORED PROCEDURE LINE= IBCNT
              CALL PROC_GETPLINE(IBCNT,0,ANSW,NCHAR,IOS)
           ELSE
C             READ FROM FILE OPENED ON NIN (FOR INTERACTIVE LOOP)
              READ(NIN,FMT='(A)',IOSTAT=IOS) ANSW
              NCHAR = LNBLNKN(ANSW)
              IF (LUNDONOW .GT. 0) THEN
C               MUST COPY INPUT LINE TO CURRENT INTERACTIVE DO-LOOP FILE
                WRITE(LUNDONOW,*) ANSW(1:NCHAR)
              ENDIF
           ENDIF

           IF (IOS .NE. 0) THEN
              WRITE(NOUT,91) WANTLABEL
 91           FORMAT(' *** DESTINATION NEVER FOUND: ',A)
              CALL ERRT(100,'FINDENDIF',NE)
              IRTFLG = 1
              RETURN
           ENDIF

C          REMOVE BLANKS FROM INPUT STRING
           CALL SHRINKQ(ANSW,MAXNAM,ANSW,NCHAR)

C          IGNORE ANY COMMENTS AT END OF INPUT STRING
           ISEMICOL = SCAN(ANSW(1:NCHAR),';!')
           IF (ISEMICOL .GT. 0) THEN
              ANSW(ISEMICOL:NCHAR) = CHAR(0)
              NCHAR = ISEMICOL - 1
           ENDIF
           IF (NCHAR .LE. 0) CYCLE

C          CONVERT INPUT STRING TO ALL UPPER CASE
           CALL SSUPCAS(ANSW(1:NCHAR))

C          SEE IF INPUT LINE CONTAINS 'THEN'
           IIFTHEN = INDEX(ANSW(1:NCHAR),'THEN')

C          FIND LAST NON-BLANK IN WANTLABEL            
           LENLB = LNBLNKN(WANTLABEL)

C          FIND LAST NON-BLANK IN ANSW            
           LENANSW = LNBLNKN(ANSW)
           LENT    = INDEX(ANSW,CHAR(0)) - 1
           IF (LENT .GT. 0) LENANSW = MIN(LENANSW,LENT)

C          -------------------------------------- SEARCHING FOR "ELSE"

           IF (WANTLABEL .EQ. 'ELSE' .AND.
     &             ANSW(1:6) .EQ. 'ELSEIF') THEN
C              HUNTING FOR 'ELSE' AND FOUND 'ELSEIF'

C              IF NOT A NESTED ELSE, RETURN TO CALLER
               IF (NEEDENDIF .LE. 0)  THEN
                   KEEPGO = .FALSE.

C                  DECREMENT PROCEDURE LINE COUNTER TO RE-READ THIS LINE
                   IF (COPT .EQ. 'B' .AND. NINSAVEOF .EQ. 0)
     &                 IBCNT = IBCNT-1
               ENDIF

           ELSEIF (WANTLABEL .EQ. 'ELSE' .AND.
     &             ANSW(1:4) .EQ. 'ELSE') THEN
C              HUNTING FOR 'ELSE' AND FOUND IT
C              IF NOT A NESTED ELSE, RETURN TO CALLER
               IF (NEEDENDIF .LE. 0)  KEEPGO = .FALSE.

           ELSEIF (WANTLABEL .EQ. 'ELSE' .AND.
     &             ANSW(1:5) .EQ. 'ENDIF') THEN
C              HUNTING FOR 'ELSE' AND FOUND ENDIF

C              IF NOT A NESTED ENDIF, RETURN TO CALLER
               IF (NEEDENDIF .LE. 0) THEN
                  IFLEVEL = IFLEVEL - 1
                  KEEPGO  = .FALSE.
               ENDIF

C              IF NESTED ENDIF KEEP READING INPUT
               NEEDENDIF = NEEDENDIF - 1

           ELSEIF (WANTLABEL .EQ. 'ELSE' .AND.
     &             ANSW(1:2) .EQ. 'IF'   .AND.
     &             IIFTHEN .GT. 5 ) THEN
C              HUNTING FOR 'ELSE' AND FOUND IF...THEN
C              THIS IS A NESTED IF,  KEEP READING INPUT
               NEEDENDIF = NEEDENDIF + 1

C          -------------------------------------- SEARCHING FOR "ENDIF"
C          DO NOT CARE ABOUT NEW LOOPS, OR LABELS IN HERE

           ELSEIF (WANTLABEL .EQ. 'ENDIF' .AND.
     &             ANSW(1:2) .EQ. 'IF'    .AND.
     &            IIFTHEN .GT. 5 ) THEN
C              HUNTING FOR 'ENDIF' AND FOUND A NESTED IF...THEN
               NEEDENDIF = NEEDENDIF + 1

           ELSEIF (WANTLABEL .EQ. 'ENDIF' .AND.
     &             ANSW(1:5) .EQ. 'ENDIF') THEN
C              HUNTING FOR 'ENDIF' AND FOUND ENDIF

C              IF NOT A NESTED ENDIF, RETURN TO CALLER
               IF (NEEDENDIF .LE. 0) THEN
                  KEEPGO  = .FALSE.
                  IFLEVEL = IFLEVEL - 1
               ENDIF

C              IF NESTED ENDIF KEEP READING INPUT
               NEEDENDIF = NEEDENDIF - 1

           ENDIF
        ENDDO

C       FOUND DESIRED 'ELSE' OR 'ENDIF' CAN RETURN NOW
        IRTFLG = 0

C       SIMULATE ECHO OF OPERATION TO RESULTS FILE
        FCHAR = ANSW

        IF (MYPID .LE. 0) WRITE(NDAT,90) ANSW(1:NCHAR)
 90     FORMAT(' .OPERATION:',5X,A)

        RETURN
	END
