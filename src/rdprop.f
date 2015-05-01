
C++*********************************************************************
C
C  RDPROP.F -- DERIVED FROM RDPR.F                DEC 2005 ARDEAN LEITH
C              GLOBAL REGISTERS                   JAN 2006 ARDEAN LEITH
C              SKIP {var] REG LOAD                APR 2006 ARDEAN LEITH
C              SKIP SYM {var] SUB                 JUN 2009 ARDEAN LEITH
C              SKIP COMMAND LINE NON VAR INIT     AUG 2009 ARDEAN LEITH
C              @*** VAR INIT BUG                  NOV 2009 ARDEAN LEITH
C              NECHO                              SEP 2012 ARDEAN LEITH
C              ! COMMENT BUG                      DEC 2013 ARDEAN LEITH
C              NO LOAD [] FOR @ LINE              MAY 2014 ARDEAN LEITH
C              REG_GET_SEL aCESSES GLOBAL BUG     APR 2015 ARDEAN LEITH
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
C  RDPROP(PROMPT,ANSRET,NCHAR,IRTFLG)
C
C  PURPOSE: OUTPUTS PROMPT
C           READS AN OPEATION STRING FROM STORED PROC. LINE OR TERMINAL
C           ECHO & SKIP LINES WHICH ONLY CONTAIN A COMMENT AT START OF LINE
C           HANDLES INTERACTIVE HELP
C           CAN ECHO LINE TO CURRENT INTERACTIVE DO-LOOP FILE.
C           CONVERTS OLD @B01[X11] PROC. ARG. FORMAT TO TO NEW: () ARG.  
C           CONVERTS OLD <> VARIABLE FORMAT TO NEW [] VARIABLE FORMAT
C           CAN INVOKE VARIABLE SUBSTITUTION FOR [string]. 
C           CONVERTS OLD X REGISTER TO  [] VARIABLE FORMAT
C               
C           RETURNS NCHAR=LENGTH OF STRING WITHOUT TRAILING BLANKS OR COMMENT. 
C           COMMENT IS LIMITED TO 160 CHAR.
C           VARIABLE VALUE RESPONSE IS LIMITED TO 160 CHAR.
C             
C           REGISTER SUBSTITUTION OCCURS IN RDPRINC
C            
C  PARAMETERS:  PROMPT     INPUT PROMPT                     (SENT)
C               ANSRET     USER RESPONSE                    (RETURNED)
C               NCHAR      LAST NON_BLANK CHAR IN           (RETURNED)
C                            ANS RESPONSE BEFORE COMMENT
C               IRTFLG     RETURN FLAG (0 IS NORMAL)        (RETURNED)
C
C  CALLED BY:   SPIDER
C
C       IF THE FIRST CHARACTER IS '@', GOTO PROCEDURE EVALUATION
C       IF THE FIRST OR SECOND CHARACTER IS NEITHER A LETTER NOR DIGIT
C          CONSIDER OPERATION AN EXPRESSION
C       IF FIRST THREE CHARACTERS ARE LETTERS AND THE FORTH IS '('
C          THEN MUST BE AN ON-LINE FUNCTION CALL. GOTO EXPRESSION EVAL.
C       IF THE OPERATION STARTS WITH A [] REGISTER, EVALUATE EXPRESSION.
C       CHAR FOLLOWED BY 2 DIGITS IS OLD STYLE BATCH (B01) CALL
C       OK TO TRANSLATE OPERATION STRING TO UPPER CASE NOW
C       IF A LABEL 'LB<DIGIT> IS FOUND, AND A DO-LOOP IS IN EFFECT ...
C       IF A LABEL IS FOUND, AND NO DO-LOOP IS IN EFFECT ... IGNORE

C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE RDPROP(PROMPT,ANSRET,NCHAR,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      COMMON /LUNDOECHO/ LUNDONOW,NDOLINE

      CHARACTER(LEN=*)   :: PROMPT, ANSRET
      CHARACTER(LEN=161) :: ANS
      CHARACTER(LEN=161) :: COMMENTSTR
      CHARACTER(LEN=2 )  :: C2,NQ12

      CHARACTER(LEN=1 )  :: CTEMP
      LOGICAL            :: LOADR

      CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SETS ICOMM AND MYPID

      IDOL   = LEN(PROMPT)
      NQ12   = CHAR(34) // CHAR(39)   ! QUOTES
      IRTFLG = 0

10    CONTINUE
C     PROMPT OUTPUT IS DELAYED IN BATCH TO IGNORE COMMENT / BLANK LINES
      IF (COPT == 'I' .AND. MYPID <= 0) THEN
         WRITE(NOUT,90,ADVANCE='NO') PROMPT(1:IDOL)
 90      FORMAT(' .',A,': ')
      ENDIF

C     READ ANSWER STRING
      IF (NIN == 1) THEN
C        UPDATE THE PROCEDURE LINE COUNTER FOR CURRENT PROCEDURE LINE
         IBCNT = IBCNT + 1

C        READ FROM CURRENT STORED PROCEDURE LINE IBCNT
         CALL PROC_GETPLINE(IBCNT,0,ANS,NCHAR,IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(101,'PROCEDURE LACKS: RE',NE)
            ANSRET = 'RE'
            NCHAR  = 2
            RETURN
         ENDIF

      ELSE
C        UPDATE THE PROCEDURE LINE COUNTER FOR CURRENT PROCEDURE LINE
         IBCNT = IBCNT + 1

C        READ FROM TERMINAL
         READ(NIN,FMT='(A)',IOSTAT=IERR) ANS
         IF (IERR .NE. 0) THEN
            WRITE(NOUT,*) ' '
            CALL ERRT(101,'RAN OUT OF INPUT',NE)
            ANSRET = 'EN'
            NCHAR  = 2
            RETURN
         ENDIF
         NCHAR = lnblnk(ANS)

      ENDIF

C     IF A LOG FILE EXISTS, WRITE  ANSWER STRING THERE
      IF (NLOG .NE. 0 .AND. mypid <= 0) THEN
         IF (NCHAR <= 0) THEN
            WRITE(NLOG,*) ' '
         ELSE
            WRITE(NLOG,*) ANS(1:NCHAR)
            NECHO = NECHO + 1
         ENDIF
      ENDIF

C     IF WRITING TO RESULTS NOT TERMINAL, WRITE  ANSWER STRING THERE
      IF (NOUT .NE. 3 .AND. NDAT == 3 .AND. mypid <= 0) THEN
         WRITE(NOUT,*) ANS(1:NCHAR)
      ENDIF

      IF (NTRACE > 0 .AND. mypid <= 0) THEN
C        WRITES OUT OPERATION IN TRACE MODE
         IF (NCHAR <= 0) THEN
            WRITE(6,*) ' '
         ELSE
            WRITE(6,*) ANS(1:NCHAR)
         ENDIF
      ENDIF

      IF (NCHAR <= 0) GOTO 10   ! SKIP BLANK LINE

      COMMENTSTR = CHAR(0)

C     REMOVE LEADING AND TRAILING NON-PRINTING CHAR. FROM ANS
C     ALSO REMOVE ANY COMMENT FROM ANS STRING

      I      = 1
      J      = 0
      DO WHILE (I <= NCHAR)
         CTEMP = ANS(I:I)
         IF (C TEMP == ';' .OR. CTEMP == '!' ) THEN
            COMMENTSTR = ANS(I:)
            EXIT
         ELSEIF (J > 0 .OR. 
     &          (CTEMP > '!' .AND. CTEMP <= '~')) THEN
            J        = J + 1
            ANS(J:J) = ANS(I:I)
         ENDIF
         I = I + 1
      ENDDO
      NCHAR    = lnblnkn(ANS(1:J))
      NCHARCOM = lnblnkn(COMMENTSTR)

      IF (NCHAR <= 0) THEN
C        NOTHING BEFORE COMMENT, MUST GET NEXT LINE
         IF (VERBOSE .AND. MYPID <= 0 .AND. COPT == 'B') THEN
C           ECHO COMMENT 
C           (PROMPT OUTPUT WAS DELAYED IN BATCH TO IGNORE BLANK LINES)
            WRITE(NOUT,90)  PROMPT(1:IDOL),COMMENTSTR(1:NCHARCOM)
 91         FORMAT(' .',A,': ',A)
         ENDIF
         GOTO 10       ! READ ANOTHER INPUT LINE
      ENDIF
      IF (NCHAR < 161) ANS(NCHAR+1:) = CHAR(0)

C     PROMPT OUTPUT DELAYED UNTIL HERE IN BATCH TO IGNORE BLANK LINES 
      IF (COPT == 'B' .AND. MYPID <= 0) THEN
         WRITE(NOUT,94) PROMPT(1:IDOL),ANS(1:NCHAR),
     &                  COMMENTSTR(1:NCHARCOM)
94       FORMAT(' .',A,': 'A,'  ',A)
      ENDIF

C     HANDLE OBSOLETE INTERACTIVE HELP
      IQUES = INDEX(ANS(:NCHAR),'?')
      IF ( COPT == 'I'  .AND. (IQUES > 0  .OR. 
     &     (INDEX(ANS(:NCHAR),'HELP') > 0) .OR.
     &     (INDEX(ANS(:NCHAR),'help') > 0)) .AND.
     &     (INDEX(ANS(IQUES+1:NCHAR),'>')== 0)) THEN
         WRITE(NOUT,*)' USE YOUR WEB BROWSER FOR SPIDER MANUAL PAGES'
C        READ ANOTHER INPUT LINE
         GOTO 10
      ENDIF

      IF (LUNDONOW > 0 .AND. MYPID <= 0) THEN
C        MUST COPY INPUT LINE TO CURRENT INTERACTIVE DO-LOOP FILE
         WRITE(LUNDONOW,*) ANS(1:NCHAR)
         NDOLINE = NDOLINE + 1
         !write(6,*) ' rdprop lundo: ',ANS(:NCHAR),':',NDOLINE !!!!
      ELSEIF ((COPT == 'I' .AND. NDOLINE > 0)) THEN
         NDOLINE = NDOLINE + 1
         !write(6,*) ' rdprop ndoline: ',ANS(:NCHAR),':',NDOLINE !!!!
      ENDIF

C     CONVERT OLD @b01[..] ARG. FORMAT TO TO NEW: () ARG.  FORMAT
      IPAT   = INDEX(ANS(1:NCHAR), '@') 
      IQUO   = SCAN(ANS(1:NCHAR), NQ12) 
      IF (IPAT == 1 .AND. IQUO == 0) THEN
         IGOBRAK  = SCAN(ANS(1:NCHAR), '([')         ! ( first is ok
         IF (IGOBRAK > IPAT .AND. ANS(IGOBRAK:IGOBRAK) == '[')THEN
C           PROBABLY USES OLD @b01[] PROC. ARG. FORMAT 
            IENDBRAK = INDEX(ANS(1:NCHAR), ']',.TRUE.)
            CALL GETNEXTTOKEN2(ANS(IGOBRAK:IENDBRAK),1,IGO,IEND)
            INOT = VERIFY(ANS(IGOBRAK+IGO-1:IEND+IGOBRAK-1), 
     &                    '[xX0123456789]')
            IF (INOT == 0) THEN
C              PROBABLY  @b01[x11...,], [[var]...
               ANS(IGOBRAK:IGOBRAK)   = '(' 
               ANS(IENDBRAK:IENDBRAK) = ')'
            ENDIF
         ENDIF
      ENDIF
            
C     SEE IF NEED TO CONVERT OLD x11 REGISTER FORMAT TO [_11]
      IX = SCAN(ANS(1:NCHAR),'xX')
      IF (IX > 0) THEN
C        CONVERT OLD x11 REGISTER FORMAT TO TO NEW: [name] FORMAT
         CALL DEXREG(ANS,NCHAR)
      ENDIF
 
C     SEE IF ANY [] REGISTER VARIABLES ON OPERATION LINE
C     DO NOT LOAD REGISTERS IF: JUST A REGISTER OPERATION STARTING: [var]
C     DO NOT LOAD REGISTERS IF: 'IF(...)  OPERATION 
      IGOBRAK = INDEX(ANS(1:NCHAR),'[')     ! HAS REGISTERS
      IGOEQ   = INDEX(ANS(1:NCHAR), '=')    ! FOR ASSIGNMENT [var]=
      IGOAT   = INDEX(ANS(1:NCHAR), '@')    ! FOR PROCEDURE CALL

      C2    = ANS(1:2)
      CALL SSUPCAS(C2)
      LOADR = ((C2 == 'IF') .AND. IGOEQ > 2)
      IF (IGOBRAK > 1 .AND. IGOAT == 0 .AND. 
     &                     (IGOEQ <= 0 .OR. LOADR)) THEN
C        LOAD [] REGISTER VARIABLES IN NSEL LIST, ZERO UNUSED REGISTERS
         !write(6,*) 'Loaded nsel',igobrak,igoeq,ans(1:nchar)

         IENDBRAK = INDEX(ANS(1:NCHAR),']',.TRUE.)
         IF (IENDBRAK > 0) THEN
            IBANK = 0
            IF (ANS(1:4) == 'RR G') IBANK = 1
            CALL REG_GET_SEL(IBANK,ANS(IGOBRAK:IENDBRAK),
     &                       .TRUE.,.TRUE.,NREG,IRTFLG)
         ENDIF
      ELSE
         !if(igobrak>1)write(6,*)'Skip nsel load',igobrak,igoeq,ans(:nchar)
         CALL REG_SET_USED(0)
      ENDIF
 

C     SET RETURNED ANSWER, TRUNCATE TO FIT LENGTH OF ANS IN CALL
      LENA = LEN(ANSRET)
      IF (NCHAR > LENA) NCHAR = LENA
      IF (NCHAR > 0)    ANSRET(1:NCHAR) = ANS(1:NCHAR)

      IRTFLG = 0
      IF (ANS(1:1) == '^' .AND. NCHAR == 1) IRTFLG = -1

      END


C     *********************** SSUPCAS_NOVAR *****************************

      SUBROUTINE SSUPCAS_NOVAR(STRING)

      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=1) :: CTEMP

      DATA IOFF/-32/

      ILEN   = LEN(STRING)
      INBRAK = 0
  
      DO  I=1,ILEN
        CTEMP = STRING(I:I)
        IF (CTEMP == '[') THEN
           INBRAK = INBRAK + 1
        ELSEIF (CTEMP == ']') THEN
           INBRAK = INBRAK -1
        ELSEIF (INBRAK <= 0 .AND. 
     &         (STRING(I:I) >= 'a' .AND. STRING(I:I) <= 'z')) THEN
          STRING(I:I) = CHAR(ICHAR(CTEMP) + IOFF)
        ENDIF
      ENDDO

      RETURN
      END


C      *********************** DE_REG ********************************

       SUBROUTINE DE_REG(CINPUT,NCHAR)

      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*), INTENT(INOUT) :: CINPUT
      INTEGER, INTENT(INOUT)          :: NCHAR
      CHARACTER(LEN=161)              :: CSUB

C        CONVERT [_11] NEW REGISTER FORMAT TO old x11 FORMAT
         I     = 1
         DO WHILE (I < (NCHAR-2))
            IF (CINPUT(I:I+1) == '[_') THEN
C              PROBABLE REGISTER START OF [_
               NDIG = VERIFY(CINPUT(I+2:NCHAR),'0123456789')
               IF (NDIG > 0) THEN
                  NDIG = NDIG - 1
               ELSE
                  NDIG = NCHAR - I
               ENDIF

               IF (NDIG > 0) THEN
                  CSUB = 'x' // CINPUT(I+2:I+NDIG+1) // CHAR(0)
                  CALL SUBCHAR(CSUB(1:NDIG+1),CINPUT,I,I+NDIG+2,
     &                          NCHAR,IRTFLG)
                  I = I + NDIG 
               ENDIF 
            ENDIF
            I = I + 1
         ENDDO
c        write(6,*) ' after de_reg: ',CINPUT(1:NCHAR)
         END





