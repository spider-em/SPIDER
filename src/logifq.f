

C++*********************************************************************
C                                                                      *
C LOGIFQ.F     ADAPTED   FOR CHAR. VARIABLES      AUG  89 ArDean Leith *
C              UNCONDITIONAL JUMP ADDED           SEPT 96 ArDean Leith *
C              IF (...) THEN IMPLEMENTED          SEPT 97 ArDean Leith *
C              EXIT ADDED                         NOV  06 ArDean Leith *
C              ARASW GLOBAL                       JUN  09 ArDean Leith *
C              [SYMVAR] SUPPORT                   APR  10 ArDean Leith *
C              <==>! SUPPORT                      DEC  11 ArDean Leith *
C             ! COMMENT DELIMITER                 DEC  11 ArDean Leith *
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
C    LOGIFQ(STRING,LABEL,JUMP,IER)
C
C    PURPOSE:     EVALUATES STRING SUCH AS: IF([x].LE.5) GOTO LB77
C                 OR: IF([x].LE.5) THEN
C                 RETURNS THE LB77, AND A LOGICAL FLAG WHETHER ONE 
C                 SHOULD FOLLOW THE GOTO TO THE LABEL/ELSE LOCATION.
C                 ALSO ACCEPTS STRING SUCH AS: IF ([X].LE.5) [y]=9
C                 AND EVALUATES THE SECOND EXPRESSION BEFORE
C                 RETURNING.
C
C    PARAMETERS:  STRING       INPUT LINE                  (SENT)
C                 LABEL        LABEL STRING                (RETURNED)
C                 JUMP         FLAG TO FOLLOW GOTO         (RETURNED)
C                 IERR         ERROR FLAG (1 IS ERROR)     (RETURNED)
C
C--*********************************************************************

      SUBROUTINE LOGIFQ(STRING,LABEL,JUMP,IER)

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'


      CHARACTER(LEN=*)   :: STRING,LABEL
      LOGICAL            :: JUMP
      INTEGER            :: IER

      CHARACTER(LEN=160) :: ST,TMPST
      CHARACTER(LEN=1)   :: NULL = CHAR(0)
      CHARACTER(LEN=1)   :: CTEMP
      INTEGER            :: RHEXP
      INTEGER            :: IEVAL,IND,MAXCH,ISEMICOL,IGOBRAK,NLP,ICHAR
      INTEGER            :: IFL,NXSTR,NNP,LHEXP,NRHP,IFUNC,IRTFLG
      INTEGER            :: IFLAG,N3,IGO,ITHEN,IEXIT,ICYCLE,NSTART,ILB
      INTEGER            :: ITAB,NE
      REAL               :: F1,F2

      CHARACTER(LEN=4)   :: COMPNEW = '<=>/'

      LOGICAL            :: BOOL(6,3) = RESHAPE(
     &     (/.FALSE.,.FALSE.,.TRUE., .TRUE., .FALSE.,.TRUE.,
     &       .TRUE., .TRUE., .TRUE., .FALSE.,.FALSE.,.FALSE.,
     &       .FALSE.,.TRUE., .FALSE.,.TRUE., .TRUE., .FALSE./),
     &     (/6,3/) )

      CHARACTER(LEN=2), PARAMETER ::COMP1(6) = 
     &                  (/'EQ','GE','LE','NE','GT','LT'/)
      CHARACTER(LEN=2), PARAMETER ::COMP2(6) = 
     &                  (/'==','>=','<=','/=','> ','< '/)

C     IEVAL WILL BE SWITCHED ON IF EXPRESSION FOUND IN 3RD POSITION
      IEVAL = 0  
      IND   = 0

C     SET NO ERROR FLAG
      IER   = 0

C     REMOVE BLANKS FROM INPUT STRING, PUT OUTPUT IN: ST
      CALL SHRINKQ(STRING,160,ST,MAXCH)

C     REMOVE ANY COMMENT FROM INPUT STRING
      ISEMICOL = SCAN(ST,';!')
      IF (ISEMICOL .GT. 0) MAXCH = ISEMICOL - 1

C     SEE IF '[' AND ']' NEED SYMBOL SUBSTITUTION, apr 2010 al
      IGOBRAK = INDEX(ST(1:MAXCH), '[') 
      IF (IGOBRAK .GT. 0) THEN
C         '[' AND ']' MAY NEED SYMBOL SUBSTITUTION E.G. [str]
          CALL SUBSYMPAR(ST(1:MAXCH),TMPST,MAXCH,0,IRTFLG)
          ST = TMPST
      ENDIF

C     LOOK FOR (
      NLP = INDEX(ST(:MAXCH),'(')
      IF (NLP .LE. 0) THEN
C        NO '(' FOUND, MAY BE PLAIN GOTO LB#
         ICHAR = INDEX(ST(1:MAXCH),'GOTOLB')
         IF (ICHAR <= 0) THEN
C           ERROR, NO 'GOTO' FOUND
            ICHAR = 1
            IFL   = 9
            GOTO 900   ! ERROR
         ENDIF
C        COPY LABEL STRING
         LABEL(1:5) = ST(ICHAR+4:ICHAR+7) // NULL
         JUMP       = .TRUE.
         RETURN
      ENDIF

C     ICHAR IS THE CURRENT POSITION IN ST
      ICHAR = NLP + 1

C     LOOK FOR FIRST PERIOD DELIMITER
      NXSTR = 2
      LHEXP = 0
      NNP   = 1

C     COMPUTE NNP: BALANCE OF ( AND ) { = NO. OF ('S MINUS NO OF )'S}
      DO 
         CTEMP = ST(ICHAR:ICHAR)
         IF (CTEMP == '.')                  EXIT ! START OF COMPARATOR
         IF (INDEX(COMPNEW(1:3),CTEMP) > 0) EXIT ! START OF COMP.: <=>
         IF ((ICHAR+1) <= MAXCH .AND.
     &        ST(ICHAR:ICHAR+1) .EQ. '/=')   EXIT ! START OF COMP.: /=

         IF (CTEMP == '(') NNP = NNP + 1
         IF (CTEMP == ')') NNP = NNP - 1
         ICHAR = ICHAR + 1
         IFL   = 2
         IF (ICHAR > MAXCH) GOTO 900   ! ERROR
         LHEXP = LHEXP + 1
      ENDDO

C     BEGINNING OF COMPARATOR FOUND (PART OF .LOGICAL. EXPRESSION)
26    IFL = 3
      IF (LHEXP .EQ. 0) GOTO 900   ! ERROR

C     EVALUATE LHEXP
      CALL EXPRQ(ST(NLP+1:),LHEXP,F1,IFLAG)
      IFL = 4
      IF (IFLAG .NE. 0) GOTO 900   ! ERROR

C     ICHAR NOW POINTS TO FIRST PERIOD OR NEW COMPARATOR
      DO IFUNC = 1,6
         IF (ST(ICHAR+1:ICHAR+2) == COMP1(IFUNC)(1:2)) THEN
            ICHAR = ICHAR + 4
            GOTO 60
         ENDIF
      ENDDO
      DO IFUNC = 1,4
         IF (ST(ICHAR:ICHAR+1)   == COMP2(IFUNC)(1:2)) THEN 
            ICHAR = ICHAR + 2
            GOTO 60
         ENDIF
      ENDDO
      DO IFUNC = 5,6
         IF (ST(ICHAR:ICHAR)     == COMP2(IFUNC)(1:1)) THEN
            ICHAR = ICHAR + 1
            GOTO 60
         ENDIF
      ENDDO
C     ERROR, LOGICAL COMPARATOR NOT IDENTIFABLE
      IFL = 5
      GOTO 900   ! ERROR


C     LOGICAL COMPARATOR IDENTIFIED
60    NRHP  = ICHAR

C     FIND RIGHT HAND EXPRESSION
      RHEXP = 0
      DO
         CTEMP = ST(ICHAR:ICHAR)
         IF (CTEMP == '(')  NNP = NNP + 1
         IF (CTEMP == ')')  NNP = NNP - 1
         IF (NNP == 0 .AND. CTEMP ==')') GOTO 64
         ICHAR = ICHAR + 1
         IFL   = 6
         IF (ICHAR > MAXCH) GOTO 900  ! ERROR
         RHEXP = RHEXP + 1
      ENDDO


64    IFL = 7
      IF (RHEXP == 0) GOTO 900
      N3 = ICHAR

C     EVALUATE RHEXP
      CALL EXPRQ(ST(NRHP:NRHP+RHEXP-1),RHEXP,F2,IFLAG)
      IFL = 8
      IF (IFLAG .NE. 0) GOTO 900   ! ERROR

C     ICHAR NOW POINTS TO THE )

      IGO   = INDEX(ST(ICHAR+1:MAXCH),'GO')
      IF (IGO .LE. 0) THEN
C        NO 'GOTO' ENCOUNTERED, DOES IT HAVE 'THEN'
         ITHEN  = INDEX(ST(ICHAR+1:MAXCH),'THEN')
         IEXIT  = INDEX(ST(ICHAR+1:MAXCH),'EXIT')
         ICYCLE = INDEX(ST(ICHAR+1:MAXCH),'CYCLE')
         IF (ITHEN > 0) THEN
C           HAS A 'THEN' INSTEAD OF A GOTO LABEL
            LABEL(1:5) = 'ELSE' // NULL

         ELSEIF (IEXIT > 0) THEN
C           HAS A 'EXIT' , MAY WANT TO FIND ENDDO
            LABEL(1:5) = 'ENDDO'

         ELSEIF (ICYCLE > 0) THEN
C           HAS A 'CYCLE', MAY WANT TO FIND ENDDO
            LABEL(1:5) = 'CYCLE'

         ELSE
C           NO 'THEN' SO SUBMIT TO EXPRESSION EVALUATOR
            LABEL(1:1) = ' '
            NSTART     = ICHAR + 1
            IEVAL      = 1
         ENDIF
      ELSE
C        SEARCH FOR "LB#"
         ILB = INDEX(ST,'GOTOLB')
         IF (ILB <= 0) THEN
C           ERROR, NO 'GOTOLB#' FOUND
            IFL = 9
            GOTO 900   ! ERROR
         ENDIF

C        ICHAR NOW POINTS TO THE LB# STRING
         ICHAR = ILB + 4

C        CHARACTERS "LB" FOUND.  COPY LABEL STRING
         LABEL(1:5) = ST(ICHAR:ICHAR+3) // NULL
      ENDIF

C     APPLY COMPARATION COMPUTATION.
      IF (F1 < F2) THEN
         ITAB = 1
      ELSEIF (F1 == F2) THEN
         ITAB = 2
      ELSE
         ITAB = 3
      ENDIF

C     SET JUMP
      JUMP = BOOL(IFUNC,ITAB)

C     LOGICAL VALUE NEGATED FOR IF...THEN JUMP
      IF (LABEL(1:4) == 'ELSE') JUMP = .NOT. JUMP

      IF ((LABEL(1:5) == 'ENDDO' .OR.
     &    LABEL(1:5) == 'CYCLE') .AND. (.NOT. JUMP)) THEN
C         'IF' FAILED SO JUST CONTINUE WITH OPERATION STREAM
          LABEL = ' '
          RETURN
      ENDIF

      IF (IEVAL .EQ. 0) RETURN
      IF (.NOT. JUMP)   RETURN

C     EVALUATE REGISTER EXPRESSION
      CALL ARASQ(ST(N3+1:N3+MAXCH-NRHP),MAXCH-NRHP,.FALSE.,IER)
      IF (IER .NE. 0) THEN
         WRITE(NOUT,904) ST(N3+1:N3+MAXCH-NRHP)
904      FORMAT(' *** ERROR EVALUATING: ',A)
         CALL ERRT(100,'LOGIFQ',NE)
      ENDIF
      RETURN



C     ERROR HANDLER
900   WRITE(NOUT,901) ST(ICHAR:MAXCH)
901   FORMAT(' *** IF STATEMENT SYNTAX ERROR STARTING AT: ',A)
      CALL ERRT(100,'LOGIFQ',NE)
      IER = 1

      END
