
C++*********************************************************************
C
C IFTORPN.F    ADAPTED FROM LOGIF.F FOR CHAR. VARIABLES   AUG 1989  al
C              UNCONDITIONAL JUMP ADDED                  SEPT 1996  al
C              IF (...) THEN IMPLEMENTED                 SEPT 1997  al
C              CHARINSIDE PARAMETERS CHANGED             JAN  2001  al
C              POLISH PARAMETERS                         DEC  2005  al
C              ! COMMENT DELIMITER                       DEC  2011  aL
C
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
C     SUBROUTINE IFTORPN(STRING,ICOMPREP,
C             IRPN1,NRPN1,VAL1, IRPN2,NRPN2,VAL2, IRPN3,NRPN3,VAL3,
C             IRTFLG)
C
C    PURPOSE:     EVALUATES STRING SUCH AS: IF(X.LE.5) P2 = F(P1)
C                 RETURNS THE LB77, RETURNS RPN AND VAL ARRAY 
C                 FOR ALL 3 EXPRESSIONS.
C
C    PARAMETERS:  STRING       INPUT LINE                     (SENT)
C                 ICOMPREP     COMPARISION INDICATOR (.GT.)   (RETURNED)
C                 IRPN?        RPN STRING                     (RETURNED)
C                 NRPN?        NO. ELEMENTS IN RPN STRING     (RETURNED)
C                 VAL?         VALUES POINTED TO BY RPN       (RETURNED)
C                 IRTFLG       ERROR FLAG (0 IS NORMAL)       (RETURNED)
C
C--*********************************************************************

      SUBROUTINE IFTORPN(STRING,ICOMPREP,
     &      IRPN1,NRPN1,VAL1, IRPN2,NRPN2,VAL2, IRPN3,NRPN3,VAL3,
     &      IRTFLG)

      INCLUDE 'CMBLOCK.INC'


      PARAMETER      (IVALEN = 40)
      PARAMETER      (IRPNLEN = 80)

      INTEGER        RHEXP
      CHARACTER *(*) STRING
      CHARACTER *80  ST
      CHARACTER *2   COMP(6),COMPRET
      DIMENSION      VAL1(IVALEN),   VAL2(IVALEN),   VAL3(IVALEN)
      INTEGER        IRPN1(IRPNLEN), IRPN2(IRPNLEN), IRPN3(IRPNLEN)
      LOGICAL        ISCHAR
      EXTERNAL       ISCHAR

      DATA COMP/'EQ', 'GE', 'LE', 'NE', 'GT', 'LT'/

C     SET ERROR FLAG
      IRTFLG   = 1

C     REMOVE BLANKS FROM INPUT STRING, PUT OUTPUT IN ST
      CALL SHRINKQ(STRING,80,ST,MAXCH)

C     REMOVE ANY COMMENT FROM INPUT STRING
      ISEMICOL = SCAN(ST,';!')
      IF (ISEMICOL .GT. 0) MAXCH = ISEMICOL - 1

C     LOOK FOR FIRST '('
      NLP = INDEX(ST(:MAXCH),'(')
      IF (NLP .LE. 0) GOTO 900

C     ICHAR IS THE CURRENT POSITION IN ST

C     FIND LEFT-HAND-EXPRESSION (FOLLOWED BY .LOGICAL. COMPARATOR)
      ICOMPREP = 0 
      DO ICHAR = NLP + 1, MAXCH - 3
C        CONTINUE STEPPING THRU EXPRESSION TILL .LOGICAL. FOUND
         IF (ST(ICHAR:ICHAR) .EQ. '.' .AND. 
     &       ISCHAR(ST(ICHAR+1:ICHAR+1)) .AND.
     &       ISCHAR(ST(ICHAR+2:ICHAR+2))) THEN

C           PERIOD FOUND, FOLLOWED BY .LOGICAL. COMPARATOR 
            ICOMPREP = 0
            DO IFUNC = 1,6
               IF (ST(ICHAR+1:ICHAR+2) .EQ. COMP(IFUNC)(1:2)) THEN
                  ICOMPREP = IFUNC
                  GOTO 60
               ENDIF
            ENDDO
         ENDIF
      ENDDO

C     ERROR IF LOGICAL COMPARATOR NOT FOUND AND IDENTIFED
60    IF (ICOMPREP .EQ. 0) GOTO 900

C     LOGICAL COMPARATOR IDENTIFIED, CONVERT LHEXP TO RPN
      CALL POLISH(0,ST(NLP+1:ICHAR-1),ICHAR-NLP-1,
     &            IRPN1,NRPN1,VAL1,NVAL1,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 900

C     FIND RIGHT HAND EXPRESSION (RHE)
      IGO  = ICHAR + 4
      IEND = INDEX(ST,')THEN') - 1

C     LOGICAL COMPARATOR IDENTIFIED, CONVERT RHE TO RPN
      CALL POLISH(0,ST(IGO:IEND),IEND-IGO+1,
     &            IRPN2,NRPN2,VAL2,NVAL2,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 900

C     FIND ASSIGNMENT WHICH FOLLOWS RIGHT-HAND-EXPRESSION 
      CALL CHARINSIDE(ST,')THEN','=',.TRUE.,.FALSE.,IGO,IEND,NCHAR)
      IF (NCHAR .LE. 0) GOTO 900

C     CONVERT OPERATIONAL EXPRESSION FOLLOWING ASSIGNMENT TO RPN
      CALL POLISH(0,ST(IEND+2:MAXCH),MAXCH-IEND-1,
     &            IRPN3,NRPN3,VAL3,NVAL3,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 900

      IRTFLG = 0
      RETURN

C     ERROR HANDLER --------------------------------------------
900   WRITE(NOUT,901) ST(ICHAR:MAXCH)
901   FORMAT(' *** IF STATEMENT SYNTAX ERROR STARTING AT: ',A)
      CALL ERRT(100,'LOGIFQ',NE)
      IRTFLG = 1
      RETURN

      END
