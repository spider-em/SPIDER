C++*************************************************** 12/11/79 1/11/81 VAX
C
C CALC.F                           REWRITTEN    MAY 98 ARDEAN LEITH                           
C                     BETTER ERROR MESSAGES   APR 2002 ARDEAN LEITH
C                     REG_GET_BYNUM           NOV 2005 ARDEAN LEITH
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
C  CALC(IRPNIN,NRPN,VAL,PIXVAL,RETVAL,IRTFLG)
C
C  PURPOSE:   EVALUATES EXPRESSIONS IN THE POSTFIX FORM, RETURNS VALUE
C
C  PARAMETERS:
C       IRPNIN      ARRAY CONTAINS POSTFIX NOTATION              (SENT)
C       NRPN        NO. OF ELEMENTS USED IN IRPNIN               (SENT)
C       VAL         ARRAY CONTAINS VALUES ASSOCIATED WITH        (SENT)
C                   INDICES IN IRPN. 
C       PIXVAL      CURRENT P1 PIXEL VALUE                       (SENT)  
C       RETVAL      CONTAINS VALUE OF EXPRESSION                 (RET.)
C       IRTFLG      ERROR FLAG                                   (RET.)
C
C  NOTES:
C       IRPN        DENOTES
C       > 200       PIXEL POINTER
C       128...200   REGISTER POINTER
C       97..114     UNITARY OPERATOR (SIN, COS, ETC)
C       41..47,96   BINARY OPERATOR (+-/^*)
C       48...88     VAL POINTER
C       0           SKIP
C       <0          BVAL POINTER
C  
C--*********************************************************************

	SUBROUTINE CALC(IRPNIN,NRPN,VAL,PIXVAL,RETVAL,IRTFLG)

        COMMON /UNITS/    LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

        PARAMETER         (IVALEN  = 40)
        PARAMETER         (IRPNLEN = 80)
        PARAMETER         (NFUNCT  = 16)

	DIMENSION         VAL(IVALEN)
C       BVAL IS USED TO AVOID OVERWRITING VAL EACH TIME CALC IS CALLED
        DIMENSION         BVAL(IVALEN)
        INTEGER           IRPNIN(IRPNLEN)
        DIMENSION         IRPN(IRPNLEN)
        LOGICAL           NOOP

C       SET DEFAULT ERROR RETURN FLAG
        IRTFLG = 1
        NOOP   = NRPN .EQ. 1

C       IRPN WILL BE DAMAGED SO THE INPUT MUST BE COPIED
	DO  N=1,NRPN
           IRPN(N) = IRPNIN(N)
        ENDDO

C       FOLLOWING LOOP IS EXECUTED FOR EACH OPERATION OR FUNCTION FOUND

2400    DO  N=1,NRPN
   
C          IRPNNOW IS CURRENT OPERATION OR VALUE POINTER           
           IRPNNOW = IRPN(N)

C          SKIP EVALUATING RPN POSITION IF IT IS ZERO
           IF (IRPNNOW .EQ. 0) GOTO 2001

           IF (IRPNNOW .GE. 201) THEN
C             PIXEL CONTENT POINTER, (IRPNNOW - 200) IS VAL POINTER
C             CONTENT OF VAL IS THE PIXEL NUMBER, (CURRENTLY 1 ONLY)
              LOC2       = LOC1
              LOC1       = IRPNNOW - 200
              BVAL(LOC1) = PIXVAL             
              IRPN(N)    = -LOC1
              NOOP       = .FALSE.
              GOTO 2400

           ELSEIF (IRPNNOW .GE. 128) THEN
C             REGISTER CONTENT POINTER, (IRPNNOW - 127) IS VAL POINTER
C             CONTENT OF VAL POSITION IS THE REGISTER NUMBER
              LOC2       = LOC1
              LOC1       = IRPNNOW - 127
              CALL REG_GET_BYNUM(INT(VAL(LOC1)),BVAL(LOC1),IRTFLG)
              NOOP       = .FALSE.
              IRPN(N)    = -LOC1
              GOTO 2400

	   ELSEIF ((IRPNNOW .GT. 40 .AND. IRPNNOW .LT. 48) .OR.
     &              IRPNNOW  .EQ. 94) THEN
C             IRPNNOW IS AN OPERATOR (+-/^*), THAT USES TWO OPERANDS
C             WHICH ARE KEPT IN VAL ARRAY.
C             LOC1 & LOC2 ARE ALREADY SET WHEN PREVIOUS VAL POINTERS 
C             WERE ENCOUNTERED.

              IF (IRPN(LOC2) .LT. 0) THEN
C                IRPNNOW POINTS TO BVAL
	         ITEMP1 = -IRPN(LOC2)
                 VALUE1 = BVAL(ITEMP1)
              ELSE
C                IRPNNOW POINTS TO VAL
	         ITEMP1 = IRPN(LOC2) - 48
                 VALUE1 = VAL(ITEMP1)
              ENDIF

              IF (IRPN(LOC1) .LT. 0) THEN
C                IRPNNOW POINTS TO BVAL
	         ITEMP2 = -IRPN(LOC1)
                 VALUE2 = BVAL(ITEMP2)
              ELSE
C                IRPNNOW POINTS TO VAL
	         ITEMP2 = IRPN(LOC1) - 48
                 VALUE2 = VAL(ITEMP2)
              ENDIF

              IF (IRPNNOW .EQ. 43) THEN
C                ADDITION, ICHAR('+') = 43
                 BVAL(ITEMP1) = VALUE1 + VALUE2

              ELSEIF (IRPNNOW .EQ. 45) THEN
C                SUBTRACTION, ICHAR('-') = 45
                 BVAL(ITEMP1) = VALUE1 - VALUE2

              ELSEIF (IRPNNOW .EQ. 42) THEN
C                MULTIPLICATION, ICHAR('*') = 42
                 BVAL(ITEMP1) = VALUE1 * VALUE2

              ELSEIF (IRPNNOW .EQ. 47) THEN
C                DIVISION, ICHAR('/') = 47
                 IF (VALUE2 .EQ. 0.) THEN
                    CALL ERRT(43,'CALC - DIVISION BY ZERO',NE)
                    RETURN
                 ENDIF
                 BVAL(ITEMP1) = VALUE1 / VALUE2

              ELSEIF (IRPNNOW .EQ. 94) THEN
C                POWER, ICHAR('^') = 94
                 BVAL(ITEMP1) = VALUE1 ** VALUE2

              ELSE
C                UNKNOWN OPERATOR
                 CALL ERRT(43,'CALC - UNKNOWN OPERATOR',NE)
                 RETURN
              ENDIF
              IRPN(LOC2) = -ITEMP1
              IRPN(LOC1) = 0
              IRPN(N)    = 0
C             START RPN EVALUATION LOOP ANEW
              GOTO 2400

           ELSEIF (IRPNNOW .LT. 0 .OR. 
     &            (IRPNNOW .GE. 48 .AND. IRPNNOW .LE. 88)) THEN
C             IRPNOW IS A NUMERICAL POINTER TO BVAL OR VAL
C             CONTENT OF BVAL OR VAL POSITION IS THE NUMBER

              LOC2 = LOC1
              LOC1 = N
              GOTO 2001

           ELSE
C             IRPNNOW IS A SINGLE OPERAND FUNCTION (LIKE: SIN, RAN, ETC)
C             THE FUNCTION USES ONE VALUE FROM VAL ARRAY.  LOC1 IS
C             ALREADY SET WHEN PREVIOUS VAL POINTER ENCOUNTERED

              IF (IRPN(LOC1) .LT. 0) THEN
C                IRPNNOW POINTS TO BVAL
	         ITEMP2 = -IRPN(LOC1) 
                 VALUE2 = BVAL(ITEMP2)
              ELSE
C                IRPNNOW POINTS TO VAL
	         ITEMP2 = IRPN(LOC1) - 48
                 VALUE2 = VAL(ITEMP2)
              ENDIF

C             CURRENTLY FUNCTIONS DENOTED BY 'a-q', SO IFUNC IS: 1..15
C             ICHAR('a') = 97
              IFUNC  = IRPNNOW - 96
              IF (IFUNC .LE. 0 .OR. IFUNC .GT. NFUNCT) THEN
C                UNKNOWN FUNCTION INDICATOR ENCOUNTERED IN IRPN ARRAY
                 CALL ERRT(43,'CALC - UNKNOWN OPERATOR',NE)
                 RETURN
              ENDIF

C             EVALUATE THE SINGLE OPERAND FUNCTION
              GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),IFUNC

C             PAD
1   	         NPOW  = INT(ALOG(VALUE2) / ALOG(2.0))
	         VALUE = 2**NPOW
                 IF (VALUE .GE. VALUE2)  THEN
	            BVAL(ITEMP2) = VALUE
                 ELSE
	            BVAL(ITEMP2) = 2 * VALUE
                 ENDIF
                 GOTO 1900

C             SINE
2                BVAL(ITEMP2) = SIN(VALUE2*0.017453293)
                 GOTO 1900

C             EXPONENTIATION
3	        BVAL(ITEMP2) = EXP(VALUE2)
                GOTO 1900
		
C             LOG10
4                IF (VALUE2 .LE. 0.) THEN
                    CALL ERRT(43,
     &                 'CAN NOT GET LOG OF NEGATIVE NUMBER)',NE)
                    RETURN
	         ENDIF
	         BVAL(ITEMP2) = ALOG10(VALUE2)
                 GOTO 1900

C             COSINE
5                BVAL(ITEMP2) = COS(VALUE2*0.017453293)
                 GOTO 1900

C             SQRT  	  
6                IF (VALUE2 .LT. 0.) THEN
                    CALL ERRT(43,
     &                 'CAN NOT GET SQRT OF NEGATIVE NUMBER)',NE)
                    RETURN
	         ENDIF
	         BVAL(ITEMP2) = SQRT(VALUE2)
                 GOTO 1900

C             NATURAL LOG
7                IF (VALUE2 .LE. 0.0) THEN
                    CALL ERRT(43,
     &                 'CAN NOT GET LOG OF NEGATIVE NUMBER)',NE)
                    RETURN
	         ENDIF
	         BVAL(ITEMP2) = ALOG(VALUE2)
                 GOTO 1900

C             INT
8	         BVAL(ITEMP2) = INT(VALUE2)
                 GOTO 1900

C             ABS
9	         BVAL(ITEMP2) = ABS(VALUE2)
                 GOTO 1900

C             ATAN
10	         BVAL(ITEMP2) = ATAN(VALUE2)*57.29578
                 GOTO 1900
	
C             ARC SIN
11	         IF (ABS(VALUE2) .GT. 1.0) THEN
	            CALL ERRT(43,
     &                 'CAN NOT GET ASIN OF NUMBER > 1.0)',NE)
                    RETURN
	         ENDIF
	         BVAL(ITEMP2) = ASIN(VALUE2) * 57.29578
                 GOTO 1900

C             ARC COS
12	         IF (ABS(VALUE2) .GT. 1.0) THEN
	            CALL ERRT(43,
     &                 'CAN NOT GET ACOS OF NUMBER > 1.0)',NE)
                    RETURN
	         ENDIF
	         BVAL(ITEMP2) = ACOS(VALUE2) * 57.29578
                 GOTO 1900

C             TANGENT
13	         BVAL(ITEMP2) = TAN(VALUE2 * 0.017453293)
                 GOTO 1900
 
C             RANDOM NUMBER UNIFORM [0,1]
14               CONTINUE
                 CALL  RANDOM_NUMBER(VALUE)
                 BVAL(ITEMP2) = VALUE
                 GOTO 1900

C             RANDOM NUMBER NORMAL(0,1)
15               CONTINUE
                 BVAL(ITEMP2) = RANN(0.0,1.0)
                 GOTO 1900

C             UNITARY SIGN CHANGE (NEGATION)
16	         BVAL(ITEMP2) = -VALUE2
C                GOTO 1900

1900          CONTINUE
C             CURRENT IRPN POSITION FINISHED
              IRPN(N)    = 0
C             IRPN SHOULD NOW POINT TO NEW BVAL LOCATION
              IRPN(LOC1) = -ITEMP2

C             START RPN EVALUATION LOOP ANEW
              GOTO 2400
           ENDIF

2001       CONTINUE
C          END OF EVALUATION LOOP -----------------------------------
        ENDDO


        IRTFLG = 0
        RETVAL = BVAL(1)
        IF (NOOP) RETVAL = VAL(1)

	END
