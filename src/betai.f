C++*********************************************************************
C
C BETAI.F 
C               
C **********************************************************************
C *                                                                         *
C **********************************************************************
C
C  BETAI(A,B,X) 
C
C  PURPOSE: RETURNS THE INCOMPLETE BETA FUNCTION  / (A,B)/X
C
C  AUTHOR: WILLIAM H PRESS ET AL., NUMERICAL RECIPES, PG 167
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        DOUBLE PRECISION FUNCTION BETAI(A,B,X)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        DOUBLE PRECISION ::  GAMMLN

        INCLUDE 'CMBLOCK.INC' 

        IF ( X < 0.0D0 .OR. X > 1.0D0) THEN
           WRITE(NOUT,*) '*** Bad argument X in BETAI:',X

        ELSEIF ( X == 0.0D0 .OR. X == 1.0D0) THEN
           BT = 0.0D0
        ELSE

C       FACTORS IN FRONT OF THE CONTINUED FRACTION.
        Y  = X
        BT = EXP( GAMMLN(A + B) - GAMMLN(A) - GAMMLN(B) +
     &               A * DLOG(Y) + B * DLOG(1. - Y))

        ENDIF

        IF ( X < ((A + 1.) / (A + B + 2.0))) THEN

C          USE CONTINUED FRACTION DIRECTLY.
           BETAI = BT * BETACF (A, B, X) / A
        ELSE

C          USE CONTINUED FRACTION AFTER MAKING THE SYMMETRY TRANSFORMATION
           BETAI = 1.0 - (BT * BETACF(B, A, 1.0 - X) / B)
        ENDIF

        END

C++*********************************************************************
C
C GAMMLN.F
C
C **********************************************************************
C
C PURPOSE: RETURNS THE VALUE ln(XX) FOR XX > 0. FULL ACCURACY IS 
C          OBTAINED FOR  0 < XX < 1, THE REFLECTION FORMULA CAN BE USED
C          FIRST.
C
C--*********************************************************************

        DOUBLE PRECISION FUNCTION GAMMLN(XX)

        DOUBLE PRECISION :: COF(6), STP, HALF, ONE, FPF, TMP, SER
        DOUBLE PRECISION :: X, XX
        INTEGER          :: J


        DATA COF, STP/76.18009173D0, -86.50532033D0, 24.01409822D0,
     &          -1.231739516D0, .120858003D-2, -.536382D-5, 
     &          2.50662827465D0/

        DATA HALF, ONE, FPF/0.5D0, 1.0D0, 5.5D0/

        X   = XX - ONE
        TMP = X + FPF
        TMP = (X + HALF) * DLOG(TMP) - TMP
        SER = ONE

        DO J = 1, 6
          X   = X + ONE
          SER = SER + COF(J) / X
        END DO

        GAMMLN = TMP + DLOG(STP * SER)

        END 


C++*********************************************************************
C
C **********************************************************************
C BETACF(A, B, X)
C
C	CONTINUED FRACTION EVALUATION BY THE RECURENCE METHOD
C       EQUATION 5.2.5  IN BOOK
C         
C--*********************************************************************

        DOUBLE PRECISION FUNCTION BETACF(A, B, X)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        INCLUDE 'CMBLOCK.INC'       

        PARAMETER (ITMAX = 100, EPS = 3. E-7)

        AM = 1.0
        BM = 1.0
        AZ = 1.0

        QAB = A + B
        QAP = A + 1
        QAM = A - 1
        BZ  = 1.0 - (QAB * X / QAP)

C	CONTINUED FRACTION EVALUATION BY THE RECURENCE METHOD
C       EQUATION 5.2.5  IN BOOK

        DO M = 1, ITMAX
          EM = M
          TEM = EM + EM
          D = EM * (B - M) * X / ((QAM + TEM) * (A + TEM))

C         ONE STEP (THE EVEN ONE) OF THE RECURENCE
          AP = AZ + D * AM
          BP = BZ + D * BM
          D  = - (A + EM) * (QAB + EM) * X / ((A + TEM) * (QAP + TEM))

C         NEXT STEP OF THE RECURRENCE (THE ODD ONE)
          APP = AP + D * AZ
          BPP = BP + D * BZ

C         SAVE THE OLD ANSWER
          AOLD = AZ

C         RENORMALIZE TO PREVENT OVERFLOW
          AM = AP / BPP
          BM = BP / BPP
          AZ = APP / BPP
          BZ = 1.0
   
C         ARE WE DONE ?
          IF(DABS(AZ - AOLD) .LT. EPS * DABS(AZ)) GOTO 1
        END DO

        WRITE(NOUT,*) 
     &        '***  IN BETACF, A or B too big, or ITMAX too small '

1       BETACF = AZ

        END
