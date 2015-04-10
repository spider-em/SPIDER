
C++*********************************************************************
C
C EDGE.FOR                                     LONG FILENAMES JAN 89 al
C                         LONGER BUG               APR 02 ArDean Leith
C                        OPFILEC                  FEB  03 ARDEAN LEITH
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
C  EDGE(LU1,LU3,LU2,NSAM,NROW)
C
C  PURPOSE:    EDGE DETECTION
C      
C  PARAMETERS:
C        LU1      LOGICAL UNIT NUMBER OF FILE
C        LU2      LOGICAL UNIT NUMBER OF FILE
C        NSAM     NUMBER OF SAMPLES
C        NROW     NUMBER OF ROWS
C       
C
C IMPLEMENTATION OF EDGE DETECTION BY WIENER FILTERING
C ORIGINAL VERSION : R.W. FRIES, RPI, TROY, NY
C MODIFIED VERSION : R. MARSHALL, 
C NY STATE DEPT. OF HEALTH, ALBANY, NY  JULY 1977
C FOR DETAILED INFORMATION CF. WRITEUP "EDGE DETECTION IN NOISY
C IMAGES USING RECURSIVE DIGITAL FILTERING" BY J.W. MODESTINO
C AND R.W. FRIES, ELECTRICAL AND SYSTEMS ENGINEERING DEPT.,
C RENSSELAER POLYTECHNIC INSTITUTE, TROY, NY.12181
C
C THE FOLLOWING FILTER TYPES ARE BEING USED:
C
C ########################## 1 ######## 2 ####### 3 ###### 4 ######
C
C XI (S/N RATIO )    =   INFINITE       10       10       3 DB
C RHO(CORRELATION)   =      -          -0.9      0.0      0.5
C LAMBDA(EDGE FREQU.)=      0.05        0.05     0.025    0.0125
C  (PER UNIT LENGTH)
C
C--*******************************************************************

      SUBROUTINE EDGE(LU1,LU3, LU2, NSAM, NROW)

      INCLUDE 'CMBLOCK.INC'

      DIMENSION X(10*NSAM)

      CHARACTER *11 SCR
      INTEGER       CUR,PREV,PREVY,CURY,NEXTY,REC
      REAL          FB11(4),FA10(4),FA11(4),FA(4),FF(3)

      DATA SCR(1:7)/'SCR999.'/

C     FILTER TYPES ####### 1 ####### 2 ###### 3 ###### 4 ######
      DATA FB11 /-0.8256, -0.9244, -0.8460, -0.7556/
      DATA FA10 /-0.2149, -0.1944, -0.2849, -0.4815/
      DATA FA11 / 0.0098, -0.1929, -0.1687,  0.0102/
      DATA FA   / 0.115,   0.084,   0.053,   0.0088/

      DATA FF / 0.2, 0.4, 0.6/

C     DEFAULT LOGICAL UNIT ASSIGNMENTS
      LUN1 = LU1
      LUN2 = LU2
      LUN3 = LU3

C     CREATE SCRATCH FILE NAME ("." IS IN SCRATCH)
      NLET = 7
      CALL LONGER(SCR,NLET,DATEXC,IRTFLG)

      CUR  = 1
      PREV = 2

C     ADDRESS ASSIGNMENTS
      IX = 0
      IXOUT = IX + 2*NSAM
      IY    = IXOUT + NSAM
      IYL   = IY + 3*NSAM
      IYR   = IYL + 2*NSAM
      ITOT  = IYR + 2*NSAM
      
      MAXIM = 0
      IFORM = 1
      CALL OPFILEC(0,.FALSE.,SCR,LUN2,'U',IFORM,NSAM,NROW,1,
     &             MAXIM,' ',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0)  RETURN

843   CALL RDPRMI(IFILTR, IDUM, NOT_USED, 'FILTER NUMBER (1-4)')

      IF (IFILTR .LE. 0 .OR. IFILTR .GE. 5) THEN
        CALL ERRT(16,'EDGE',NE)
        GOTO 843
      ENDIF

      B11 = FB11(IFILTR)
      A10 = FA10(IFILTR)
      A11 = FA11(IFILTR)
      A   = FA(IFILTR)

      IFTH = 2
      CALL RDPRI1S(IFTH,NOT_USED,
     &   'THRESHOLD (1)LOW, (2)MEDIUM, (3)HIGH',IRTFLG)
      FAC    = FF(IFTH)
      THRESH = 0.
      B10    = -.5*(B11+1.)

C     INITIALIZE PREVIOUS FILTER OUTPUTS TO ZERO

      NSAM0 = 10*NSAM
      DO  I = 1,NSAM0
           X(I) = 0.0
      ENDDO

C     INITIALIZE CURRENT INPUT LINE
      CALL REDLIN (LUN1,X((CUR-1)*NSAM+1),NSAM,NROW)

C     MOVE DOWN THROUGH THE IMAGE LINE BY LINE

      DO I=NROW-1,2,-1

C       MAKE OLD CURRENT LINE NEW PREVIOUS LINE

        IHOLD = CUR
        CUR   = PREV
        PREV  = IHOLD
        NCUR  = (CUR-1)*NSAM
        NPREV = (PREV-1)*NSAM

C       READ IN NEW CURRENT LINE

        CALL REDLIN(LUN1, X(NCUR+1), NSAM, I)
 
C       MOVE THROUGH LINE POINT BY POINT

        DO  J=2,NSAM-1

C          RIGHT MOVING FILTER SECTION

           JCUR = NCUR+J
           JPREV = NPREV+J
           JPREV1 = JPREV-1
           JCUR1 = JCUR-1
      X(IYR+JCUR) = -A10*(X(IYR+JCUR1)+X(IYR+JPREV)) + B10*(X(JPREV)+
     1     X(JCUR1)) + B11*X(JPREV1) - A11*X(IYR+JPREV1) + X(JCUR)

C          LEFT MOVING FILTER SECTION

           JL = NSAM - J + 1
           JLCUR = NCUR + JL
           JLCUR1 = JLCUR + 1
           JLPREV = NPREV + JL
           JLPRE1 = JLPREV + 1

           X(IYL+JLCUR) = -A10*(X(IYL+JLCUR1)+X(IYL+JLPREV)) + X(JLCUR)
     1       + B10*(X(JLCUR1)+X(JLPREV)) + B11*X(JLPRE1)
     2        - A11*X(IYL+JLPRE1)
C          END POINT BY POINT LOOP
        ENDDO


C       SUM LEFT AND RIGHT FILTER OUTPUTS

        DO  J=2,NSAM-1
           X(IY+J) = X(IYR+NCUR+J) + X(IYL+NCUR+J+1)
        ENDDO

C       WRITE OUT OUTPUT OF DOWN MOVING FILTER SECTIONS

        CALL WRTLIN(LUN2, X(IY+1), NSAM, I)

C       END LINE BY LINE LOOP MOVING DOWN

        ENDDO


C       INITIALIZE PREVIOUS FILTER OUTPUTS
C       YL CHANGED TO IYL JUNE 93 al to avoid alpha compiler error message
C       (LOOKS LINE THIS WAS NEVER CORRECT!)

        DO J=IYL+1,NSAM0
           X(J) = 0
        ENDDO


1125    DO L=1,NSAM
          X(IY+NSAM+L) = 0
          X(IY+2*NSAM+L) = 0
        ENDDO

        X(IXOUT+1) = 0
        X(IXOUT+NSAM-1) = 0
        X(IXOUT+NSAM) = 0

C       INITIALIZE POINTERS FOR FINAL EDGE OUTPUT

        PREVY = 1
        CURY = 2
        NEXTY = 3

C       REREAD FIRST INPUT LINE
C       ISET INITIALZED TO 0 TO AVOID ALPHA COMPILER ERROR JUNE 93 al
        ISET = 0
        IF(ISET.NE.3) CALL REDLIN(LUN1, X((CUR-1)*NSAM+1), NSAM, 1)

C       REC+1 IS FILTER OUTPUT POINTER
C       REC IS EDGE OUTPUT POINTER
C       MOVE THROUGH IMAGE LINE BY LINE UP

        DO  REC=2,NROW-2

C       UPDATE POINTERS FOR CURRENT, PREVIOUS AND NEXT LINES
C       OF FILTER OUTPUT

        IHOLD = PREVY
        PREVY = CURY
        CURY = NEXTY
        NEXTY = IHOLD
        
        NCURY = (CURY-1)*NSAM
        NPREVY = (PREVY-1)*NSAM
        NNEXTY = (NEXTY-1)*NSAM
        IF(ISET.EQ.3) GO TO 1550

C       UPDATE POINTERS FOR CURRENT AND PREVIOUS INPUT LINES

        IHOLD = CUR
        CUR = PREV
        PREV = IHOLD
       NCUR = (CUR-1)*NSAM
      NPREV = (PREV-1)*NSAM

C       READ IN CURRENT LINE

        CALL REDLIN(LUN1, X(NCUR+1), NSAM, REC)

C       MOVE THROUGH LINE POINT BY POINT

        DO  J=2,NSAM-1
      JCUR = NCUR + J
      JCUR1 = JCUR - 1
      JPREV = NPREV + J
      JPREV1 = JPREV - 1

C       RIGHT MOVING FILTER SECTION

        X(IYR+JCUR) =-A10*(X(IYR+JCUR1)+X(IYR+JPREV))-A11*X(IYR+JPREV1)
     1  + X(JCUR) + B10*(X(JCUR1)+X(JPREV)) + B11*X(JPREV1)

C       LEFT MOVING FILTER SECTION

        JL = NSAM-J+1

      JLCUR = NCUR + JL
      JLCUR1 = JLCUR + 1
      JLPREV = NPREV + JL
      JLPRE1 = JLPREV + 1
        X(IYL+JLCUR) =-A10*(X(IYL+JLCUR1)+X(IYL+JLPREV)) 
     1 -A11*X(IYL+JLPRE1)
     2 + X(JLCUR) +B10*(X(JLCUR1)+X(JLPREV)) + B11*X(JLPRE1)

C       END OF POINT BY POINT LOOP

        ENDDO

C       READ IN OUTPUT OF DOWN MOVING FILTER SECTION

        CALL REDLIN(LUN2, X(IY+NNEXTY+1), NSAM, REC+1)

C       ADD UP AND DOWN MOVING FILTER SECTIONS

        DO J=2,NSAM-1
           X(IY+NNEXTY+J) = A*(X(IYR+NCUR+J) + X(IYL+NCUR+J+1)
     1      + X(IY+NNEXTY+J))
        ENDDO

        CALL WRTLIN(LUN2, X(IY+NNEXTY+1), NSAM, REC)
        
        IF(ISET .EQ.2) GO TO 1600

C       MOVE THROUGH POINT BY POINT AND DETERMINE IF 
C       EDGE ELEMENT IS PRESENT

1550    IF(ISET.EQ.3) CALL REDLIN(LUN2,X(IY+NNEXTY+1),NSAM,REC+1)

        DO I=2,NSAM-2
      INEXT = NNEXTY + I
      INEXTN = INEXT - 1
      INEXTP = INEXT + 1
      IPREVY = NPREVY + I
      IPREVP = IPREVY + 1
      IPREVN = IPREVY - 1
      ICURP  = NCURY + I + 1
      ICURN  = ICURP - 2
      DIF1   = ABS(X(IY+INEXTN)-X(IY+IPREVP))
      DIF2   = ABS(X(IY+INEXT)-X(IY+IPREVY))*1.414
      DIF3   = ABS(X(IY+INEXTP)-X(IY+IPREVP))
      DIF4   = ABS(X(IY+ICURN)-X(IY+ICURP))*1.414

        DIF = 0
        DIF = AMAX1(DIF,DIF1,DIF2,DIF3,DIF4)

        ZSUM = X(IY+ICURN)    + X(IY+NCURY+I) + X(IY+ICURP)
     1         + X(IY+INEXTN)
     2         + X(IY+INEXT)  + X(IY+INEXTP) + X(IY+IPREVN) 
     3         + X(IY+IPREVY) + X(IY+IPREVP)

           ZSUM =ZSUM*ZSUM*FAC

C          INITIALIZE TO NO EDGE

           X(IXOUT+I) = 0

C          IF AN EDGE SET TO MAXIMUM VALUE

           IF(DIF*DIF-ZSUM .GT. THRESH) X(IXOUT+I) = 2.
        ENDDO

C       END POINT BY POINT LOOP


C       WRITE OUT EDGE INFORMATION
        CALL WRTLIN(LUN3, X(IXOUT+1), NSAM, REC)

1600    CONTINUE
        ENDDO

        DO I=2,NSAM-2
          X(IXOUT+I) = 0
        ENDDO

        IF (ISET .EQ.2) GO TO 2150
        CALL WRTLIN(LUN3, X(IXOUT+1), NSAM, 1)
        CALL WRTLIN(LUN3, X(IXOUT+1), NSAM, NROW)
        CALL WRTLIN(LUN3, X(IXOUT+1), NSAM, NROW-1)
        
        IF (ISET .EQ.3) GO TO 2200


2150    CALL WRTLIN(LUN2, X(IXOUT+1), NSAM, 1)
        CALL WRTLIN(LUN2, X(IXOUT+1), NSAM, NROW-1)
        CALL WRTLIN(LUN2, X(IXOUT+1), NSAM, NROW)

2200    CLOSE(LUN2)

C       END LINE BY LINE LOOP

        RETURN
        END




