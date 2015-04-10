
C++*********************************************************************
C
C GRAPHS.F
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
C       GRAPHS(NDEV,YVAL,NPT,NSET,IDEV,FACT,IRTFLG)
C
C       DISPLAYS NSET ONE-DIMENSIONAL FUNCTIONS STORED
C       IN AN ARRAY YVAL(NPT) ON THE LINE PRINTER  (ANCIENT CODE!!)
C
C       NDEV    OUTPUT LUN
C       YVAL    FLOATING POINT ARRAY OF DIMENSION NSET*NPT CONTAINING
C               NSET FUNCTIONS
C       NPT     DIMENSION OF EACH FUNCTION ( MUST ALL HAVE SAME DIM.)
C       NSET    NUMBER OF FUNCTIONS TO BE DISPLAYED ON A SINGLE
C               GRAPH < OR = 10
C       IDEV    0 IS TERMINAL OUTPUT 1 IS LINE PRINTER OUTPUT
C       FACT    Y SCALE MAGNIFCATION FACTOR
C       IRTFLG  ERROR FLAG (0 IS NORMAL)
C
C  DESCRIPTION:  THE RANGE BETWEEN MAXIMUM AND MINIMUM OF THE WHOLE
C                ARRAY YVAL IS MAPPED ONTO ONE PRINTER PAGE, WHILE
C                THE NPT SAMPLING POINTS ARE SPREAD ACROSS THE FULL
C                WIDTH OF THE PAGE.  EACH OF THE NSET FUNCTIONS IS
C                POINTED WITH A DIFFERENT SYMBOL, STARTING
C                WITH '0', '1', ETC.
C
C--*******************************************************************

      SUBROUTINE GRAPHS(NDEV,YVAL,NPT,NSET,IDEV,FACT,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      CHARACTER*1 IPLOT(128),ICHAR(10),ICH
      DIMENSION   YVAL(1),ISTORE(5)

      DATA ICHAR/'0','1','2','3','4','5','6','7','8','9'/

C     (YSIZE + 1) IS Y DIMENSION OF GRAPH  
C     DEFAULT IS TERMINAL OUTPUT  70 X 18
      IXSIZE    = 70
      YSIZE     = 18.0

      IF (IDEV .NE. 0) THEN
C        OUTPUT TO LINE PRINTER  128 X 50
         IXSIZE  = 128
         YSIZE   = 50.0
      ENDIF

      IYSIZE = YSIZE + 0.5
      KK     = NSET * NPT

C     FIND MIN/MAX OF YVAL(1...NSET*NPT)
      TOP   = YVAL(1)
      BOT   = YVAL(1)
      DO   I = 1,KK
         IF (YVAL(I) .GT. TOP) TOP = YVAL(I)
         IF (YVAL(I) .LT. BOT) BOT = YVAL(I)
      ENDDO

      IF (TOP .LE. BOT) THEN
         WRITE(NDEV,*) '*** ZERO Y RANGE. GRAPH ABANDONED'
         IRTFLG = 1
         RETURN
      ENDIF

      T1 = BOT + (TOP - BOT) / FACT
      WRITE(NDEV,*) 
     &' --------------------------- GRAPH -------------------------'

      WRITE(NDEV,13) BOT,TOP, BOT,T1
13    FORMAT(/'  Y RANGE ',G12.3,' TO ',G12.3,
     &        ' (ACTUAL ',G12.4,' TO ',G12.3,')')
      WRITE(NDEV,*) ' '

C     SCALE Y VALUES 
      S = FACT * YSIZE / (TOP - BOT)
      DO  I = 1,KK
         YVAL(I) = S * (YVAL(I) - BOT)
      END DO

C     UNSURE WHY THIS WAS IN HERE al
      T1  = -S * BOT + 0.5
      IF (T1 .LT. 0.0) T1 = T1 - 1.0
      TIX = 32000.0
      IF (T1 .LT. TIX) TIX = T1
      IX  = TIX

      DX =  FLOAT(IXSIZE-1) / FLOAT(NPT-1)

      DO IY = IYSIZE, 0, -1
C        PUT BACKGROUND CHAR(S). IN THIS LINE OF GRAPH
         ICH = ' ' 
         IF (IY.EQ.IYSIZE .OR. IY.EQ.0 .OR. IY.EQ.IX) ICH = '-'
         DO  I = 1,IXSIZE
            IPLOT(I) = ICH
         ENDDO
         IPLOT(1)      = 'I'
         IPLOT(IXSIZE) = 'I'

C        CHECK IF YVAL IS AT THIS LINE OF GRAPH
         X = 1.0
         DO  I = 1,NPT
            K = I
            DO  J = 1,NSET
              IF (IFIX(YVAL(K)+.5).EQ.IY) IPLOT(IFIX(X+.5)) = ICHAR(J)
              K = K + NPT
            ENDDO
            X = X + DX
         ENDDO

C        OUTPUT THIS LINE OF THE GRAPH
         WRITE(NDEV,8) (IPLOT(I),I=1,IXSIZE)
8        FORMAT('   ',128A)
      ENDDO

C     PLOT X AXIS LABELS
      CALL GRAPHAX(NDEV,IPLOT,ICHAR,IXSIZE,IRTFLG)

      WRITE(NDEV,*) ' '

C     RESTORE ORIGINAL YVAL. VALUES
      S = 1.0 / S
      DO I=1,KK
         YVAL(I) = YVAL(I) * S + BOT
      END DO

      IRTFLG = 0

      END

C     ---------------- GRAPHAX -------------------------------------

      SUBROUTINE GRAPHAX(NDEV,IPLOT,ICHAR,IXSIZE,IRTFLG)



      INCLUDE 'CMBLOCK.INC'

      CHARACTER*1 IPLOT(128),ICHAR(10)
      DIMENSION   ISTORE(5)

C     LABEL  X AXIS, WITH SCALES 

	XO    = -2.0
	FMAXQ = FMAX / 100.0
	FMINQ = FMIN / 100.0
	DO  IXO=1,4
          IF (FMAXQ - FMINQ .GE. 2.) GO TO 200
          XO    = XO + 1.0
          FMAXQ = FMAXQ * 10.0
          FMINQ = FMINQ * 10.0
	ENDDO

200	FJ    = (IXSIZE - 1.0) / (FMAXQ - FMINQ)
   	IMINH = FMINQ
	IF (FLOAT(IMINH) .LT. FMINQ) IMINH = FMINQ + 1.0
	FMINV = FLOAT(IMINH)
	IMAXH = FMAXQ
	FMAXV = FLOAT(IMAXH)
   	STEP  = (FMAXV-FMINV)/4.
	DO  I=1,IXSIZE
  	   IPLOT(I) = ' '
	ENDDO

	DO  I=0,4
          BT = FMINV + FLOAT(I) * STEP
          J  = INT((BT - FMINQ) * FJ) + 1.5
          IF (J .GT. IXSIZE) J = IXSIZE
          IF (J .LT. 1) J = 1
          IPLOT(J)    = '.'
          ISTORE(I+1) = J
	ENDDO

	WRITE(NDEV,260)(IPLOT(I),I=1,IXSIZE)
260	FORMAT('   ',128A1)

	FMAXV = FMAXV / (10.0 ** XO)
	FMINV = FMINV / (10.0 ** XO) + 0.00001
	IF (FMINV .LT. 0.0) FMINV = FMINV  - 0.00002
	STEP = STEP / (10.0 ** XO)
	DO  IA=1,5
          FVALUE = FMINV + FLOAT(IA-1) * STEP
          FVALUE = FVALUE / 100.0
          DO  IB=1,9
            IC = IB-5
            IF (ISTORE(IA) .LE. 4)   IC = IB - ISTORE(IA)
            IF (ISTORE(IA) .GT. (IXSIZE -4)) THEN 
               IC = IC - ISTORE(IA) + (IXSIZE-4)
            ENDIF

            IF (IB .EQ .5) THEN
               IPLOT(ISTORE(IA)+IC) = '.'

            ELSEIF (IB .EQ. 1) THEN
               IPLOT(ISTORE(IA)+IC) = '+'
               IF (FVALUE .LT. 0) IPLOT(ISTORE(IA)+IC) = '-'
            ELSE
               IVALUE = FVALUE
               DO ID=1,10
                 IF (IABS(IVALUE) .EQ. ID-1) THEN
                    IPLOT(ISTORE(IA)+IC) = ICHAR(ID)
                 ENDIF
               ENDDO
               FVALUE = FVALUE - FLOAT(IVALUE)
               FVALUE = FVALUE * 10.0
            ENDIF
	 ENDDO
      ENDDO
      WRITE(NDEV,260) (IPLOT(I),I=1,IXSIZE)

      IRTFLG = 0

      END
