C++*********************************************************************
C
C DEFO1001.F
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
C   DEFO1001(NUM,NP,KP,NA,NX,SPMAX)
C
C   USING LEAST SQUARE METHOD TO DETERMINE DEFOCUS, AMPLITUDE CONTRAST
C   (WEIGHTED FIT; CONSTRAIN REQUIRE A2 IS THE SAME FOR DEFOCUS SERIES)
C   X(K,A) = PI*(0.5*CS*LAMBDA**3*K**4-DZ*LAMBDA*K**2)-OFFSET
C   X(K,A) = PI*(0.5*CS*LAMBDA**3*K**4-A1*LAMBDA*K**2)-A2
C
C   PARAMETERS:
C	NUM      NUMBER OF IMAGES
C       NP(I)    NUMBER OF MINIMUM IN EACH IMAGE
C	KP(I,J)  ARRAY OF SP. FREQ. POINTS OF MINIMUM
C	NA(I,J)  ARRAY OF ABBERATION
C	NX       IMAGE X DIMENSION
C	SPMAX    MAX OF SP. FREQ.
C	NA       NUMBER OF ABBERATION IN UNIT OF PI
C       
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

	SUBROUTINE DEFO1001(NUM,NP,KP,NA,NX,SPMAX)

        INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 
	
	DIMENSION NP(*)
	REAL      KP(20,20),NA(20,20)

        REAL      A1(20)
	EQUIVALENCE (A0(1),A1(1))
        COMMON  A0(20),B(20),C(20,20),DA1(10),DXA1(10),NM(30),
     &          XA1(30,30),XA2(30),XM(30),Y1(20,10),Y2(20,10),
     &          WEIGHT(20,10),F2(512)

	REAL                   :: KM,KS,LAMBDA,KF
	CHARACTER (LEN=1)      :: CHO1
	CHARACTER (LEN=1)      :: NULL = CHAR(0)
        CHARACTER (LEN=MAXNAM) :: OUTNAME

	REAL, PARAMETER        :: PI = 3.141592654
	INTEGER, PARAMETER     :: LUN2 = 20

	AA0 = -100.0  
 
	DO I=1,30
           NM(I) = 0
	ENDDO

C       INPUT EM PARAMETERS

	CALL RDPRM1S(LAMBDA,NOT_USED,'WAVELENGTH LAMBDA [A]',IRTFLG)

	CALL RDPRM1S(CS,NOT_USED,'SPHERICAL ABERRATION CS [MM]',IRTFLG)

        call flushresults

        IF (CS < 0.0001) CS = 0.0001

	CS = CS * 1.0E07
	KM = SPMAX
	KS = KM / FLOAT(NX)

C       GET VALUE OF Y1  
	DO  I=1,NUM
	   DO  J=1,NP(I)
	      Y1(I,J) = PI * NA(I,J)
	   ENDDO
	ENDDO

C       CALCULATE  DEFOCUS (ROUGH)
	DO I=1,20
	   A0(I) = 0
	   B(I)  = 0

	   DO J=1,20
	      C(I,J) = 0
	   ENDDO
	ENDDO

	IF (NUM == 1 .AND. NP(1) == 1 ) THEN
	   CALL RDPRM1S(A2,NOT_USED,'AMPLITUDE CONTRAST [RAD]',IRTFLG)

	   KF = KP(1,1) * KS
	   A0(1) = -(Y1(1,1) + A2 - 0.5 * CS * PI * LAMBDA**3 * KF**4) /
     &              (PI * LAMBDA * KF**2)
	   GOTO 130

	ELSE
	   DO  I=1,NUM
	      DO  J=1,NP(I)
	         KF             = KP(I,J) * KS
	         C(I,I)         = C(I,I) - PI * (LAMBDA*KF**2)**2
	         C(NUM+1,I)     = C(NUM+1,I) - PI*LAMBDA*KF**2
	         C(I,NUM+1)     = C(I,NUM+1) - LAMBDA*KF**2
	         C(NUM+1,NUM+1) = C(NUM+1,NUM+1) - 1.0
	         B(I)           = B(I) - (PI*(0.5*CS*LAMBDA**3*KF**4)-
     &                             Y1(I,J))*(LAMBDA*KF**2)
	         B(NUM+1)       = B(NUM+1) -
     &                            (PI*(0.5*CS*LAMBDA**3*KF**4)-Y1(I,J))
	      ENDDO
	   ENDDO

	   CALL MATINV(C,NUM+1,DET)

	   DO  I=1,NUM+1
	      DO  J=1,NUM+1
	         A0(I) = A0(I) + C(I,J) * B(J)
	      ENDDO
	   ENDDO
	ENDIF

	A2 = A0(NUM+1)

	IF ( A2 > 0.3 .OR. A2 < 0) THEN
C          USING DEEPEST GRADIENT METHOD TO CALCULATE AMPLITUDE

C          CALCULATE THE WEIGHT OF EACH POINTS
	   DO  I=1,NUM
	      DO  J=1,NP(I)
	         KF          = KP(I,J) * KS
	         WEIGHT(I,J) = PI*SQRT((2. * CS * LAMBDA**3 * KF**3)**2 +
     & 	                       (2. * A0(I)*LAMBDA * KF)**2) * 2. * KS
C	         WRITE(NOUT,*) WEIGHT(I,J)
	      ENDDO
	   ENDDO

	   DO  K=1,31

	      DO I=1,NUM
	         A1(I) = A0(I)
	      ENDDO

	      A2 = FLOAT(K) * 0.01
	      IF (K == 31) THEN
	         A2 = AM2
	      ENDIF

C             SET ITERATION STEP 
	      NSTEP = 0

C             CALCULATE VALUE OF Y2 
	      DO  I=1,NUM
	         DO  J=1,NP(I)
	            KF      = KP(I,J) * KS
	            Y2(I,J) = PI * (0.5 * CS * LAMBDA**3 * KF**4 -
     &                        A1(I)*LAMBDA*KF**2) - A2
	         ENDDO
	      ENDDO

C             SET INITIAL VALUE X**2(X,A)=SUM((Y1(I)-Y2(I))**2) */
	      X1 = 0
	      DO  I=1,NUM
	         DO J=1,NP(I)
                    X1 = X1+(Y1(I,J)-Y2(I,J))**2/WEIGHT(I,J)**2
	         ENDDO
	      ENDDO


C             CALCULATE THE VALUE OF Y2(I) 
999	      DO I=1,NUM
	         DA1(I) = 0.001 * A1(I)
	      ENDDO

	      DA2 = 0.001 * A2
	      DO  I=1,NUM
	         DO  J=1,NP(I)
	            KF      = KP(I,J)*KS
 	            Y2(I,J) = PI * (0.5*CS*LAMBDA**3*KF**4-
     &                        A1(I)*LAMBDA*KF**2)-A2
	         ENDDO
	      ENDDO

C             CALCULATE DX**2 / DA1 

	      DO  L=1,NUM
	         DO  I=1,NUM
	            IF(L == I) THEN
	               X = A1(L) + 0.1 * DA1(L)
	            ELSE
	               X = A1(I)
	            ENDIF

	            DO  J=1,NP(I)
	               KF      = KP(I,J) * KS
	               Y2(I,J) = PI * (0.5*CS*LAMBDA**3*KF**4 -
     &                           X * LAMBDA*KF**2)-A2
	            ENDDO
	         ENDDO

	         X2 = 0
	         DO  I=1,NUM
	            DO  J=1,NP(I)
	               X2 = X2+(Y1(I,J)-Y2(I,J))**2/WEIGHT(I,J)**2
	            ENDDO
	         ENDDO

	         DXA1(L) = (X2-X1)/(0.1*DA1(L))
	      ENDDO

C             CALCULATE DX**2 / DA2 
 	      DO  I=1,NUM
	         X = A1(I)
	         DO  J=1,NP(I)
	            KF      = KP(I,J) * KS
	            Y2(I,J) = PI * (0.5*CS*LAMBDA**3*KF**4 -
     &                        X * LAMBDA*KF**2)-(A2+0.1*DA2)
	         ENDDO
	      ENDDO

	      X2 = 0
	      DO  I=1,NUM
	         DO  J=1,NP(I)
	            X2 = X2 + (Y1(I,J) - Y2(I,J))**2 / WEIGHT(I,J)**2
	         ENDDO
	      ENDDO

	      DXA2 = (X2-X1) / (0.1*DA2)

C             CALCULATE THE SUM = SUM((DX**2/DAI*DAI)**2) 
	      SUM = 0

	      DO I=1,NUM
	         SUM = SUM + (DXA1(I)*DA1(I))**2
	      ENDDO

	      SUM = SQRT((DXA2*DA2)**2 + SUM)

	      DO I=1,NUM
                 A1(I) = A1(I) - DXA1(I) * DA1(I)**2 / SUM
	      ENDDO

           A2 = A2 - DXA2 * DA2**2 / SUM

C          CRITERIA FOR ITERATION 
	   DO  I=1,NUM
	      X = A1(I)
	      DO  J=1,NP(I)
                 KF      = KP(I,J)*KS
	         Y2(I,J) = PI*(0.5*CS*LAMBDA**3*KF**4 -
     &                     X*LAMBDA*KF**2)-A2
	      ENDDO
	   ENDDO

	   X2 = 0
	   DO  I=1,NUM
	      DO  J=1,NP(I)
	         X2 = X2 + (Y1(I,J)-Y2(I,J))**2 / WEIGHT(I,J)**2
	      ENDDO
	   ENDDO

	   IF ( (X1-X2) < 0) THEN

	       DO I=1,NUM
	          A1(I) = A1(I) + 0.5 * DXA1(I) * DA1(I)**2 / SUM
	       ENDDO

	       A2 = A2+0.5*DXA2*DA2**2 / SUM
C	       WRITE(NOUT,*) 'INITIAL A2=',FLOAT(K)*0.01,'STEP',NSTEP
C	       WRITE(NOUT,*) 'A1=',(A1(I),I=1,NUM),'OFFSET(RAD)=',A2,'X**2=',X1

	       AM2      = A2/0.01
	       NAM2     = INT(AM2+0.5)
	       NM(NAM2) = NM(NAM2)+1
	       XA2(K)   = A2

	       DO I=1,NUM
	          XA1(K,I) = A1(I)
	       ENDDO

	       XM(K) = X1

	    ELSE
C               SET PARAMETERS FOR NEXT STEP
	        X1    = X2
	        NSTEP = NSTEP + 1
	        GOTO 999
	   ENDIF
	ENDDO

C       FIND CONVERGE POINTS
	NN = 0
	DO I=1,30
	   IF (NM(I) > 5) NN = NN + 1
	ENDDO

	IF (NN > 1) THEN
C          THERE ARE TWO CONVERGE POINTS
           X0 = 999999999.0
           DO I=1,30
              IF (NM(I) > 5 .AND. XM(I) < X0) THEN
                 NN = I
                 X0 = XM(I)	
              ENDIF
           ENDDO

	   ELSE
C             THERE IS ONLY ONE CONVERGE POINT
              NN = 0
	      DO I=1,30
	        IF( NM(I) > NN .OR. (NM(I) == NN .AND.
     $           XM(I) < XM(NN))) NN = I
	      ENDDO
           ENDIF

	   A2 = XA2(NN)
	   DO I=1,NUM
	      A0(I) = XA1(NN,I)
	      A1(I) = A0(I)
	   ENDDO

	   X1 = XM(NN)

	ENDIF

        WRITE(NOUT,*) ' '	
130	WRITE(NOUT,*) ' DEFOCUS:', (A0(I), I=1,NUM)

	WRITE(NOUT,140) A2
140	FORMAT('  AMPLITUDE CONTRAST:',F10.6,/)
	
C       GENERATE A FILTER FILE FROM CTF W/O ENVELOPE FUNCTION
C	GATE VALUE IS 0.08

	GATE = 0.08

	CALL RDPRMC(CHO1,NUMC,.TRUE.,
     &     'DO YOU WANT TO GENERATE A FILTER? (Y/N)' ,NULL,IRTFLG)

	IF (CHO1 == 'Y') THEN
           WRITE(NOUT,*)'GENERATE FILTER IN THE INPUT FILE SEQUENCE'

           DO  I=1,NUM
              X = A0(I)
              DO  J=1,NX
                 KF = FLOAT(J)*KS
                 X1 = PI*(0.5*CS*LAMBDA**3*KF**4-X*LAMBDA*KF**2)-A2
                 F1 = SIN(X1)
                 IF (F1 > 0) THEN
                    F2(J) = -1.
                 ELSE
                    F2(J) = 1.
                 ENDIF
                 IF (F1 > GATE) THEN
                    F2(J) = -1.
                 ELSE
                    IF (F1 < -1.*GATE) THEN
                       F2(J) = 1.
                    ELSEIF( ABS(X1) > PI/2.) THEN
                       F2(J) = 0.
                    ENDIF
                 ENDIF
              ENDDO

              WRITE(NOUT,'(A,I0)')'  FILE #: ',I	
        
              IFORM = 1
              MAXIM = 0
              CALL OPFILEC(0,.TRUE.,OUTNAME,LUN2,'U',IFORM,NX,1,1,
     &                    MAXIM,'OUTPUT',.FALSE.,IRTFLG)
              CALL WRTLIN(LUN2,F2,NX,1)
              CLOSE(LUN2)
           ENDDO
	ENDIF
		
	END
