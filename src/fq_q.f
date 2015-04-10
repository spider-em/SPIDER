
C++*********************************************************************
C
C FQ_Q.F                                        12/22/94  
C                RDPRAF REMOVED                 DEC 2005 ARDEAN LEITH 
C                COSINE FILTER                  JUL 2012 G. KISHCHENKO 
C                INCORE WITHOUT LUNS            OCT 2012 ARDEAN LEITH 
C                FREQ + PIXELS                  NOV 2012 G. KISHCHENKO 
C                PARM2 BUG                      DEC 2012 ARDEAN LEITH 
C                FREQ UNIT CUTOFF = 1           AUG 2013 ARDEAN LEITH 
C                RAISED SINC                    FEB 2014 ARDEAN LEITH 
C                PARM2 BUG                      NOV 2014 ARDEAN LEITH 
C        
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014, Health Research Inc.,                         *
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
C FQ_Q.F(LUN,LUNO,B,LSD,N2X,N2Y,NX,NY,IOPT)
C
C PURPOSE: QUICK FILTERING OF REAL-SPACE IMAGE FILE BY FFT
C
C PARAMETERS:
C        LUN         I/O UNIT OF REAL-SPACE FILE TO BE FILTERED
C        LUNO        I/O UNIT OF REAL-SPACE OUTPUT FILE 
C        B           BUFFER 
C        NX,NY       DIMENSIONS OF REAL-SPACE FILE
C        N2X=2*NX    AT LEAST
C        N2Y=2*NY    "    "
C        IOPT        TYPE OF FILTER
C
C NOTE:  APPEARS TO HAVE UNDOCUMENTED AND UNTESTED ELLIPTICAL 
C        FILTRATION FOR OPTIONS: 1,2,3,4.   THIS WILL GIVE BAD 
C        ERRORS IF PERSON ENTERS MORE THAN ONE VALUE ON THE INPUT
C        PARAMETER LINE AS IT IS USED FOR ELLIPSES. al nov 2014
C
C        'FF' (ffilts.f)            SETS X1 TO: (NX  / 2)**2 BUT
C        'FQ' (four_fq.f or fq_q.f) SETS X1 TO: (NXF / 2)**2 WHERE
C             NX  IS X DIMENSION OF POSSIBLY PADDED IMAGE
C             NXF IS SLIGHTLY LARGER DUE TO MIXED RADIX FOURIER PAD
C        SO THEY GIVE SLIGHTLY DIFFERENT RESULTS.  I SUSPECT THAT
C        'FF' IS ACTUALLY CORRECT?
C
C23456789012345678901234567890123456789012345678901234567890123456789012  
C--*******************************************************************

        SUBROUTINE FQ_BUF(IOPTT,BFPS,PARM1T,PARM2T,TEMPT,
     &                    B, LSD,N2X,N2Y, NX,NY, IRTFLG)
        
        IMPLICIT NONE
	INCLUDE 'CMBLOCK.INC'

        INTEGER          :: IOPTT	
	REAL             :: BFPS(4)
	REAL             :: FP, FS
	REAL             :: PARM1T,PARM2T,TEMPT
	REAL             :: B(LSD,N2Y)
	INTEGER          :: LSD,N2X,N2Y
	INTEGER          :: NX,NY,IRTFLG

	REAL             :: PARM1,PARM2,TEMP,PARM,PARM22,X1,Y1
	REAL             :: F,F2,FPE,FSE,ORDT
	DOUBLE PRECISION :: AVE
	REAL             :: ORD,IPS,AA,PARMT,EPS,AVG
        INTEGER          :: IOPT,J,I,INV,IX,IY,NR2
	
        REAL, PARAMETER  :: PI = 3.14159265358979323846

        IOPT  = IOPTT
        PARM1 = PARM1T
        PARM2 = PARM2T
        TEMP  = TEMPT

	IF (N2X .NE. NX .AND. N2Y .NE. NY)  THEN
C          BORDER PADDING

	   AVE = (SUM(B(1:NX,1))   + SUM(B(1:NX,NY)) +
     &	          SUM(B(1,2:NY-1)) + SUM(B(NX,2:NY-1)) ) / 
     &		  REAL(2*(NX+NY)-4)

c$omp      parallel do private(i,j)
	   DO J=1,N2Y
	      DO I=NX+1,N2X
	         B(I,J) = AVE
	      ENDDO
	   ENDDO

c$omp      parallel do private(i,j)
	   DO J=NY+1,N2Y
	      DO I=1,NX
	         B(I,J) = AVE
	      ENDDO
	   ENDDO
	ENDIF

C       FORWARD FFT, 2XFFT PADDED
	INV = 1
	CALL FMRS_2(B,N2X,N2Y,INV)
	IF (INV == 0) THEN
	   IRTFLG = 1
	   RETURN
	ENDIF

C       BUTTERWORTH FILTER ***********************************

	IF (IOPT == 7 .OR. IOPT == 8 .OR.  
     &      IOPT == 9 .OR. IOPT == 10)  THEN

	   EPS = 0.882
	   AA  = 10.624
	   ORD = 2.0 * ALOG10(EPS / SQRT(AA**2-1.0) )

	   IF (BFPS(3) == 0.0 .AND. BFPS(4) == 0.0) THEN
	      ORD   = ORD / ALOG10(BFPS(1) / BFPS(2))
	      PARM1 = BFPS(1) / (EPS)**(2./ORD)
	   ELSE
C             BUTTERWORTH FILTER ELLIPTIC FILTER:
C             LOW-PASS  IOPT=11,  HIGH-PASS IOPT=12
	      IOPT = IOPT + 4
           ENDIF

	ELSE

	   IF (PARM1 <  0.0 .OR. PARM1 > 1.0) PARM1 = 0.5*PARM1/(NX/2)
	   IF (PARM2 == 0.0)                  PARM2 = PARM1
	   IF (PARM2 <  0.0 .OR. PARM2 > 1.0) PARM2 = 0.5*PARM2/(NY/2)

	   IF (IOPT == 5 .OR. IOPT == 6)  THEN

C             FERMI DISTRIBUTION FILTER ********************

C             EXPONENTIAL FOR HIGH-PASS OPTION
	      IF (IOPT == 6) TEMP = -TEMP
	   ENDIF
	ENDIF

	NR2    = N2Y / 2
	X1     = FLOAT(N2X/2)**2
	Y1     = FLOAT(NR2)  **2
	PARM   = PARM1**2
	PARM22 = PARM2**2

C       KEEP ZERO TERM FOR HIGH PASS OPTIONS
	AVG = B(1,1)

c$omp   parallel do private(i,j,ix,iy,f,fpe,fse,ordt,parmt,f2)
	DO J=1,N2Y
	   IY = (J-1)
	   IF (IY > NR2) IY = IY-N2Y

	   DO I=1,LSD,2
	      IX = (I-1)/2

	      IF (IOPT == 1) THEN
C                LOWPASS *************************************
                 IF (0.25*(FLOAT(IX*IX)/X1/PARM +
     &                     FLOAT(IY*IY)/Y1/PARM22) > 1.0) THEN
	            B(I,J)   = 0.0
	            B(I+1,J) = 0.0
	         ENDIF

	      ELSEIF (IOPT == 2) THEN	
C                HIGH PASS ***********************************
	         IF ( (IX.NE.0 .OR. IY.NE.0) .AND.
     &                0.25*(FLOAT(IX*IX)/X1/PARM + 
     &                      FLOAT(IY*IY)/Y1/PARM22) <= 1.0) THEN
	            B(I,J)   = 0.0
	            B(I+1,J) = 0.0
	         ENDIF

	      ELSEIF(IOPT == 3)  THEN
C                GAUSSIAN LOW PASS ***************************
	         F = 0.125*(FLOAT(IX*IX)/X1/PARM +
     &                      FLOAT(IY*IY)/Y1/PARM22)
	         IF (F < 16.0)  THEN
	            F        = EXP(-F)
                    B(I,J)   = B(I,J)  *F
                    B(I+1,J) = B(I+1,J)*F
	         ELSE
                    B(I,J)   = 0.0
                    B(I+1,J) = 0.0
	         ENDIF

	      ELSEIF (IOPT==4)  THEN	
C                GAUSSIAN HIGH PASS **************************

	         IF (IX .NE. 0 .OR. IY .NE. 0)  THEN
	            F = 0.125*(FLOAT(IX*IX)/X1/PARM +
     &                         FLOAT(IY*IY)/Y1/PARM22)
	            IF (F < 16.0)  THEN
	               F        = 1.0 - EXP(-F)
                       B(I,J)   = B(I,J)  *F
                       B(I+1,J) = B(I+1,J)*F
	            ENDIF
	         ENDIF

	      ELSEIF (IOPT == 5 .OR. IOPT == 6)  THEN
C                FERMI DISTRIBUTION FILTER *******************
	      
	         F = (0.5*SQRT(FLOAT(IX*IX)/X1 +
     &                         FLOAT(IY*IY)/Y1)-PARM1) / TEMP
	         F        = AMIN1(AMAX1(F,-10.0), 10.0)
                 F        = (1.0/(1.0+EXP(F)))

                 B(I,J)   = B(I,J)  *F
                 B(I+1,J) = B(I+1,J)*F

 	      ELSEIF (IOPT == 7) THEN
C                BUTTERWORTH LOWPASS FILTER ******************

 	         F        = 0.5*SQRT(FLOAT(IX*IX)/X1 +
     &                               FLOAT(IY*IY)/Y1)

 	         F        = SQRT(1.0/(1.0+(F/PARM1)**ORD))
                 B(I,J)   = B(I,J)  *F
                 B(I+1,J) = B(I+1,J)*F

 	      ELSEIF (IOPT == 8) THEN
C                BUTTERWORTH HIGHPASS FILTER *****************
	
                 IF (IX.NE.0 .OR. IY.NE. 0) THEN
 	            F = 0.5*SQRT(FLOAT(IX*IX)/X1 +
     &                           FLOAT(IY*IY)/Y1)
 	            F = (1.0-SQRT(1.0/(1.0+(F/PARM1)**ORD)))

                    B(I,J)   = B(I,J)*F
                    B(I+1,J) = B(I+1,J)*F
 	         ENDIF


 	      ELSEIF (IOPT == 9) THEN
C                RAISED COSINE LOWPASS FILTER ******************

	         F = 0.5*SQRT(FLOAT(IX*IX)/X1 +
     &                        FLOAT(IY*IY)/Y1)
C	         F = (F-BFPS(1)) / (BFPS(2)-BFPS(1))
	         F = (F-FP) / (FS-FP)

                 IF (F < 0) THEN
	            F2 = 1
                 ELSEIF (F > 1) THEN
	            F2 = 0
                 ELSE
	            F2 = 0.5 * (COS(PI*F)+1)
	         ENDIF

                 B(I,J)   = B(I,J)  *F2
                 B(I+1,J) = B(I+1,J)*F2

	      ELSEIF (IOPT == 10) THEN
C                RAISED COSINE HIGHPASS FILTER ******************

	         F = 0.5*SQRT(FLOAT(IX*IX)/X1 +
     &                        FLOAT(IY*IY)/Y1)
	         F = (F-BFPS(1)) / (BFPS(2)-BFPS(1))

                 IF (F < 0) THEN
                    F2 = 0
                 ELSEIF (F > 1) THEN
	            F2 = 1
                 ELSE
	            F2 = 0.5 * (-COS(PI*F)+1)
	         ENDIF

                 B(I,J)   = B(I,J)  *F2
                 B(I+1,J) = B(I+1,J)*F2

	      ELSEIF (IOPT == 11) THEN 
C                BUTTERWORTH ELLIPTIC LOWPASS FILTER *********	
C                CALCULATE EFFECTIVE FP AND FS IN A GIVEN 
C                DIRECTION ON THE PLANE

                 IF (IX.NE.0 .OR. IY.NE.0) THEN
	            FPE = ATAN2(BFPS(1)*SQRT(FLOAT(IY*IY)/Y1),
     &                          BFPS(3)*SQRT(FLOAT(IX*IX)/X1))
                    FPE = SQRT((BFPS(1)*COS(FPE))**2 + 
     &                         (BFPS(3)*SIN(FPE))**2)

	            FSE = ATAN2(BFPS(2)*SQRT(FLOAT(IY*IY)/Y1),
     &                          BFPS(4)*SQRT(FLOAT(IX*IX)/X1))
                    FSE = SQRT((BFPS(2)*COS(FSE))**2 + 
     &                         (BFPS(4)*SIN(FSE))**2)

	            ORDT     = ORD/ALOG10(FPE/FSE)
	            PARMT    = FPE/(EPS)**(2./ORDT)
	            F        = 0.5*SQRT(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1)
	            F        = SQRT(1.0/(1.0+(F/PARMT)**ORDT))
                    B(I,J)   = B(I,J)  *F
                    B(I+1,J) = B(I+1,J)*F
	         ENDIF

	      ELSEIF (IOPT == 12) THEN
C                BUTTERWORTH ELLIPTIC HIGHPASS FILTER *********	

                 IF (IX .NE. 0 .OR. IY.NE. 0) THEN
	            FPE = ATAN2(BFPS(1)*SQRT(FLOAT(IY*IY)/Y1),
     &                          BFPS(3)*SQRT(FLOAT(IX*IX)/X1))
                    FPE = SQRT((BFPS(1)*COS(FPE))**2 +
     &                         (BFPS(3)*SIN(FPE))**2)

	            FSE = ATAN2(BFPS(2)*SQRT(FLOAT(IY*IY)/Y1),
     &                          BFPS(4)*SQRT(FLOAT(IX*IX)/X1))
                    FSE = SQRT((BFPS(2)*COS(FSE))**2 + 
     &                         (BFPS(4)*SIN(FSE))**2)

	            ORDT     = ORD / ALOG10(FPE/FSE)
	            PARMT    = FPE / (EPS)**(2./ORDT)
	            F        = 0.5*SQRT(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1)
	            F        = (1.0-SQRT(1.0/(1.0+(F/PARMT)**ORDT)))
                    B(I,J)   = B(I,J)  *F
                    B(I+1,J) = B(I+1,J)*F
	         ENDIF

	      ELSEIF (IOPT == 13) THEN
C                RAISED SINC WINDOW **************************
	         F = 0.5 * SQRT(FLOAT(IX*IX)/X1/PARM +
     &                          FLOAT(IY*IY)/Y1/PARM22)
                 IF (F <= 0.0001) THEN
	            F2 = 1
                 ELSEIF (F >= 1.0) THEN
	            F2 = 0
                 ELSE
	            F2 = SIN(PI*F)/(PI*F)
	         ENDIF
                 B(I,J)   = B(I,J)  *(1+9*F2)
                 B(I+1,J) = B(I+1,J)*(1+9*F2)
              ENDIF
	   ENDDO
	ENDDO

C       RESTORE ZERO TERM FOR HIGH PASS OPTIONS
	IF (IOPT == 2 .OR. IOPT == 4 .OR. IOPT == 6 .OR. IOPT == 8) 
     &     B(1,1) = AVG

C       REVERSE FFT, 2X PADDED
	INV = -1
	CALL FMRS_2(B,N2X,N2Y,INV)

        IRTFLG = 0

        END



C       ----------------------  FQ_Q --------------------------------


        SUBROUTINE FQ_Q(IOPT,LUN,LUNO, B, LSD,N2X,N2Y, NX,NY, IRTFLG)
        
	INCLUDE 'CMBLOCK.INC'
	
	INTEGER          :: IOPT
	INTEGER          :: LUN,LUNO
	REAL             :: B(LSD,N2Y)
	INTEGER          :: LSD,N2X,N2Y
	INTEGER          :: NX,NY,IRTFLG

	DOUBLE PRECISION :: AVE
	REAL             :: BFPS(4)
	REAL             :: FP, FS
	
        REAL, PARAMETER  :: PI = 3.14159265358979323846

C       TO SET THEM TO SOMETHING.
	PARM1 = 0.0
	PARM2 = 0.0

        IF (LUN > 0) THEN
C          READ  IMAGE
	   DO I=1,NY
 	      CALL  REDLIN(LUN,B(1,I),NX,I)
	   ENDDO
        ENDIF

C       BORDER PADDING
	IF (N2X .NE. NX .AND. N2Y .NE. NY)  THEN
	   AVE = (SUM(B(1:NX,1))   + SUM(B(1:NX,NY)) +
     &	          SUM(B(1,2:NY-1)) + SUM(B(NX,2:NY-1)) ) / 
     &		  REAL(2*(NX+NY)-4)

c$omp      parallel do private(i,j)
	   DO J=1,N2Y
	      DO I=NX+1,N2X
	         B(I,J) = AVE
	      ENDDO
	   ENDDO

c$omp      parallel do private(i,j)
	   DO J=NY+1,N2Y
	      DO I=1,NX
	         B(I,J) = AVE
	      ENDDO
	   ENDDO
	ENDIF

C       FORWARD FFT
	INV = 1
	CALL FMRS_2(B,N2X,N2Y,INV)
	IF (INV == 0) THEN
	   IRTFLG = 1
	   RETURN
	ENDIF

	IF (IOPT == 7  .OR. IOPT == 8 .OR.  
     &      IOPT == 9  .OR. IOPT == 10)  THEN
C          BUTTERWORTH FILTER OR  RAISED COSINE FILTER ***************

	   NMAX = 4
	   BFPS = 0.0
	   CALL RDPRA(
     &        'LOWER & UPPER LIMITING FREQ. (IN FREQ OR PIXEL UNITS)',
     &        NMAX,0,.FALSE.,BFPS,NGOT,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

	   EPS = 0.882
	   AA  = 10.624
	   ORD = 2.0 * ALOG10(EPS / SQRT(AA**2-1.0) )

	   IF (BFPS(1) > 1.0) THEN
              FP = BFPS(1) / NX
              FS = BFPS(2) / NX
           ELSE
	      FP = BFPS(1)
	      FS = BFPS(2)
	   ENDIF

	   IF (BFPS(3) == 0.0 .AND. BFPS(4) == 0.0) THEN
C             BUTTERWORTH CIRCULAR FILTER:

	      ORD   = ORD / ALOG10(FP / FS)
	      PARM1 = FP  / (EPS)**(2./ORD)

	   ELSE
C             BUTTERWORTH ELLIPTIC FILTER:
C             LOW-PASS  IOPT=11,  HIGH-PASS IOPT=12
	      IOPT = IOPT + 4
           ENDIF

	ELSE
           PARM1T = 0.25
           PARM2T = -9999999
  	   CALL RDPRM2S(PARM1T,PARM2T,NOT_USED,
     &        'FILTER RADIUS (IN FREQUENCY OR PIXEL UNITS)',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           PARM1 = PARM1T
	   IF (PARM1T < 0.0 .OR. PARM1T > 1.0) PARM1 = 0.5*PARM1/(NX/2)

	   IF (PARM2T == -9999999) THEN
              PARM2 = PARM1

	   ELSEIF  (PARM2T < 0.0 .OR. PARM2T > 1.0) THEN
              PARM2 = 0.5 * PARM2T/(NY/2)
           ENDIF

	   IF (IOPT == 5 .OR. IOPT == 6)  THEN

C             FERMI DISTRIBUTION FILTER ********************
	      CALL RDPRM1S(TEMP,NOT_USED,
     &                     'TEMPERATURE (0=CUTOFF)',IRTFLG)

C             EXPONENTIAL FOR HIGH-PASS OPTION
	      IF (IOPT == 6) TEMP = -TEMP
	   ENDIF
	ENDIF

	NR2    = N2Y / 2
	X1     = FLOAT(N2X/2)**2
	Y1     = FLOAT(NR2)  **2
	PARM   = PARM1**2
	PARM22 = PARM2**2

C       KEEP ZERO TERM FOR HIGH PASS OPTIONS
	AVG = B(1,1)

c$omp   parallel do private(i,j,ix,iy,f,fpe,fse,ordt,parmt,f2)
	DO J=1,N2Y

	   IY = (J-1)
	   IF (IY > NR2) IY = IY-N2Y

	   DO I=1,LSD,2
	      IX = (I-1)/2

	      IF (IOPT == 1) THEN
C                LOWPASS *************************************
                 IF (0.25*(FLOAT(IX*IX)/X1/PARM +
     &                     FLOAT(IY*IY)/Y1/PARM22) > 1.0) THEN
	            B(I,J)   = 0.0
	            B(I+1,J) = 0.0
	         ENDIF

	      ELSEIF (IOPT == 2) THEN	
C                HIGH PASS ***********************************
	         IF ( (IX.NE.0 .OR. IY.NE.0) .AND.
     &                0.25*(FLOAT(IX*IX)/X1/PARM + 
     &                      FLOAT(IY*IY)/Y1/PARM22) <= 1.0) THEN
	            B(I,J)   = 0.0
	            B(I+1,J) = 0.0
	         ENDIF

	      ELSEIF(IOPT == 3)  THEN
C                GAUSSIAN LOW PASS ***************************
	         F = 0.125*(FLOAT(IX*IX)/X1/PARM +
     &                      FLOAT(IY*IY)/Y1/PARM22)
	         IF (F < 16.0)  THEN
	            F        = EXP(-F)
                    B(I,J)   = B(I,J)  *F
                    B(I+1,J) = B(I+1,J)*F
	         ELSE
                    B(I,J)   = 0.0
                    B(I+1,J) = 0.0
	         ENDIF

	      ELSEIF (IOPT==4)  THEN	
C                GAUSSIAN HIGH PASS **************************

	         IF (IX .NE. 0 .OR. IY .NE. 0)  THEN
	            F=0.125*(FLOAT(IX*IX)/X1/PARM +
     &                       FLOAT(IY*IY)/Y1/PARM22)
	            IF (F < 16.0)  THEN
	               F        = 1.0 - EXP(-F)
                       B(I,J)   = B(I,J)  *F
                       B(I+1,J) = B(I+1,J)*F
	            ENDIF
	         ENDIF

	      ELSEIF (IOPT == 5 .OR. IOPT == 6)  THEN
C                FERMI DISTRIBUTION FILTER *******************
	      
	         F = (0.5*SQRT(FLOAT(IX*IX)/X1 +
     &                         FLOAT(IY*IY)/Y1)-PARM1) / TEMP
	         F        = AMIN1(AMAX1(F,-10.0), 10.0)
                 F        = (1.0/(1.0+EXP(F)))

                 B(I,J)   = B(I,J)  *F
                 B(I+1,J) = B(I+1,J)*F

 	      ELSEIF (IOPT == 7) THEN
C                BUTTERWORTH LOWPASS FILTER ******************

 	         F        = 0.5*SQRT(FLOAT(IX*IX)/X1 +
     &                               FLOAT(IY*IY)/Y1)

 	         F        = SQRT(1.0/(1.0+(F/PARM1)**ORD))
                 B(I,J)   = B(I,J)  *F
                 B(I+1,J) = B(I+1,J)*F

 	      ELSEIF (IOPT == 8) THEN
C                BUTTERWORTH HIGHPASS FILTER *****************
	
                 IF (IX.NE.0 .OR. IY.NE. 0) THEN
 	            F = 0.5 * SQRT(FLOAT(IX*IX)/X1 +
     &                             FLOAT(IY*IY)/Y1)
 	            F = (1.0-SQRT(1.0/(1.0+(F/PARM1)**ORD)))

                    B(I,J)   = B(I,J)*F
                    B(I+1,J) = B(I+1,J)*F
 	         ENDIF

 	      ELSEIF (IOPT == 9) THEN
C                RAISED COSINE LOWPASS FILTER ******************

	         IF (BFPS(1) > 1.0) THEN
                    FP = BFPS(1)/NX
                    FS = BFPS(2)/NX
                 ELSE
	            FP = BFPS(1)
	            FS = BFPS(2)
	         ENDIF

	         F = 0.5 * SQRT(FLOAT(IX*IX)/X1 +
     &                          FLOAT(IY*IY)/Y1)
C	         F = (F-BFPS(1)) / (BFPS(2)-BFPS(1))
	         F = (F-FP) / (FS-FP)

                 IF (F < 0) THEN
	            F2 = 1
                 ELSEIF (F > 1) THEN
	            F2 = 0
                 ELSE
	            F2 = 0.5 * (COS(PI*F)+1)
	         ENDIF

                 B(I,J)   = B(I,J)  *F2
                 B(I+1,J) = B(I+1,J)*F2

	      ELSEIF (IOPT == 10) THEN
C                RAISED COSINE HIGHPASS FILTER ******************


	         IF (BFPS(1) > 1.0) THEN
                    FP = BFPS(1) / NX
                    FS = BFPS(2) / NX
                 ELSE
                    FP = BFPS(1)
	            FS = BFPS(2)
	         ENDIF

	         F = 0.5 * SQRT(FLOAT(IX*IX)/X1 +
     &                          FLOAT(IY*IY)/Y1)
	         F = (F-FP) / (FS-FP)

                 IF (F < 0) THEN
                    F2 = 0
                 ELSEIF (F > 1) THEN
	            F2 = 1
                 ELSE
	            F2 = 0.5 * (-COS(PI*F)+1)
	         ENDIF

                 B(I,J)   = B(I,J)  *F2
                 B(I+1,J) = B(I+1,J)*F2

	      ELSEIF (IOPT == 11) THEN 
C                BUTTERWORTH ELLIPTIC LOWPASS FILTER *********	
C                CALCULATE EFFECTIVE FP AND FS IN A GIVEN 
C                DIRECTION ON THE PLANE

                 IF (IX .NE. 0 .OR. IY .NE. 0) THEN
	            FPE = ATAN2(BFPS(1)*SQRT(FLOAT(IY*IY)/Y1),
     &                          BFPS(3)*SQRT(FLOAT(IX*IX)/X1))
                    FPE = SQRT((BFPS(1)*COS(FPE))**2 + 
     &                         (BFPS(3)*SIN(FPE))**2)

	            FSE = ATAN2(BFPS(2)*SQRT(FLOAT(IY*IY)/Y1),
     &                          BFPS(4)*SQRT(FLOAT(IX*IX)/X1))
                    FSE = SQRT((BFPS(2)*COS(FSE))**2 + 
     &                         (BFPS(4)*SIN(FSE))**2)

	            ORDT     = ORD/ALOG10(FPE/FSE)
	            PARMT    = FPE/(EPS)**(2./ORDT)
	            F        = 0.5*SQRT(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1)
	            F        = SQRT(1.0/(1.0+(F/PARMT)**ORDT))
                    B(I,J)   = B(I,J)  *F
                    B(I+1,J) = B(I+1,J)*F
	         ENDIF

	      ELSEIF (IOPT == 12) THEN
C                BUTTERWORTH ELLIPTIC HIGHPASS FILTER *********	

                 IF (IX .NE. 0 .OR. IY.NE. 0) THEN
	            FPE = ATAN2(BFPS(1)*SQRT(FLOAT(IY*IY)/Y1),
     &                          BFPS(3)*SQRT(FLOAT(IX*IX)/X1))
                    FPE = SQRT((BFPS(1)*COS(FPE))**2 +
     &                         (BFPS(3)*SIN(FPE))**2)

	            FSE = ATAN2(BFPS(2)*SQRT(FLOAT(IY*IY)/Y1),
     &                          BFPS(4)*SQRT(FLOAT(IX*IX)/X1))
                    FSE = SQRT((BFPS(2)*COS(FSE))**2 + 
     &                         (BFPS(4)*SIN(FSE))**2)

	            ORDT     = ORD / ALOG10(FPE/FSE)
	            PARMT    = FPE / (EPS)**(2./ORDT)
	            F        = 0.5*SQRT(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1)
	            F        = (1.0-SQRT(1.0/(1.0+(F/PARMT)**ORDT)))
                    B(I,J)   = B(I,J)  *F
                    B(I+1,J) = B(I+1,J)*F
	         ENDIF

	      ELSEIF (IOPT == 13) THEN
C                RAISED SINC WINDOW **************************
	         F = 0.5 * SQRT(FLOAT(IX*IX)/X1/PARM +
     &                          FLOAT(IY*IY)/Y1/PARM22)
                 IF (F <= 0.0001) THEN
	            F2 = 1
                 ELSEIF (F >= 1.0) THEN
	            F2 = 0
                 ELSE
	            F2 = SIN(PI*F)/(PI*F)
	         ENDIF
                 B(I,J)   = B(I,J)  *(1+9*F2)
                 B(I+1,J) = B(I+1,J)*(1+9*F2)
              ENDIF
	   ENDDO
	ENDDO

C       RESTORE ZERO TERM FOR HIGH PASS OPTIONS
	IF (IOPT == 2 .OR. IOPT == 4 .OR. 
     &      IOPT == 6 .OR. IOPT == 8) 
     &      B(1,1) = AVG

C       REVERSE FFT 
	INV = -1
	CALL FMRS_2(B,N2X,N2Y,INV)

	IF (LUNO > 0) THEN
C          WRITE  IMAGE
           DO I=1,NY
 	      CALL  WRTLIN(LUNO,B(1,I),NX,I)
	   ENDDO
        ENDIF

        IRTFLG = 0

        END
