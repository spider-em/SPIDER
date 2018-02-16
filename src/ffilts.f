C++*******************************************************************
C
C FFILTS.F               USED OPFILE               NOV 00 ARDEAN LEITH
C                        ADDED BFACTOR             OCT 01 BILL BAXTER
C                        OPFILEC                   FEB 03 ARDEAN LEITH
C                        GAUSSIAN BUG              FEB 04 PP
C                        SAMPLED ADDED             MAR 07 C. RENKEN
C                        SQRT2M1 in BUTTERWORTH    JAN 18 ARDEAN LEITH
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
C  FFILTS(LUN,NX,NY,NZ,NXO)

C  PARAMETERS:
C        LUN         LOGICAL UNIT NUMBER OF FOURIER FILE TO BE FILTERED
C        LUN         LOGICAL UNIT NUMBER OF FOURIER FILE OUTPUT
C        NX,NY,NZ    DIMENSIONS OF FOURIER FILE
C        NXO         X DIMENSION OF ORIGINAL REAL-SPACE FILE
C
C--*******************************************************************

	SUBROUTINE FFILTS(LUN,LUNO,NX,NY,NZ,NXO)

	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

	COMMON                 B(1)
	EQUIVALENCE            (CB,B)
	COMPLEX             :: CB(1)
        CHARACTER *(MAXNAM) :: FILNAM
	CHARACTER*1         :: NULL = CHAR(0)

        REAL                :: SQRT2M1

	INTEGER, PARAMETER  :: LUNF = 27

951     WRITE(NOUT,1009)
 1009   FORMAT
     &     ('  1 - LOW-PASS,              2 - HIGH-PASS'         ,/,
     &      '  3 - GAUSS.  LOW-PASS,      4 - GAUSS.  HIGH-PASS' ,/,
     &      '  5 - FERMI                  6 - FERMI'             ,/,
     &      '  7 - BUTTER. LOW-PASS,      8 - BUTTER. HIGH-PASS' ,/,
     &      '  9 - REMEZ,                10 - B FACTOR',/,
     &      ' 11 - SAMPLED SPACE')

        IOPT = 1
        CALL RDPRI1S(IOPT,NOT_USED,'FILTER TYPE (1-11)',IRTFLG)
	IF (IOPT <1 .OR. IOPT > 11) THEN 
           CALL  ERRT(102,'UNKNOWN FILTER TYPE',IOPT)
           GOTO  951
        ENDIF


	IF (IOPT == 10)  THEN
C          B FACTOR FILTER *************************************************
           CALL BFACT(LUN,LUNO,NX,NY,NZ,NXO)
	   RETURN
	ENDIF


	IF (IOPT == 9)  THEN
C          REMEZ FILTER *************************************************
           MAXIM  = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUNF,'O',IFORM,NS1,NR1,NSL1,
     &                   MAXIM,'FILTER',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

	   IF (IFORM > 0 .OR.
     & 	      NS1.NE.NX .OR. NR1.NE.NY .OR. NSL1.NE.NZ) THEN
              CALL  ERRT(101,'INCONSISTENT IMAGE SIZE FOR FILTER',NE)
              CLOSE(LUNF)
              RETURN
	   ENDIF

	   NS2 = NX/2

	   DO K=1,NZ
	     DO J=1,NY
		NR = (K-1)*NY+J
		CALL REDLIN(LUN,CB,NX,NR )
		CALL REDLIN(LUNF,CB(NS2+1),NX,NR)

		DO I=1,NX/2
		   CB(I)=CB(I)*CB(I+NS2)
		ENDDO

		CALL WRTLIN(LUNO,CB,NX,NR)
	     ENDDO
	   ENDDO

           CLOSE(LUNF)
	   RETURN
	ENDIF


	IF (IOPT == 7 .OR. IOPT == 8)  THEN
C          BUTTERWORTH FILTER ****************************************
           EPS = 0.882
           AA  = 10.624

           CALL RDPRM2S(FP,FS,NOT_USED,
     &        'LOWER & UPPER LIMITING FREQ. (IN FREQ UNITS)', IRTFLG)

           ORD   = 2.  * ALOG10(EPS / SQRT(AA**2-1.0))
           ORD   = ORD / ALOG10(FP  / FS)

	   IF (FP > 0.5)  FP = FP / NXO   
	   IF (FS > 0.5)  FS = FS / NXO   

           PARM1 = FP / (EPS)**(2./ORD)

           GO TO 5768
	ENDIF


	IF (IOPT == 11)THEN

C          SAMPLED SPACE FILTER *********************************
           CALL RDPRI1S(NUMPRJ,NOT_USED,
     &                  'NUMBER OF PROJECTIONS',IRTFLG)
	   GOTO 5768	
        ENDIF



C       OTHER FILTERS ***********************************************

   	CALL RDPRM1S(PARM1,NOT_USED, 
     &        'FILTER RADIUS (IN FREQUENCY OR PIXEL UNITS)',IRTFLG)

	IF (PARM1 < 0.0 .OR. PARM1 > 0.5) PARM1 = PARM1 / NXO


	IF (IOPT == 5 .OR. IOPT == 6)  THEN
C          FERMI DISTRIBUTION FILTER ********************************
	   CALL RDPRM1S(TEMP,NOT_USED,
     &                 'TEMPERATURE (0=CUTOFF)',IRTFLG)

C          EXPONENTIAL FOR HIGH-PASS OPTION
	   IF (IOPT == 6) TEMP = -TEMP
	ENDIF


5768	NS2 = NX/2
	NR2 = NY/2
	NL2 = NZ/2

	X1  = FLOAT(NXO/2)**2
	Y1  = FLOAT(NR2)**2

        !write(6,*) ' nx,nxo,x1: ',nx,NXO,x1

	IF (NZ == 1) THEN
	   Z1   = 1.0
	ELSE
	   Z1   = FLOAT(NL2)**2
	ENDIF

	PARM    = PARM1**2

        SQRT2M1 = (SQRT(2.0)-1)

	DO K=1,NZ
	   IZ = K-1
	   IF (IZ  >  NL2)  IZ = IZ-NZ

	   DO J=1,NY
	      IY = J-1
	      IF (IY  >  NR2)  IY = IY-NY
	      NR = J+(K-1)*NY
	      CALL  REDLIN(LUN,B,NX,NR )

	      DO  2  I=1,NS2
	         IX=I-1
		 IF (IOPT  ==  11) GOTO 800

	         GOTO(100,200,300,400,500,500,600,700),IOPT


C       LOWPASS ******************************************************
100	IF (0.25*(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1+FLOAT(IZ*IZ)/Z1)
     &	        >  PARM)  CB(I) = CMPLX(0.0,0.0)
	GOTO  2

C       HIGH PASS ****************************************************
200	IF((IX.NE.0 .OR. IY.NE.0 .OR. IZ.NE.0) .AND.
     &      0.25*(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1+FLOAT(IZ*IZ)/Z1)
     &	       .LE. PARM)  CB(I) = CMPLX(0.0,0.0)
	GOTO  2

C       GAUSSIAN LOW PASS ********************************************
300  	F =0.125*(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1+FLOAT(IZ*IZ)/Z1)/PARM
	IF (F .LT. 16.0)  THEN
	   CB(I) = CB(I) * EXP(-F)
	ELSE
	   CB(I) = CMPLX(0.0,0.0)
	ENDIF
	GOTO  2

C       GAUSSIAN HIGH PASS ********************************************
400	IF (IX.NE.0 .OR. IY.NE.0 .OR. IZ.NE.0)  THEN
  	   F=0.125*(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1+FLOAT(IZ*IZ)/Z1)/
     &         PARM
	   IF (F .LT. 16.0)  THEN
	      CB(I) = CB(I) * (1.0-EXP(-F))
	   ENDIF
	ENDIF
	GOTO  2

C       FERMI DISTRIBUTION FILTER *************************************
500	ARG = (0.5*SQRT(
     &	   FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1+FLOAT(IZ*IZ)/Z1)-PARM1)/TEMP
	ARG = AMIN1(AMAX1(ARG,-10.0),10.0)
	IF (IOPT == 6.AND.IX.NE.0.AND.IY.NE.0.AND.IZ.NE.0) GOTO 2
	CB(I) = CB(I) * (1.0/(1.0+EXP(ARG)))
	GOTO 2

C       BUTTERWORTH  LOWPASS FILTER **********************************
600	ARG = 0.5*SQRT(
     &	   FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1+FLOAT(IZ*IZ)/Z1)
	CB(I) = CB(I) * SQRT(1.0/(1.0+SQRT2M1*(ARG/PARM1)**ORD))
	GOTO 2

C       BUTTERWORTH HIGHPASS FILTER **********************************
700     IF (IX.NE.0 .OR. IY.NE.0 .OR. IZ.NE.0) THEN
           ARG = 0.5*SQRT(
     &	      FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1+FLOAT(IZ*IZ)/Z1)
	   CB(I) = CB(I)*(1.0-SQRT(1.0/(1.0+SQRT2M1*(ARG/PARM1)**ORD)))
	ENDIF
        GOTO 2

C          SAMPLED SPACE FILTER **************************************
800        ARG = SQRT(FLOAT(IX*IX + IY*IY + IZ*IZ))
	   IF (ARG  ==  0.0) ARG = 1.0

           F = 3*NUMPRJ*((ARG+0.5)**2 - (ARG-0.5)**2)/
     &            (4*((ARG+0.5)**3 - (ARG-0.5)**3))

	   IF (F  >  1) F = 1.0
           CB(I) = CB(I) * F
 


2	   CONTINUE
	   CALL WRTLIN(LUNO,B,NX,NR)

	   ENDDO
	ENDDO

	END
