C++*********************************************************************
C
C TRAFC3.F
C              OPFILEC                             FEB 03  ARDEAN LEITH
C              CMLIMIT & PROMPTS                   NOV 10  ARDEAN LEITH
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
C
C changed 9/5/94 to include complete sign reversal and intelligible
C Gaussian parameter jf
C  FIXED, PP 3/2/95
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE TRAFC3(LUN)

         INCLUDE 'CMBLOCK.INC'
	 INCLUDE 'CMLIMIT.INC'

	 COMPLEX                :: B
         COMMON        B(1)

         CHARACTER(LEN=MAXNAM)  ::  FILNAM
         REAL                   :: LAMBDA,KM,SIGN
      	 CHARACTER(LEN=1)       :: NULL = CHAR(0)  

         DATA PI/3.1415926/


         CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
         IF (IRTFLG .EQ. -1) RETURN

	CALL RDPRM(CS,NOT_USED,'SPHERICAL ABERRATION CS[MM]')
           IF (CS < 0.0001)    CS = 0.0001

         CALL RDPRM2(DZ,LAMBDA,NOT_USED, 'DEFOCUS [A], LAMBDA [A]')

         CALL RDPRMI(NSAM,NDUM,NOT_USED,
     &               'NUMBER OF SPATIAL FREQ. POINTS')

         CALL RDPRM(KM,NOT_USED,'MAXIMUM SPATIAL FREQUENCY [1/A]')

         CALL RDPRM2(Q,DS,NOT_USED,
     &       'SOURCE SIZE [1/A], DEFOCUS SPREAD [A]')

         CALL RDPRM2(DZA,AZZ,NOT_USED,'ASTIGMATISM [A], AZIMUTH [DEG]')

         CALL RDPRM2(WGH,ENV,NOT_USED,
     &      'AMPL CONTRAST RATIO [0-1], GAUSSIAN ENV HALFW [1/A]')
 
         ENV = 1. / ENV**2

         CALL  RDPRM(SIGN,NOT_USED,'SIGN [+1 OR -1]')

         IFORM  = -7
         NROW   = NSAM
         NSLICE = NROW
	 IF (MOD(NSAM,2) .EQ. 0)  THEN
            IFORM = -22
	    LSM   = NSAM+2
	 ELSE
	    IFORM = -21
	    LSM   = NSAM+1
	 ENDIF

         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,LSM,NROW,NSLICE,
     &               MAXIM,' ',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         SC = KM / FLOAT(NSAM/2)

C        IE=0 SELECTS TRANSFER FUNCTION OPTION IN SUBROUTINE TFD
	 IE  = 0
	 WGH = ATAN(WGH / (1.0-WGH))
	 CS  = CS * 1.E7

         NS2 = NSAM/2
         NR2 = NROW/2
         NL2 = NSLICE/2
C
	DO K=1,NSLICE
	   IZ = K-1
	   IF (IZ .GT. NL2) IZ = IZ-NSLICE
	   DO J=1,NROW
	      IY = J-1
	      IF (IY .GT. NR2) IY = IY-NROW
	      DO I=1,LSM
	         IX = (I-1)/2

                 AK = SQRT(FLOAT(IX*IX)+FLOAT(IY*IY)+FLOAT(IZ*IZ))*SC
                 IF (AK .NE. 0.0) THEN
                    AZ=0.0
                 ELSE
                    AZ = PI/2.
                 ENDIF
C                AZ  = ATAN2(0.0,AK)
                 AZR = AZZ * (PI / 180.)
                 DZZ = DZ + DZA / 2*SIN(2*(AZ-AZR))

                 CALL TFD(TF,CS,DZZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)

                 B(IX+1) = CMPLX(TF*SIGN,0.0)
              ENDDO
              CALL WRTLIN(LUN,B,LSM,J+(K-1)*NROW)
           ENDDO
	 ENDDO

         END
