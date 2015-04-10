C++*********************************************************************
C
C TRAFC.F
C                                   USED OPFILE NOV  00 ARDEAN LEITH
C                           RECTANGULAR OUTPUTS OCT  01 BILL BAXTER
C                                   OPFILEC     FEB  03 ARDEAN LEITH
C                                   WANT_CT     MAR  04 ARDEAN LEITH
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
C  TRAFC.F
C
C CHANGED 9/5/94 TO INCLUDE COMPLETE SIGN REVERSAL AND INTELLIGIBLE
C          GAUSSIAN PARAMETER JF
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE TRAFC(LUN,WANT_CT)

	 INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC' 
 
         CHARACTER(LEN=MAXNAM)   ::  FILNAM
         COMMON /COMMUN/FILNAM

         REAL          LAMBDA,KM
         COMMON        B(1)
         COMPLEX       B
         CHARACTER     NULL
         LOGICAL       WANT_CT

         PARAMETER (QUADPI = 3.1415926535897932384626)

         NULL = CHAR(0)

         CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

	CALL RDPRM(CS,NOT_USED,'SPHERICAL ABERRATION CS[MM]')
           IF (CS < 0.0001)    CS = 0.0001

         CALL RDPRM2(DZ,LAMBDA,NOT_USED,
     &        'DEFOCUS(ANGSTROEMS), LAMBDA(ANGSTROEMS)')
 
         CALL RDPRMI(NSAM,NDUM,NOT_USED,'DIMENSIONS OF OUTPUT ARRAY')

         CALL RDPRM(KM,NOT_USED,'MAXIMUM SPATIAL FREQUENCY[A-1]')

         CALL RDPRM2(Q,DS,NOT_USED,
     &        'SOURCE SIZE[A-1], DEFOCUS SPREAD[A]')

         CALL RDPRM2(DZA,AZZ,NOT_USED,'ASTIGMATISM[A], AZIMUTH[DEG]')

         IF (WANT_CT) THEN
            CALL RDPRM(WGH,NOT_USED,'AMPL CONTRAST RATIO [0-1]')
            ENV = 0.0
         ELSE
            CALL RDPRM2(WGH,ENV,NOT_USED,
     &      'AMPL CONTRAST RATIO [0-1], GAUSSIAN ENV HALFW [FOU UNITS]')
            IF (ENV .NE. 0.0) THEN
               ENV = 1./ENV**2
            ENDIF
         ENDIF

         CALL  RDPRM(SIGN,NOT_USED,'SIGN (+1 or -1.)')

         IF (MOD(NSAM,2) .EQ. 0)  THEN
            IFORM = -12
            LSM   = NSAM+2
         ELSE
            IFORM = -11
            LSM   = NSAM+1
         ENDIF

         IXC    = NSAM/2+1
         IF (NDUM.EQ.0)  THEN
            NROW   = NSAM
            IYC    = IXC
         ELSE
            NROW   = NDUM
            IYC    = NROW/2+1
         ENDIF

         NSLICE = 1
         MAXIM  = 0
         IRTFLG = 0
         CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,LSM,NROW,NSLICE,
     &                MAXIM,' ',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        SC=KM/FLOAT(NSAM/2)
         SCX = 2.0 / NSAM
         SCY = 2.0 / NROW

         IE =0
C        IE=0 SELECTS TRANSFER FUNCTION OPTION IN SUBROUTINE TFD
         WGH = ATAN(WGH/(1.0-WGH))
         CS  = CS*1.E7

         DO  K=1,NROW
            KY = K-1
            IF (K.GT.IYC) KY = KY-NROW
            DO  I=1,LSM,2
               KX = (I-1)/2

C              Changed AK to handle rectangular images
C              AK = SQRT(FLOAT(KY)**2 + FLOAT(KX)**2)*SC
               AK = KM * SQRT((KX*SCX)**2 + (KY*SCY)**2)

C              AZ = QUADPI/2.
	       IF (KX.NE.0) THEN
                  AZ = ATAN2(FLOAT(KY),FLOAT(KX)) + QUADPI/2.
	       ELSE
	          AZ =  QUADPI/2.
               ENDIF
	       
	       
	       
               AZR = AZZ*(QUADPI/180.)
               DZZ = DZ+DZA/2*SIN(2*(AZ-AZR))

               CALL TFD(TF,CS,DZZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)
               IF (WANT_CT) THEN
                  IF (TF .GE. 0.0) THEN
                     B(KX+1) = CMPLX(1.0, 0.0) * SIGN
                  ELSE
                     B(KX+1) = CMPLX(-1.0, 0.0) * SIGN
                  ENDIF
               ELSE
                  B(KX+1) = CMPLX(TF*SIGN, 0.0)
               ENDIF
            ENDDO
            CALL WRTLIN(LUN,B,LSM,K)
         ENDDO

         END

