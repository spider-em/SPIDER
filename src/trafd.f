C++*********************************************************************
C
C TRAFD.F
C                                   USED OPFILE NOV 00 ARDEAN LEITH
C                                   OPFILEC     FEB  03 ARDEAN LEITH
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
C   TRAFD(LUN)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE TRAFD(LUN)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC' 
 
         CHARACTER(LEN=MAXNAM)   ::  FILNAM

         COMMON          B
         COMMON /COMMUN/ FILNAM
         CHARACTER       NULL,ANS,Y,D,E,S
         REAL            B(512),LAMBDA,KM

         DATA Y/'Y'/,D/'D'/,E/'E'/,S/'S'/,PI/3.14159/

         NULL = CHAR(0)

         CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

	CALL RDPRM(CS,NOT_USED,'SPHERICAL ABERRATION CS[MM]')
           IF (CS < 0.0001)    CS = 0.0001
         CALL RDPRM2(DZ,LAMBDA,NOT_USED,
     &      'DEFOCUS(ANGSTROEMS), LAMBDA(ANGSTROEMS)')
         CALL RDPRMI(NSAM,NDUM,NOT_USED,'NUMBER OF SP.FREQ.PTS')
         CALL RDPRM(KM,NOT_USED,'MAXIMUM SPATIAL FREQUENCY[A-1]')
         CALL RDPRM2(Q,DS,NOT_USED,
     &      'SOURCE SIZE[A-1], DEFOCUS SPREAD[A]')
         CALL RDPRM2(DZA,AZZ,NOT_USED,'ASTIGMATISM[A], AZIMUTH[DEG]')
         CALL RDPRM2(WGH,ENV,NOT_USED,
     &      'AMPL CONTRAST RATIO [0-1], GAUSSIAN ENV HALFW [FOU UNITS]')
         ENV    = 1./ENV**2

         IFORM  = 1
         NROW   = NSAM
	 NSLICE = 1
         MAXIM  = 0
         CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         SC=KM/FLOAT(NSAM/2)
         CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &       '(D)IFFRACTOGRAM / (E)NVELOPE / (S)TRAIGHT',NULL,IRTFLG)
         IE = 0
         IF (ANS.EQ.E) IE=1

	 WGH = ATAN(WGH/(1.0-WGH))
	 CS = CS*1.E7

         NS1 = (NSAM/2+1)
         DO  I=1,NROW
            DO  K=1,NSAM

               AK = SQRT(FLOAT(K-NS1)**2 + FLOAT(I-NS1)**2)*SC
               AZ = PI/2.
C               IF (K.EQ.NS1) GOTO 5
               AZ  = ATAN2(FLOAT(I-NS1),FLOAT(K-NS1)) + PI/2.
5              AZR = AZZ*(PI/180.)
               DZZ = DZ+DZA/2*SIN(2*(AZ-AZR))

               CALL TFD(B(K),CS,DZZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)
               IF (ANS .NE. S) B(K)=B(K)*B(K)
            ENDDO
            CALL WRTLIN(LUN,B,NSAM,I)
	 ENDDO
         END
