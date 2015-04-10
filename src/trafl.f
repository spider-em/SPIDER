
C++*********************************************************************
C
C TRAFL.F
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
C   TRAFL
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

           SUBROUTINE TRAFL

           INCLUDE 'CMBLOCK.INC'
           INCLUDE 'CMLIMIT.INC'

           COMMON        B

           CHARACTER     NULL,ANS
           REAL          LAMBDA,KM

           PARAMETER     (NDLI=3)
           DIMENSION     DLIST(NDLI)
          
           LUN9 = 70
           NULL = CHAR(0)

	CALL RDPRM(CS,NOT_USED,'SPHERICAL ABERRATION CS[MM]')
           IF (CS < 0.0001)    CS = 0.0001

           CALL RDPRM2(DZ,LAMBDA,NOT_USED,
     &              'DEFOCUS[A], LAMBDA[A]')

           CALL RDPRMI(NSAM,NDUM,NOT_USED,
     &              'NUMBER OF SPATIAL FREQUENCY PTS')

           CALL RDPRM(KM,NOT_USED,
     &              'MAXIMUM SPATIAL FREQUENCY[A-1]')

           CALL RDPRM2(Q,DS,NOT_USED,
     &              'SOURCE SIZE[A-1], DEFOCUS SPREAD[A]')

           CALL RDPRM2(WGH,ENV,NOT_USED,
     &   'AMPL CONTRAST RATIO [0-1], GAUSSIAN ENV. HALFW. [FOU. UNITS]')
           ENV    = 1. / ENV**2

           SC = KM / FLOAT(NSAM / 2)
           CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &       '(D)IFFRACTOGRAM / (E)NVELOPE / (S)TRAIGHT',NULL,IRTFLG)

           IE = 0
           IF (ANS .EQ. 'E') IE = 1
           WGH = ATAN(WGH/(1.0-WGH))
           CS = CS*1.E7

           NS1 = (NSAM/2+1)
           DO K=NS1,NSAM+1

              AK = (NS1-K)*SC
              CALL TFD(B,CS,DZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)
              IF (ANS .NE. 'S') B = B * B
              DLIST(1) = K - NS1 + 1
              DLIST(2) = B
              DLIST(3) = REAL(K-NS1) / NSAM
              CALL SAVD(LUN9,DLIST,NDLI,IRTFLG)
           ENDDO

           CLOSE(LUN9)
           CALL  SAVDC

           END


