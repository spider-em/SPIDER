
C ++********************************************************************
C                                                                      *
C DISTEST                                                             *
C                                                                      *
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
C                                                                      *
C  DISTEST                                                             *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE DISTEST(M0,NSUM0,MGR0,MGR1,FA1S0,ALL,N)

       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       DIMENSION ALL(M0,M0,1),N(1)

       M1=M0-1
       MGR=MGR1
       ES=NSUM0-MGR0
       FS=FA1S0
       EK1=MGR0-1
       EM=M1
       N1=NSUM0/MGR0-1

C       WRITE(LUN51,100) MGR0,N1,N1,MGR0
C 100   FORMAT(//' VARIABLE     F-HARTLEY     M BARTLETT-BOX     SIGN',
C     & 9X,'C-COCHRAN'/11X,I2,' DF1',I4,' DF2',27X,I4,
C     & ' DF AND',I4,' GROUPS'/)

       DO  I=1,M1
          XMIN=1.E30
          XMAX=-1.E30
          SUM=0.
          XSUM=0.
          XLOG=0.
          DO 3 J=1,MGR
             K=N(J)
             IF(K.EQ.0) GOTO 3
             IF(K.EQ.1) K=2
             X=ALL(I,I,J)
             XSUM=XSUM+X
             X=X/(FLOAT(K-1))
             IF(X.GT.XMAX) XMAX=X
             IF(X.LT.XMIN) XMIN=X
             SUM=SUM+X
             XLOG=DBLE(K-1)*DLOG(X)+XLOG
 3        CONTINUE

          XMM=ES*DLOG(XSUM/ES)-XLOG
          A1=(FS-1.0/ES)/(3.0*EK1)
          F2=(EK1+2.0)/(A1*A1)
          B1=F2/(1.0-A1+(2.0/F2))
          F=(F2*XMM)/(EK1*(B1-XMM))
          N1=EK1
          IF(F2.GT.32000.0) F2=32000.0
          N2=F2

C         ALPH=ALPHAINT(F,N1,N2)
          COCH=XMAX/SUM
          HAR=XMAX/XMIN
       ENDDO

       END
