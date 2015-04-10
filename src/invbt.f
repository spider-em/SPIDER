
C ++********************************************************************
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE INVBT(X,P,Q,IFAULT,BETAIN)

      IMPLICIT  REAL*8  (A-H,O-Z)
      LOGICAL*1 INDEX
      DATA  ACU/0.1D0/

      BETAIN=X
C
C
      IFAULT=1
      IF(P.LE.0.0.OR.Q.LE.0.0)RETURN
      IFAULT=2
      IF(X.LT.0.0.OR.X.GT.1.0)RETURN
      IFAULT=0
C            IF(X.LT.ACU.OR.X.GT.(1.0-ACU))  RETURN
C
C
      PSQ=P+Q
      CX=1.0-X
      IF(P.GE.PSQ*X)GO TO 1
      XX=CX
      CX=X
      PP=Q
      QQ=P
      INDEX=.TRUE.
      GO TO 2
    1 XX=X
      PP=P
      QQ=Q
      INDEX=.FALSE.
    2 TERM=1.0
      AI=1.0
      BETAIN=1.0
      WNS=QQ+CX+PSQ
C
C
      RX=XX/CX
    3 TEMP=QQ-AI
      IF(ABS(WNS).LT.ACU)   RX=XX
    4 TERM=TERM*TEMP*RX/(PP+AI)
      BETAIN=BETAIN+TERM
      TEMP=ABS(TERM)
      IF(TEMP.LE.ACU.AND.TEMP.LE.ACU*BETAIN)   GO TO 5
      AI=AI+1.0
      WNS=WNS-1.
      IF(WNS.GE.0.)GO TO 3
      TEMP=PSQ
      PSQ=PSQ+1.0
      GO TO 4
C
C
    5 BET1=GAM(P)+GAM(Q)-GAM(P+Q)
      BET1=PP*DLOG(XX)+(QQ-1.0)*DLOG(CX)-BET1
      BET1=DSIGN(1.0D0,BET1)*DMIN1(709.0D0,DABS(BET1))
      BETAIN=BETAIN*DEXP(BET1)/PP
      IF(INDEX)   BETAIN=1.0-BETAIN
      END
