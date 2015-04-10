
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

        DOUBLE PRECISION FUNCTION CHISQ(X,N)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)

       CHISQ=1.0
       IF(X.LE.0.0.OR.N.LE.0) RETURN
       CHISQ=0.0
       X1=REAL(N)+3.85*SQRT(REAL(N))+17.587
       IF(X.GT.X1) RETURN
       X1=X/2.0
       IF(X1.GT.70.0) X1=70.0
       Y=EXP(-X1)
       SUM=0.
       XIL=1.0
       IF(2*(N/2).EQ.N) GOTO 10
       Y=Y*.797885
       T=1.0/(1.0+0.33267*X)
       T=T*(0.937298*T-.1201676)*T+0.4361836
       IF(N.EQ.1) GOTO 2
       M=(N-1)/2
       XIL=1.0/SQRT(X)
       DO  I=1,M
       XIL=XIL*X/(REAL(2*I-1))
       SUM=SUM+XIL
       ENDDO
 2     CHISQ=Y*(T+SUM)
       RETURN
 10    CONTINUE
       IF(N.EQ.2) GOTO 12
       M=(N-2)/2
       DO  I=1,M
       XIL=XIL*X*0.5/REAL(I)
       SUM=SUM+XIL
       ENDDO
  12   CHISQ=Y*(1.0+SUM)
       END
