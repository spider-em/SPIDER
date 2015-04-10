
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

       SUBROUTINE DIST1(NMAX,MD,XMEAN,N,LIN,E)
      
       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       DIMENSION XMEAN(NMAX,MD),N(NMAX),LIN(90)
       DIMENSION E(NMAX)

       STP=1./15.
       X=-3.0
       Y=3.0
       DO  I=1,90
       L=0
       XMAX=1.E30
       DO 4 J=1,NMAX
       IF(N(J).EQ.0) GOTO 4
       Z=(X-XMEAN(J,1))**2+(Y-XMEAN(J,2))**2-E(J)
       IF(Z.GE.XMAX) GOTO 4
       XMAX=Z
       L=J
 4     CONTINUE
       X=X+STP
       LIN(I)=L
       ENDDO
       RETURN
       END
