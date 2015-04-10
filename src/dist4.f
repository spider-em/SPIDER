
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

       SUBROUTINE DIST4(M1,NMAX,MD,KG0,D,V,TMEAN,
     & AR,N,E,IHISTI,X)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)
       DIMENSION N(1)
       DIMENSION AR(NMAX,MD),D(M1,M1),V(M1),TMEAN(1),X(1),E(NMAX)
       DIMENSION  IHISTI(NMAX,NMAX)

       ZMIN=1.0E30
       DO 1 I=1,NMAX
       IF(N(I).EQ.0) GOTO 1
       DO  J=1,MD
       X(J)=0.
       ENDDO
       DO  J=1,M1
       Z=V(J)-TMEAN(J)
       DO  K=1,MD
       X(K)=X(K)+Z*D(J,K)
       ENDDO
       ENDDO
       Z=0.
       DO  K=1,MD
       Z=Z+(X(K)-AR(I,K))**2
       ENDDO
       Z=Z-E(I)
       IF(Z.GT.ZMIN) GOTO 1
       ZMIN=Z
       L=I
 1     CONTINUE

       IHISTI(L,KG0)=IHISTI(L,KG0)+1

       END
