
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

        SUBROUTINE CLSS(M1,MD,NMAX,TMEAN,D,AR,E,JV
     &  ,GR,VART,MVAR,V2,L)

        DOUBLE PRECISION TMEAN(9),D(9,2),AR(3,2),
     &   E(3),V(10),VV(10),X(2),Z,ZMIN
        INTEGER*2 M1,MD,NMAX,JV(9)
        DIMENSION VART(MVAR)

C
        V(1)=GR
C       V(2)=VARI
C       V(3)=SKEW
C       V(4)=AKURT
C       V(5)=ENTP
C       V(6)=AVAV
C       V(7)=AVVR
C       V(8)=SDAV
C       V(9)=SDVR

        DO LKS =1,MVAR
        V(LKS+1)=VART(LKS)
        END DO
C
        V(10)=V2
C
        DO KL=1,M1
        JK=JV(KL)
        VV(KL)=V(JK)
        ENDDO
        ZMIN=1.0E30
       DO 1 II=1,NMAX
       DO  J=1,MD
       X(J)=0.
       ENDDO
       DO  J=1,M1
       Z=VV(J)-TMEAN(J)
       DO  K=1,MD
       X(K)=X(K)+Z*D(J,K)
       ENDDO
       ENDDO
       Z=0.
       DO  K=1,MD
       Z=Z+(X(K)-AR(II,K))**2
       ENDDO
       Z=Z-E(II)
       IF(Z.GT.ZMIN) GOTO 1
       ZMIN=Z      
       L=II
1	CONTINUE
	V(1)=FLOAT(L)
        END

