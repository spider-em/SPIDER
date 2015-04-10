
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

       SUBROUTINE MAHAL(M0,MGR0,W,A,G,XMEAN,N,DIS)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)
	INTEGER*4  M1,N2
       DIMENSION W(M0,M0),A(M0,M0),G(MGR0,MGR0),DIS(MGR0,MGR0)
       DIMENSION XMEAN(M0,1),N(1)

         MGR=MGR0-1
       M1=M0-1
       NSUM=0
       DO  I=1,MGR
       IF(N(I).LT.1) GO TO 6
       NSUM=NSUM+N(I)-1
 6     DO  J=1,MGR
       DIS(I,J)=0.
       G(I,J)=0.
       ENDDO
       ENDDO
       X=1./REAL(NSUM)
       DO  K=1,M1
       DO  L=1,M1
       A(K,L)=W(K,L)*X
       ENDDO
       ENDDO
       CALL MTNV(A,M0,DET)
       DO 1 I=1,MGR
       EN1=N(I)
       IF(EN1.LE.0.5) GOTO 1
       DO 4 J=I,MGR
       IF(J.EQ.I) GOTO 4
       EN2=N(J)
       IF(EN2.LE.0.5) GOTO 4
       X=0.0
       DO  K=1,M1
       Y=XMEAN(K,I)-XMEAN(K,J)
       DO  L=K,M1
       Z=Y*A(K,L)*(XMEAN(L,I)-XMEAN(L,J))
       IF(L.GT.K) Z=Z+Z
       X=X+Z
       ENDDO
       ENDDO
       DIS(I,J)=X
       DIS(J,I)=X
       X=X*EN1*EN2/(EN1+EN2)
       N2=NSUM-M1+1
       X=X*REAL(N2)/(REAL(NSUM)*REAL(M1))
       G(J,I)=ALPHAINT(X,M1,N2)
 4     CONTINUE
 1     CONTINUE
       END
