
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


        SUBROUTINE MYWR0(LUN51,K,M,J,DET0,JV,NV,XMEAN,ALL,XMIN,XMAX,W)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)
	INTEGER*4 LUN51
       DIMENSION ALL(M,M,1),W(M,M)
        DIMENSION XMEAN(M,1),XMIN(M,1),XMAX(M,1)
       DIMENSION JV(1)
         CHARACTER*4  NV(1)

       M1=M-1
       X=J-1
       IF(J.EQ.1) X=1.0
       DO  I=1,M1
       DO  L=I,M1
       X0=ALL(I,L,K)/X
       W(I,L)=X0
       W(L,I)=X0
       ENDDO
       ENDDO
       WRITE(LUN51,10)
 10   FORMAT(//' VARIABLE',9X,'MEAN',9X,'ST.DEV.',9X,'MIN',9X,'MAX'/)
       DO  I=1,M1
       X0=SQRT(W(I,I))
       WRITE(LUN51,11)JV(I),NV(I),XMEAN(I,K),X0,XMIN(I,K),XMAX(I,K)
       ENDDO
 11    FORMAT(1X,I4,1X,A4,3X,G13.4,4X,G12.4,2G12.4)
       CALL MTNV(W,M,DET)
       DET0=DET
       WRITE(LUN51,12) DET
 12    FORMAT(//' DISPERSION DETERMINANT = ' ,G10.4)
       END
