C ++********************************************************************
C                                                                      *
C MTNV                                                                *
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
C MTNV                                                               *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE MTNV(A,M0,DET0)

       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       DIMENSION A(M0,M0)

       M=M0-1
       DET=1.0
       DO 1 J=1,M
       PVT=A(J,J)
       DET=DET*PVT
       A(J,J)=1.0
       DO  K=1,M
          A(J,K)=A(J,K)/PVT
       ENDDO
       DO 1 K=1,M
         IF(K.EQ.J) GOTO 1
       T=A(K,J)
       A(K,J)=0.0
       DO  L=1,M
          A(K,L)=A(K,L)-A(J,L)*T
       ENDDO
 1     CONTINUE
       DET0=DET
       END
