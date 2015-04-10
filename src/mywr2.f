
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

        SUBROUTINE MYWR2(MGR0,M0,N0,W,T,A)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)
	INTEGER*4  N1,N2
c	INTEGER*4 LUN50,LUN51
c       DIMENSION JV(1)
       DIMENSION W(M0,M0),A(M0,M0),T(M0,M0)
c         CHARACTER*4  NV(1)
       N1=MGR0-1
       N2=N0-MGR0
       M1=M0-1
C        WRITE(LUN51,10) N1,N2
C 10     FORMAT(//42H VARIABLE   AMONG MEAN SQ   WITHIN MEAN SQ  ,
C     & '    F-RATIO     ETA SQUARE',7X,'SIGN'/41X,I3,
C     & ' DF1',I4,' DF2'/)
       DO  J=1,M1
       DO  I=J,M1
       X=T(J,I)-W(J,I)
       A(J,I)=X
       A(I,J)=X
       ENDDO
       ENDDO
       DO  J=1,M1
       X=A(J,J)
       ETASQ=X/T(J,J)
       AMS =X/REAL(N1)
       WMS=W(J,J)/REAL(N2)
       F=AMS/WMS
       ALPH=ALPHAINT(F,N1,N2)
C 1     WRITE(LUN51,11) JV(J),NV(J),AMS,WMS,F,ETASQ,ALPH
C  11   FORMAT(1X,I4,1X,A4,3X,G11.4,5X,G11.4,F11.2,8X,F7.4,5X,F9.3)
	ENDDO	
        END
