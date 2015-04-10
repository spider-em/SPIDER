C ++********************************************************************
C                                                                      *
C  MTPR                                                               *
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
C  MTPR                                                                *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE MTPR(M0,M)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)
c	INTEGER*4 LUN50,LUN51
c       DIMENSION J(M)
c       DIMENSION R(M0,M0)
c         CHARACTER*4 VNAM(M)
      DATA L1/9/

       M=M0-1
      I1=0
      I2=0
      I3=0
       JSEC=0
 9    I1=I2+1
      I2=I1+L1
      IF(I2.GT.M) I2=M
 13   JSEC=JSEC+1
C      WRITE(LUN51,17) JSEC
C  17   FORMAT(/' SECTION',I5)
C      WRITE(LUN51,40)(VNAM(I),I=I1,I2)
C 40   FORMAT(//3X,10(7X,A4))
C      WRITE(LUN51,27)(J(I),I=I1,I2)
C 27    FORMAT(3X,10I11)
C 30   FORMAT(1X,I4,1X,A4,1X,10(1X,1PE10.3))
      DO  I=I1,M
      I3=I3+1
      IF(I3.GT.I2) I3=I2
C       WRITE(LUN51,30)J(I),VNAM(I),(R(I,K),K=I1,I3)
      ENDDO
      IF(I2.LT.M) GOTO 9
      RETURN
       END
 
