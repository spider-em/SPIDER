
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

       SUBROUTINE   COREL(T,A,B,M,M1,NSUM)

       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       DOUBLE PRECISION  R,DF,TT,BETAI
       DIMENSION  T(M,M),A(M,M),B(M,M)
       PARAMETER (TINY=1.D-20)

       DO    I=1,M1
          X=SQRT(T(I,I))
          DO    J=1,I
             A(I,J)=T(I,J)/X/SQRT(T(J,J))
          ENDDO
       ENDDO
       DF=NSUM-2
       DO  2  I=1,M1
          DO  2  J=1,I
              R=A(I,J)*A(I,J)
             IF(R.GT.0.999)  GOTO  3
             TT=R*DSQRT(DF/(((1.-R)+TINY)*((1.+R)+TINY)))
             B(I,J)=BETAI(0.5*DF,0.5D0,DF/(DF+TT**2))
             GOTO  2
 3           B(I,J)=0.0
 2     CONTINUE
       END
