
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

       SUBROUTINE RYS1(X,Y,KG,MAP)

       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       DIMENSION MAP(30,42)

       I=21.-Y*7.0
       J=45.+X*15.
       IF(I.GT.42) I=42
       IF(I.LE.0) I=1
       IF(J.GT.90) J=90
       IF(J.LE.0) J=1
       K=(J-1)/3+1
       L=(J-K*3)+3
       M=MAP(K,I)
       IF(L-2) 1,2,3
 1     N=M/1024
       IF(N.EQ.0) GOTO 10
       IF(N .EQ. KG) RETURN
       N=(M-N*1024)+16384
       GOTO 100
 10    N=M+KG*1024
       GOTO 100
 2     N=M/32-(M/1024)*32
        IF(N.EQ.0) GOTO 20
       IF(N.EQ.KG) RETURN
       N=(M-N*32)+512
        GOTO 100
 20    N=M+KG*32
       GOTO 100
 3     N=M-(M/32)*32
       IF(N.EQ.0) GOTO 30
       IF(N.EQ.KG) RETURN
       N=(M-N)+16
       GOTO 100
 30    N=M+KG
 100   MAP(K,I)=N
       RETURN
       END
