
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

       SUBROUTINE RYS3(I0,LINE,MAP,NAM)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)
       DIMENSION MAP(30,42)
        CHARACTER*1   IX,LINE(90),NAM(12)
       DATA IX /'*'/

       J=I0
       DO 1 I=1,30
       N=MAP(I,J)
       K=3*I
       M=N/1024
       IF(M.EQ.0) GOTO 3
        IF(M.GT.12) GOTO 2
       LINE(K-2)=NAM(M)
       GOTO 3
 2     LINE(K-2)=IX
 3     M=(N-M*1024)/32
       IF(M.EQ.0) GOTO 5
       IF(M.GT.12) GOTO 4
        LINE(K-1)=NAM(M)
       GOTO 5
 4     LINE(K-1)=IX
 5     M=N-(N/32)*32
       IF(M.EQ.0) GOTO 1
        IF(M.GT.12) GOTO 6
       LINE(K)=NAM(M)
       GOTO 1
  6    LINE(K)=IX
 1     CONTINUE
       END
