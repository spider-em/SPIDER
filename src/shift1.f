
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
C       THIS SUBROUTINE WILL SHIFT AN ARRAY, A, CIRCULARLY BY IS.
C
C       CALL SHIFT1(A,b,N,IS)
C         A     REAL ARRAY OF DIMENSION N TO BE SHIFTED
c         b     output real array
C         N     DIMENSION
C         IS    INTEGER SPECIFYING SHIFT
C
C--*******************************************************************

      SUBROUTINE SHIFT1(MA,mo,M,IS)

      REAL MA(m),mo(m)

      ISH = MOD(IS,M)
      IF (ISH.LT.0) ISH = ISH+M
	m1=-ish-1+m
      DO  I = 1,m
	mo(i)=ma(mod(i+m1,m)+1)
      ENDDO
      END
