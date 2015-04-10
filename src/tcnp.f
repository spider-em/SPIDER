
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

	SUBROUTINE TCNP(A,B,POL,N,C,WORK)

	IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
	DIMENSION POL(2),C(2),WORK(4)

	IF(N-1)2,1,3
1 	POL(1)=C(1)
2 	RETURN
3 	POL(1)=C(1)+C(2)*B
	POL(2)=C(2)*A
	IF(N-2)2,2,4
4 	WORK(1)=1.
	WORK(2)=B
	WORK(3)=0.
	WORK(4)=A
	XD=A+A
	X0=B+B
	DO  J=3,N
	P=0.
	DO  K=2,J
	H=P-WORK(2*K-3)+X0*WORK(2*K-2)
	P=WORK(2*K-2)
	WORK(2*K-2)=H
	WORK(2*K-3)=P
	POL(K-1)=POL(K-1)+H*C(J)
 	P=XD*P
	ENDDO
	WORK(2*J-1)=0.
	WORK(2*J)=P
 	POL(J)=C(J)*P
	ENDDO
	END
