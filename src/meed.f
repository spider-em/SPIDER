
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

	SUBROUTINE MEED(Q,A,NSAM,L,Y,K)

	DIMENSION Q(NSAM,NSAM),A(K),Y(NSAM,NSAM)

	LH=L/2
	K21=K/2+1
	DO    J=1,NSAM
	DO    I=1,NSAM
	LB=0
	DO    MJ=-LH,LH
	MJM=MOD(J+MJ+NSAM-1,NSAM)+1
	DO    MI=-LH,LH
	LB=LB+1
	A(LB)=Q(MOD(I+MI+NSAM-1,NSAM)+1,MJM)
	ENDDO
	ENDDO
	CALL   FSORT(A,K)
	Y(I,J)=A(K21)
	ENDDO
	ENDDO
	END
