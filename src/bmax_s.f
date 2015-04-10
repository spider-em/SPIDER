
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

	SUBROUTINE  BMAX_S(BCKE,NMAT,IPCUBE,NN,BMAX)

	PARAMETER  (NMPR=64)
	DIMENSION    BCKE(NMAT),IPCUBE(5,NN)
	DIMENSION   A(NMPR),LN(NMPR),LT(NMPR)

	NT=NN/NMPR
	IF (NT.EQ.0.OR.NMPR.EQ.1)  THEN
	   CALL  BMAX_C(BCKE,NMAT,IPCUBE,NN,BMAX)
	ELSE
	   DO   I=1,NMPR
	      LT(I)=(I-1)*NT+1
	      LN(I)=AMIN0(I*NT,NN)-LT(I)+1
	   ENDDO
c$omp parallel do  private(i),schedule(static)
	   DO   I=1,NMPR
	      CALL  BMAX_C(BCKE,NMAT,IPCUBE(1,LT(I)),LN(I),A(I))
	   ENDDO
	   BMAX=A(1)
	   DO   I=2,NMPR
	      BMAX=AMAX1(BMAX,A(I))
	   ENDDO
	ENDIF
	END
