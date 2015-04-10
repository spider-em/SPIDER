
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

	SUBROUTINE  BCKC0(CUBE,LTC,DM,B,NSAM,IPCUBE,NN)

        DIMENSION  DM(9)
	DIMENSION  CUBE(LTC),B(NSAM)
	INTEGER  IPCUBE(5,NN)
	COMMON /PAR/  LDPX,LDPY,LDPZ,LDPNMX,LDPNMY

c$omp parallel do  private(i,xb,xbb,j,iqx,dipx)
	DO    I=1,NN
	   XB=(IPCUBE(3,I)-LDPX)*DM(1)+(IPCUBE(5,I)-LDPZ)*DM(3)
	   DO    J=IPCUBE(1,I),IPCUBE(2,I)
	      XBB=(J-IPCUBE(1,I))*DM(1)+XB
	      IQX=IFIX(XBB+FLOAT(LDPNMX))
	      DIPX=XBB+LDPNMX-IQX
	      CUBE(J)=CUBE(J)+B(IQX)+DIPX*(B(IQX+1)-B(IQX))
	   ENDDO
	ENDDO
C     1                 +(1.0-DIPX)*B(IQX)
C     2                 +     DIPX *B(IQX+1)
	END
