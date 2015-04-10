
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

	SUBROUTINE  PREPSL_2(NSAM,NSLICE,NN,NMAT,IPCUBE,RI)

	INTEGER  IPCUBE(5,*)
	LOGICAL  FIRST
	COMMON /PAR/  LDPX,LDPY,LDPZ,LDPNMX,LDPNMY

C
C IPCUBE: 1 - beginning
C         2 - end
C         3 - ix
C         4 - iy
C         5 - iz
C                   
	R=RI*RI
	NN=0
	NMAT=0
C
	DO    I1=1,NSLICE
	T=I1-LDPZ
	XX=T*T
	FIRST=.TRUE.
	DO    I3=1,NSAM
	NMAT=NMAT+1
	T=I3-LDPX
        RC=T*T+XX
	IF(FIRST)  THEN
	IF(RC-R)  80  ,80,14
80	FIRST=.FALSE.
	NN=NN+1
	IPCUBE(1,NN)=NMAT
	IPCUBE(2,NN)=NMAT
	IPCUBE(3,NN)=I3
	IPCUBE(4,NN)=1
	IPCUBE(5,NN)=I1	
	ELSE
	IF(RC.le.R)  IPCUBE(2,NN)=NMAT
	ENDIF
14	CONTINUE
	ENDDO
16	CONTINUE
20	CONTINUE
	ENDDO
	END
