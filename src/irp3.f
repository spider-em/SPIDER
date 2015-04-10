
C ++********************************************************************
C                                                                      *
C IRP3                                                                 *
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
C  IRP3                                                                *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE IRP3(Q1,Q2,NSAM,NROW,NSLICE,NSAM1,NROW1,NSLICE1,LUN2)

	DIMENSION Q1(NSAM,NROW,NSLICE),Q2(NSAM1)
	DOUBLE PRECISION PX,PY,PZ,RX,RY,RZ
	DOUBLE PRECISION TMP1,TMP2,TMP3

C       REMAINING CASES
	RX = DBLE((FLOAT(NSAM-1))-0.0001)/FLOAT(NSAM1-1)
	RY = DBLE((FLOAT(NROW-1))-0.0001)/FLOAT(NROW1-1)
	RZ = DBLE((FLOAT(NSLICE-1))-0.0001)/FLOAT(NSLICE1-1)
	PZ = 1.0
	DO IZ=1,NSLICE1
           PY   = 1.0
           IOZ  = PZ
           DZ   = DMAX1(PZ-IOZ,1.0D-5)
           TMP3 = (1.0D0-DZ)
           DO IY=1,NROW1
              PX   = 1.0
              IOY  = PY
              DY   = DMAX1(PY-IOY,1.0D-5)
              TMP2 = (1.0D0-DY)
	      DO IX=1,NSAM1
	         IOX  = PX
	         DX   = DMAX1(PX-IOX,1.0D-5)

	         TMP1 = (1.0D0-DX)
	         Q2(IX)=
     &		   TMP1 * TMP2 * TMP3
     &			* Q1(IOX,IOY,IOZ)
     &		+   DX * TMP2 * TMP3
     &			* Q1(IOX+1,IOY,IOZ)
     &	  	+ TMP1 *   DY *(1.0D0-DZ)
     &			* Q1(IOX,IOY+1,IOZ)
     &		+ TMP1 * TMP2 * DZ
     &			* Q1(IOX,IOY,IOZ+1)
     &		+   DX *   DY *TMP3
     &			* Q1(IOX+1,IOY+1,IOZ)
     &		+   DX * TMP2 * DZ 
     &			* Q1(IOX+1,IOY,IOZ+1)
     &		+ TMP1 *   DY * DZ 
     &			* Q1(IOX,IOY+1,IOZ+1)
     &		+   DX *   DY * DZ 
     &			* Q1(IOX+1,IOY+1,IOZ+1)

	         PX = PX+RX
	      ENDDO
	      CALL WRTLIN(LUN2,Q2,NSAM1,IY+(IZ-1)*NROW1)
	      PY = PY + RY
	   ENDDO
	   PZ = PZ + RZ
	ENDDO

	END
