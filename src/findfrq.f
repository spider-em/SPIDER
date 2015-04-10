
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

	SUBROUTINE  FINDFRQ(NX,NY,NZ,NDIM,NROW,NSLICE,IX,IY,IZ)

	IF(NX.EQ.1) THEN
	IX=0
		IF(NY.EQ.1)  THEN
		IY=0
			IF(NZ.EQ.1)  THEN
			IZ=0
			RETURN
			ELSEIF(NZ.EQ.2)  THEN
			IZ=NSLICE/2
			RETURN
			ELSE
			IZ=(NZ-1)/2
			ENDIF
		ELSEIF(NY.EQ.2)  THEN
		IY=NROW/2
			IF(NZ.EQ.1)  THEN
			IZ=0
			RETURN
			ELSEIF(NZ.EQ.2)  THEN
			IZ=NSLICE/2
			RETURN
			ELSE
			IZ=(NZ-1)/2
			ENDIF
		ELSE
		IY=(NY-1)/2
		IZ=NZ-1
		IF(IZ.GT.NSLICE/2)  IZ=IZ-NSLICE
		RETURN
		ENDIF
	ELSEIF(NX.EQ.2)  THEN
	IX=NDIM/2
		IF(NY.EQ.1)  THEN
		IY=0
			IF(NZ.EQ.1)  THEN
			IZ=0
			ELSEIF(NZ.EQ.2)  THEN
			IZ=NSLICE/2
			ELSE
			IZ=(NZ-1)/2
			ENDIF
			RETURN
		ELSEIF(NY.EQ.2)  THEN
		IY=NROW/2
			IF(NZ.EQ.1)  THEN
			IZ=0
			ELSEIF(NZ.EQ.2)  THEN
			IZ=NSLICE/2
			ELSE
			IZ=(NZ-1)/2
			ENDIF
			RETURN
		ELSE
		IY=(NY-1)/2
		IZ=NZ-1
		IF(IZ.GT.NSLICE/2)  IZ=IZ-NSLICE
		RETURN
		ENDIF
	ELSE
	IX=(NX-1)/2
	IY=NY-1
	IF(IY.GT.NROW/2)  IY=IY-NROW
	IZ=NZ-1
	IF(IZ.GT.NSLICE/2)  IZ=IZ-NSLICE
	ENDIF
	END
