C ++********************************************************************
C
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
C
C PURPOSE: TO SET CERTAIN POINTS IN A FILE
C
C	MPOINT(LUN,NSAM,NROW,RP)
C
C	LUN :       LOGICAL UNIT NUMBER
C	NSAM,NROW : FILE DIMENSIONS
C	RP:         VALUE TO BE SUBSTITUTED OR USED
C
C **********************************************************************

	SUBROUTINE MPOINT(LUN,NSAM,NROW,RP)

	COMMON BUF(1)

	CALL RDPRMI(NPT,NDUM,NOT_USED,'NUMBER OF POINTS')
	IF (NPT .LT. 1) RETURN

	DO I=1,NPT
	   CALL RDPRMI(IX,IY,NOT_USED,
     &             'COORDINATES FOR EACH POINT')
	   IF (IX.LE.0 .OR. IX.GT.NSAM) GOTO 100
	   IF (IY.LE.0 .OR. IY.GT.NROW) GOTO 100

	   CALL REDLIN(LUN,BUF,NSAM,IY)
	   BUF(IX) = RP
	   CALL WRTLIN(LUN,BUF,NSAM,IY)
100	CONTINUE
	ENDDO
	END
