
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
C  MRECTL(LUN,NSAM,NROW,RP,IDIM)
C
C  PURPOSE:     DRAW RECTANGLE IN IMAGE
C
C  PARAMETERS:  LUN          I/O UNIT                             (SENT)
C               NSAM,NROW    DIMENSIONS                           (SENT)
C               RP           FILL VALUE                           (SENT)
C               IDIM                                              (SENT)
C **********************************************************************

	SUBROUTINE MRECTL(LUN,NSAM,NROW,RP,IDIM)

	REAL  :: BUF(NROW)
	
	CALL RDPRMI(IX,IY,NOT_USED,
     &      'COORDINATES OF UPPER LEFT POINT')

	CALL RDPRMI(IXOFF,IYOFF,NOT_USED, 'X & Y OFFSETS')
	IF (IXOFF .LT. 0 .OR. IYOFF .LT. 0) THEN
           CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',NF)
	   RETURN
        ENDIF

	IYSTRT = MAX(1,IY)
	IYEND  = MIN(NROW,IY+IYOFF)
	IF (IYSTRT.GT.NROW .OR. IYEND.LE.0)  THEN
           CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',NF)
	   RETURN
        ENDIF

	IXSTRT = MAX(1,IX)
	IXEND  = MIN(NSAM,IX+IXOFF)
	IF (IXSTRT.GT.NSAM .OR. IXEND.LE.0)  THEN
           CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',NF)
	   RETURN
        ENDIF

C       LOOP OVER ALL ROWS
	DO I=IYSTRT,IYEND
	   CALL REDLIN(LUN,BUF,NSAM,I)
	   IF  ((IDIM .EQ. 2) .OR.
     &          (I  .EQ.IYSTRT .AND. IYSTRT.EQ.IY) .OR.
     &          (I  .EQ.IYEND  .AND. IYEND.EQ.(IY+IYOFF))) THEN
              BUF(IXSTRT:IXEND) = RP
           ELSE
	      IF (IX .EQ. IXSTRT)        BUF(IX)       = RP
	      IF ((IX+IXOFF) .EQ. IXEND) BUF(IX+IXOFF) = RP
	   ENDIF

	   CALL WRTLIN(LUN,BUF,NSAM,I)
        ENDDO

	END
