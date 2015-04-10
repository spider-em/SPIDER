C++*********************************************************************
C
C MIRROR.F   REFORMATTED, NX                     MAY 2012 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  MIRROR(LUN1,LUN2, NX,NY,NZ)
C
C PURPOSE: MIRROR IMAGE/VOLUME AROUND CENTER ON ONE OF XYZ AXES
C
C--*********************************************************************

	SUBROUTINE MIRROR(LUN1,LUN2,NX,NY,NZ)

        INCLUDE 'CMLIMIT.INC'
        COMMON /IOBUF/ BUF(NBUFSIZ)

	INTEGER           :: LUN1,LUN2,NX,NY,NZ

	INTEGER           :: NCHAR,K,I,J
	CHARACTER (LEN=1) :: XY
	LOGICAL           :: EVEN
	CHARACTER (LEN=1) :: NULL = CHAR(0)

	IF (NZ .LE. 1)  THEN
           CALL RDPRMC(XY,NCHAR,.TRUE.,
     &        'MIRROR AT (Y) OR (X) AXIS (DEF=Y)',NULL,IRTFLG)
	ELSE
           CALL RDPRMC(XY,NCHAR,.TRUE.,
     &        'MIRROR AT (Z), (Y) OR (X) AXIS (DEF=Y)',NULL,IRTFLG)
	ENDIF
        IF (IRTFLG .NE. 0) RETURN

	IF (XY == 'X') THEN
	   EVEN = MOD(NY,2) == 0

	   DO K=1,NZ

	      IF (EVEN)  THEN
	         DO I=1,NY
	            CALL REDLIN(LUN1,BUF,NX,I+(K-1)*NY)
	            CALL WRTLIN(LUN2,BUF,NX,
     &                       MOD(NY+1-I,NY)+1+(K-1)*NY)
	         ENDDO
	      ELSE
	         DO I=1,NY
	            CALL REDLIN(LUN1,BUF,NX,I+(K-1)*NY)
	            CALL WRTLIN(LUN2,BUF,NX,NY-I+1+(K-1)*NY)
	         ENDDO
	      ENDIF
	   ENDDO

	ELSEIF (XY == 'Z')  THEN
	   EVEN = MOD(NZ,2) == 0

	   DO K=1,NZ
	      IF (EVEN)  THEN
	         DO I=1,NY
                    CALL REDLIN(LUN1,BUF,NX,I+(K-1)*NY)
                    CALL WRTLIN(LUN2,BUF,NX,
     &                         I+MOD(NZ+1-K,NZ)*NY)
	         ENDDO
	      ELSE
	         DO I=1,NY
	            CALL REDLIN(LUN1,BUF,NX,I+(K-1)*NY)
	            CALL WRTLIN(LUN2,BUF,NX,I+(NZ-K)*NY)
	         ENDDO
	      ENDIF
	   ENDDO

	ELSE
	   EVEN = MOD(NX,2) == 0

	   DO I=1,NY*NZ
	      CALL REDLIN(LUN1,BUF,NX,I)

	      IF (EVEN)  THEN
	         DO J=1,NX
	           BUF(NX+MOD(NX+1-J,NX)+1) = BUF(J) ! 61-->1,120-->2
	         ENDDO
	      ELSE
	         DO J=1,NX
	            BUF(NX+NX+1-J) = BUF(J)          ! 120-->1,119-->2
	         ENDDO
	      ENDIF

	      CALL WRTLIN(LUN2,BUF(NX+1),NX,I)
	   ENDDO
	ENDIF

	END
