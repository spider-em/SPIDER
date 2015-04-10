
C ++************************************************************************
C
C PICKSV
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
C                                                                          *
C  $$  PICKSV
C
C      PURPOSE: PICKS A X OR Y SLICE OUT OF A 3-D VOLUME FILE
C
C        0         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

	SUBROUTINE PICKSV(LUN3,LUN2,NSAM,NROW,NSLICE)

        INCLUDE 'CMBLOCK.INC'
	COMMON BUF(1)

	IF (FCHAR(4:4) .NE. 'X') THEN

          CALL RDPRMI(IROW,IDUM,NOT_USED,'Y-SLICE NUMBER')
          IF (IROW .LE. 0 .OR. IROW .GT. NROW) THEN
            CALL ERRT(24,'PICKSV',NE)
            RETURN
          ENDIF

          DO ISLICE=1,NSLICE
            IRECT = (ISLICE-1)*NROW+IROW
            CALL REDLIN(LUN3,BUF,NSAM,IRECT)
            CALL WRTLIN(LUN2,BUF,NSAM,ISLICE)
          ENDDO

        ELSE

          CALL RDPRMI(ISAM,IDUM,NOT_USED,'X-SLICE NUMBER')
          IF (ISAM .LE. 0 .OR. ISAM .GT. NSAM) THEN
            CALL ERRT(24,'PICKSV',NE)
            RETURN
          ENDIF

          DO ISLICE=1,NSLICE
            IRECT = (ISLICE-1)*NROW
            DO I=1,NROW
              CALL REDLIN(LUN3,BUF,NSAM,IRECT+I)
              BUF(NSAM+I) = BUF(ISAM)
            ENDDO
            CALL WRTLIN(LUN2,BUF(NSAM+1),NROW,ISLICE)
          ENDDO
        ENDIF

	END
