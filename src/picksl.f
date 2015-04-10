
C++*******************************************************************
C
C $$ PICKSL.FOR
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
C $$ PICKSL:   PICK SLICE FROM 3-D DENSITY FILE AND STORE IT IN 
C              A 2-D FILE.
C
C    PICKSL(LUN2,LUN3,NSAM,NROW,NSLICE)
C
C	  LUN2		        LOGICAL UNIT NUMBER 
C 	  LUN3		        LOGICAL UNIT NUMBER
C	  NSAM,NROW,NSLICE	DIMENSIONS OF 3-D FILE
C
C--*******************************************************************

      SUBROUTINE PICKSL(LUN3,LUN2,NSAM,NROW,NSLICE)

      INCLUDE 'CMBLOCK.INC'
      COMMON BUF(1)

      CALL RDPRMI(ISLICE,IDUM,NOT_USED,'SLICE NUMBER')
      IF (ISLICE .LE .0 .OR. ISLICE .GT. NSLICE) THEN
         CALL ERRT(24,'PICKSL',NE)
         RETURN
      ENDIF

      NREC1 = (ISLICE-1)*NROW+1
      NREC2 = NREC1 + NROW-1

      DO I = NREC1,NREC2
        CALL REDLIN(LUN3,BUF,NSAM,I)
        CALL WRTLIN(LUN2,BUF,NSAM,I-NREC1+1)
      ENDDO

      END

