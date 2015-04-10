
C++*********************************************************************
C
C SIZCHK.F -- NEW DEC 2010    AUTHOR: ARDEAN LEITH
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
C SIZCHK(UNUSED,NX1,NY1,NZ1,ITYPE1,
C               NX2,NY2,NZ2,ITYPE2,IRTFLG)
C
C    PURPOSE:  CHECK ON IO IMAGE SIZE.
C
C    PARAMETERS:
C        UNUSED                  FOR FUTURE                    (SENT)
C        NX1,NY1,NZ1             PREVIOUS DIMENSIONS           (SENT)
C        NX, NY, NZ              CURRENT  DIMENSIONS           (SENT)
C        ITYPE1                  PREVIOUS FILE TYPE            (SENT)
C        ITYPE2                  CURRENT  FILE TYPE            (SENT)
C        IRTFLG                  ERROR RETURN FLAG (0=NORMAL)  (RET)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE SIZCHK(UNUSED,NX1,NY1,NZ1,ITYPE1,
     &                           NX2,NY2,NZ2,ITYPE2,IRTFLG)

        IMPLICIT NONE

        CHARACTER (LEN=*) :: UNUSED
        INTEGER           :: NX1,NY1,NZ1,ITYPE1
        INTEGER           :: NX2,NY2,NZ2,ITYPE2,IRTFLG

        IRTFLG = 1
        IF (NX1 .NE. NX2) THEN
            CALL ERRT(102,'X DIMENSION NOT SAME AS FIRST',NX2)
            RETURN

        ELSEIF (NY1 .NE. NY2) THEN
            CALL ERRT(102,'Y DIMENSION NOT SAME AS FIRST',NY2)
            RETURN

        ELSEIF (NZ1 .NE. NZ2) THEN
C           EARLIER FILE HAS DIFFERING DIMENSIONS
            CALL ERRT(102,'Z DIMENSION NOT SAME AS FIRST',NZ2)
            RETURN

        ELSEIF (ITYPE1 .NE. ITYPE2) THEN
C           EXISTING STACK FILE FORMAT NOT SAME AS PREVIOUS FORMAT
            CALL ERRT(102,
     &             'IMAGE TYPE MUST BE SAME AS FIRST FILE',ITYPE2)
            RETURN

        ENDIF

        IRTFLG = 0

        END
