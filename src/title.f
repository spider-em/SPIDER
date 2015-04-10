
C++*********************************************************************
C
C TITLE.F  -- CREATED NOV 9 1987     ARDEAN LEITH
C             LONG FILE NAMES FEB 89 ARDEAN LEITH
C             REWRITTEN JUNE 1999    ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C     TITLE(LUN,CTITLE,LENTIT,SAYIT,IRTFLG)
C
C     PURPOSE:   ALTER TITLE IN HEADER OF SPIDER FILE
C
C     PARAMETERS: LUN       LOGICAL UNIT NUMBER OF TITLE FILE
C                 CTITLE    NEW TITLE
C                 LENTIT    NEW TITLE LENGTH
C                 SAYIT     LOGICAL FLAG TO ECHO NEW TITLE
C                 IRTFLG    UNUSED
C
C--*******************************************************************

        SUBROUTINE TITLE(LUN,CTITLE,LENTIT,SAYIT,IRTFLG)

	INCLUDE 'CMBLOCK.INC' 

        CHARACTER *(*) CTITLE
        LOGICAL        SAYIT

C       PUT NEW TITLE IN CURRENT HEADER OBJECT
        CALL LUNSETTITLE(LUN,CTITLE,IRTFLG)

C       REPLACE CURRENT HEADER OBJECT BACK IN THE CURRENT FILE
        CALL LUNWRTCURHED(LUN,IRTFLG)

        IF (SAYIT) THEN
           LENT = MIN(LENTIT,40)
           WRITE(NOUT,90) CTITLE(1:LENT)
90         FORMAT('  NEW TITLE: ',A)

           IF (LENTIT > 40) THEN
C             WRITE NEXT 60 CHAR OF TITLE
              LENT = MIN(LENTIT-40,60)
              WRITE(NOUT,91) CTITLE(40:40+LENT)
91            FORMAT('            ',A)

              IF (LENTIT > 100) THEN
C                WRITE LAST 60 CHAR OF TITLE
                 LENT = MIN(LENTIT-100,60)
                 WRITE(NOUT,91) CTITLE(100:100+LENT)
              ENDIF
           ENDIF
        ENDIF

        END

