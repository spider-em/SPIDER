
C++*********************************************************************
C
C $$ LONGER.FOR -- OCT 88 ArDean Leith
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
C $ LONGER
C
C   PURPOSE:  CONCATENATES TWO STRINGS WITH ERROR CHECKING FOR LENGTHS
C
C **********************************************************************

	SUBROUTINE LONGER(STRING,NLET,SUFFIX,IRTFLG)

 
        CHARACTER *(*) STRING,SUFFIX
        CHARACTER*1 NULL
 
        NULL=CHAR(0)

        LENSUF = INDEX(SUFFIX,NULL) - 1
        IF (LENSUF .LT. 0) LENSUF = LEN(SUFFIX)

        LENSTR = LEN(STRING)
        IF ((NLET + LENSUF + 1) .GT. LENSTR) THEN
C          TOO LONG
           IRTFLG = 1

        ELSEIF (LENSUF .EQ. 0) THEN
C          SUFFIX HAS NULL LENGTH
           STRING(NLET+1:NLET+1) = NULL
           IRTFLG = 2

        ELSE
           STRING(NLET+1:NLET+LENSUF+1) = SUFFIX(1:LENSUF) // NULL
           NLET = NLET + LENSUF
           IRTFLG = 0
        ENDIF

        RETURN
        END

