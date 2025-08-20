C++*********************************************************************
C
C LNBLNKN.F   -- NEW FEB 1999                       AUTHOR: ARDEAN LEITH
C                ADDED TAB TRAP NOV 00                      ARDEAN LEITH                        
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
C    LNBLNKN(STRING) 
C
C    PURPOSE:       FIND POSITION OF LAST NON-BLANK, PRINTING CHAR IN 
C                   STRING
C
C    PARAMETERS :   STRING    CHAR. STRING
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       INTEGER FUNCTION LNBLNKN(STRING)

       CHARACTER *(*) STRING

C      FIND LOCATION OF LAST NON-NULL, PRINTING CHARACTER IN STRING
       LNBLNKN = lnblnk(STRING)

C      Check for null character inside string
       LENN    = INDEX(STRING,CHAR(0)) 
       IF (LENN .GE. 1 .AND. LENN .LE. LNBLNKN) LNBLNKN = LENN - 1

C      If lnblnkn is 0, then string(lnblnkn:lnblnkn) below will be invalid
       IF (LNBLNKN .LE. 0) THEN
          LNBLNKN = 0
          RETURN
       ENDIF

C      PRINT *, 'lnblnkn.f : 54: LNBLNKN=', LNBLNKN
C      IF (STRING(LNBLNKN:LNBLNKN) .LE. CHAR(32)) THEN
       IF (ICHAR(STRING(LNBLNKN:LNBLNKN)) .LE. 32) THEN
C         LAST POSITION IS STILL NON PRINTING
          LENB = LNBLNKN
          DO IDX = LENB,1,-1
C            IF (STRING(LNBLNKN:LNBLNKN) .GT. CHAR(32)) THEN
             IF (ICHAR(STRING(idx:idx)) .GT. 32) THEN
                LNBLNKN = IDX
                RETURN
             ENDIF
          ENDDO
          LNBLNKN = 0 
       ENDIF

       RETURN
       END

