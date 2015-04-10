
C++*********************************************************************
C
C $$ ISDIGI.FOR -- CREATED AUG 89 al
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
C  $$ ISDIGI(TCHAR)
C
C     PURPOSE:   TO DETERMINE IF A CHARACTER IS A DIGIT (0...9) OR NOT
C
C--*********************************************************************

      LOGICAL FUNCTION ISDIGI(TCHAR)

      CHARACTER TCHAR  

      IF (TCHAR .GE. '0' .AND. TCHAR .LE. '9') THEN
        ISDIGI = .TRUE.
      ELSE
        ISDIGI = .FALSE.
      ENDIF

      END
