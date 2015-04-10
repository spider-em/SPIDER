

C++*********************************************************************
C
C SSUPCAS.F  -- CREATED AUG 86 
C
C **********************************************************************
C *  AUTHOR:  ArDean Leith                                                 *
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
C    SSUPCAS(STRING)
C
C    PURPOSE:  CONVERTS LOWER CASE STRINGS TO UPPER CASE
C
C    PARAMETERS:             STRING     STRING TO BE CONVERTED
C
C    CALLED BY:              SSINP0
C 
C--*********************************************************************

      SUBROUTINE SSUPCAS(STRING)

 
      CHARACTER *(*) STRING

      DATA IOFF/-32/

      ILEN = LEN(STRING)
      
      DO  I=1,ILEN
        IF (STRING(I:I) .GE. 'a' .AND. STRING(I:I) .LE. 'z') THEN
          STRING(I:I) = CHAR(ICHAR(STRING(I:I)) + IOFF)
        ENDIF
      ENDDO

      RETURN
      END

