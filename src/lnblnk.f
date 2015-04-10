 
C++*********************************************************************
C
C LNBLNK.FOR -- CREATED NOV 93 AL
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
C  LNBLNK(STRING)
C
C  PURPOSE: RETURNS LENGTH OF STRING TO LAST NON-BLANK, NON-NULL CHARACTER
C
C  PARAMETERS:    STRING    TEST STRING
C
C  CALLED BY:     
C    
C  CALLS:
C
C  NOTE: USE lnblnk function on unix machines instead
C
C        1         2         3         4         5         6         7       
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

#if defined (SP_IBMSP3) || defined (__ia64)
        INTEGER FUNCTION LNBLNK(STRING)

	INCLUDE 'CMBLOCK.INC'

        CHARACTER *(*) STRING
        CHARACTER      NULL

        NULL = CHAR(0)

        IEND = LEN(STRING)
        DO I = IEND,1,-1
           IF (STRING(I:I) .NE. ' ' .AND. STRING(I:I) .NE. NULL .AND.
     &         STRING(I:I) .NE. CHAR(9)) THEN
C             FOUND A NON-BLANK CHARACTER IN STRING
              LNBLNK = I
              RETURN
           ENDIF
        ENDDO

        LNBLNK = 0
        RETURN

        END
 
#else 
 
C   THIS ROUTINE IS A SYSTEMS CALL IN UNIX AND IS NOT NEEDED THERE
C   RENAME THE OBJECT FILE TO AVOID COLLISION
 
       SUBROUTINE sp_lnblnk
 
       COMMON /UNITS/LUNC,NIN,NOUT
 
       WRITE(NOUT,*) 'NOTIFY PROGRAMMER OF BAD LNBLNK CALL.)'
       RETURN
       END
 

#endif
 
