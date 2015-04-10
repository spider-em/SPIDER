
C++*********************************************************************
C
C MYGETENV
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
C   MYGETENV(INPUT,OUTPUT,MESSAGE,IRTFLG)
C
C   PURPOSE:   RETRIEVES AN ENVIRONMENTAL VARIABLE
C
C   PARAMETERS:   INPUT     NAME OF ENVIRONMENTAL VARIABLE   (SENT)
C                 OUTPUT    VALUE OF ENVIRONMENTAL VARIABLE  (RETURNED)
C                 NCHAR     NUMBER OF CHARACTERS IN OUTPUT   (RETURNED)
C                 MESSAGE   ERROR MESSAGE                    (SENT)
C                 ITRFLG    ERROR FLAG   (0 IS NORMAL)       (RETURNED)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE MYGETENV(INPUT,OUTPUT,NCHAR,MESSAGE,IRTFLG)

	INCLUDE 'CMBLOCK.INC' 

        CHARACTER *(*) INPUT, OUTPUT, MESSAGE

        IRTFLG = 0

C       READ IN THE ENVIRONMENT VARIABLE 'INPUT'.

#if defined (SP_IFC) || defined(__INTEL_COMPILER) 
        CALL GETENVQQ(INPUT,OUTPUT)
#else
        CALL GETENV(INPUT,OUTPUT)
#endif

        NCHAR = lnblnk(OUTPUT)

        IF (NCHAR <= 0 .OR. NCHAR > 74) THEN
           IF (MESSAGE .NE. CHAR(0)) THEN
             WRITE(NOUT,*) '*** UNDEFINED ENVIRONMENTAL VARIABLE:',INPUT
             WRITE(NOUT,*) 'PUT DEFINITION IN YOUR STARTUP FILE.  E.G.'
             WRITE(NOUT,*) 'FOR C SHELL, ADD FOLLOWING TO .cshrc FILE '
             WRITE(NOUT,*) 'setenv ',INPUT,' ',MESSAGE
           ENDIF
           IRTFLG = 1
        ENDIF

        RETURN
        END

