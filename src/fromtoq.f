
C++*********************************************************************
C
C  FROMTOQ.F               ADAPTED FROM FROMTO.FOR FOR CHAR. AUG 89 al
C                   CHARINSIDE PARAMETERS CHANGED          JAN 2001 AL
C                   () = []                                SEP 2002 AL
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
C FROMTOQ(STRING,NCHAR,ILIST,RLIST,NLIST,NPARG)
C
C PARAMETERS:    STRING     INPUT STRING          (MAY BE ALTERED)
C                NCHAR      CHARS. IN STRING      (MAY BE ALTERED)
C                ILIST      INPUT LIST            (RETURNED)
C                NLIST      NO. IN LISTS          (RETURNED)
C                NPARG      MAX. LIST LENGTH      (SENT)
C
C--*********************************************************************

       SUBROUTINE FROMTOQ(STRING,NCHAR,ILIST,NLIST,NPARG)

        INCLUDE 'CMBLOCK.INC'

        CHARACTER *(*) STRING

        INTEGER        ILIST(NPARG)

	NLIST = 0

C       NS, NE ARE START AND END OF SUBSTRING CONTAINING REGISTERS
        CALL CHARINSIDE(STRING(1:NCHAR),'[',']',.TRUE.,.FALSE.,
     &                  NS,NE,NLEN)
 
	IF (NS .LE. 1 .OR. NE .LE. NS) THEN
C          NS, NE ARE START AND END OF SUBSTRING CONTAINING REGISTERS
           CALL CHARINSIDE(STRING(1:NCHAR),'(',')',.TRUE.,.FALSE.,
     &                     NS,NE,NLEN)
	   IF (NS .LE. 1 .OR. NE .LE. NS) RETURN
        ENDIF

C       RETRIEVE REGISTER LIST FROM [] OR () BRACKETED STRING
        CALL CHKSTR(STRING(NS:NE),NLEN,'IT',ILIST,DUM,NPARG,
     &              NLIST,IRTFLG)

C       CHECK FOR ERRORS
	IF (IRTFLG .NE. 0) THEN
	   WRITE(NOUT,*) '*** INCORRECT ARGUMENTS SENT TO PROCEDURE: ',
     &                   STRING(1:NCHAR)

	   CALL ERRT(16,'FROMTOQ',NE)
	ENDIF

C       ALTER NCHAR TO STOP FURTHER PROCESSING OF BRACKETED SUBSTRING
        IF (NS .GT. 2) NCHAR = LNBLNKN(STRING(1:NS-2))

        RETURN
        END

