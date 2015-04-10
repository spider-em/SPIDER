C++*********************************************************************
C
C REGPARSE.F                                  NEW AUGUST 00 ARDEAN LEITH
C                            GLOBAL REG. SUPPORT   NOV 2005 ARDEAN LEITH
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
C REGPARSE(STRING,IREG,IGO,IEND,IRTFLG)
C
C PURPOSE:       EXTRACTS REGISTER NUMBER FROM A STRING CONTAINING X????...
C       RETURNS: REG. NUMBER AND LOCATION OF X.... STRING (IGO..IEND)
C
C PARAMETERS:    STRING        INPUT STRING TO BE EVALUEATED        SENT
C                IGO           START OF REGISTER STRING             RET.
C                IEND          END OF REGISTER STRING               RET.
C                needx         NEEDS X AT START OF REG STRING      SENT
C                IRTFLG        ERROR RETURN FLAG (0 IS NORMAL)      RET.
C
C--*********************************************************************

        SUBROUTINE REGPARSE(STRING,IREG,IGO,IEND,NEEDX,IRTFLG)

C       EXTRACTS REGISTER NUMBER FROM A STRING CONTAINING X????...
C       RETURNS: REG. NUMBER AND LOCATION OF X.... STRING (IGO..IEND)

        CHARACTER *(*) STRING
        LOGICAL        NEEDX

        IRTFLG = 1
        
        LENT = LEN(STRING)
        IGO  = INDEX(STRING,'X')
        IF (IGO .LE. 0) IGO  = INDEX(STRING,'x')

        IF (IGO .LE. 0 .AND. NEEDX) RETURN

        IEND = VERIFY(STRING(IGO+1:),'0123456789')
        IF (IEND .EQ. 0) THEN
            IEND = LENT
        ELSE 
           IEND = IEND
        ENDIF

        READ(STRING(IGO+1:IEND),*,IOSTAT=IERR) IREG
        IF (IERR .NE. 0) RETURN

        IRTFLG = 0
        RETURN
        END
