
C++*********************************************************************
C
C  GETLBNO.F  
C                              3 DIGITS OK         DEC 2005 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  GETLBNO(STRING,LBNO,IRTFLG)
C
C  PURPOSE:      GET LABEL NUMBER FROM STRING
C        
c  PARAMETERS:   STRING  STRING CONTAINING LB NUMBER (ONLY)     (SENT)
C                LBNO    LABEL NUMBER                       (RETURNED)
C                IRTFLG  ERROR FLAG (UNUSED)                (RETURNED)
C
C--*********************************************************************

	SUBROUTINE GETLBNO(STRING,LBNO,IRTFLG)

        CHARACTER(LEN=*)  :: STRING

        CHARACTER(LEN=80) :: TEMP1
        LOGICAL           :: ISDIGI

        INTEGER           :: LBNO,IRTFLG
        INTEGER           :: NCHAR1,IGO,NCHAR


        IRTFLG = 0        ! UNUSED?
        LBNO   = -1       ! ERROR RETURN

C       FIND LAST NON-BLANK            
        NCHAR1 = lnblnk(STRING)
        
C       WHAT IF IT IS ZERO?
        IF (NCHAR1 .LE. 0) RETURN

        TEMP1 = STRING(1:NCHAR1)
        CALL SSUPCAS(TEMP1)

        IGO = INDEX(TEMP1,'LB')
        IF (IGO <= 0) RETURN

        NCHAR = 1
        IF (ISDIGI(TEMP1(IGO+3:IGO+3))) NCHAR = 2
        IF ((IGO+4) < NCHAR1 .AND. ISDIGI(TEMP1(IGO+4:IGO+4))) THEN
           NCHAR = 3
        ENDIF 

        !write(3,*)' In getlbno: igo,temp1,string:', igo, temp1,string

        READ(TEMP1(IGO+2:IGO+2+NCHAR-1),8000) LBNO
8000    FORMAT(I6)

        !write(3,*)' In getlbno: string,lbno:', temp1,string

	END



