
C++*********************************************************************
C
C SHRINKQ.F                 ADAPTED FROM SHRINK.FOR FOR CHAR. AUG 89 al
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
C    SHRINKQ:      SUBROUTINE TO SHRINK STRING BY IGNORING BLANKS AND
C                  TAB CHARACTERS
C
C    PARAMETERS:   INSTR      INPUT STRING TO BE SHRANK       (INPUT)
C                  LENIN      LENGTH OF INPUT STRING.  KEPT FOR
C                             COMPATIBILIY WITH OLD CALLS     (INPUR)
C                  OUTSTR     OUPUT SHRUNKEN STRING           (RETURNED)
C                  LENOUT     LENGTH OF SHRUNKEN STRING       (RETURNED)
C
C **********************************************************************

	SUBROUTINE SHRINKQ(INSTR,LENIN,OUTSTR,LENOUT)

        CHARACTER(LEN=*) :: INSTR,OUTSTR
        CHARACTER(LEN=1) :: CTEMP

        LENMAX = LEN(OUTSTR)
        LENS   = LENIN
        IF (LENS .EQ. 0) LENS = LEN(INSTR)

        LENOUT = 0
	DO  I=1,LENS
          CTEMP = INSTR(I:I)
          IF ((CTEMP .NE. ' ' .AND. CTEMP .NE. CHAR(9)) .AND. 
     &         LENOUT .LT. LENMAX) THEN
             LENOUT                = LENOUT + 1
             OUTSTR(LENOUT:LENOUT) = CTEMP
          ENDIF
	ENDDO

        IF (LENOUT. LT. LENMAX) THEN
C          PUT BLANKS AT END OF OUTSTR
           OUTSTR(LENOUT+1:LENMAX) = ' '
        ENDIF

	RETURN
	END
