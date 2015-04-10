
C++*********************************************************************
C
C RDPRMI.F                 USED RDPRINC             FEB 99  Ardean Leith 
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
C RDPRMI(I1,I2,NOT_USED,PROMPT)
C 
C PURPOSE:        READS ONE OR TWO INTEGERS FROM INPUT
C
C    PARAMETERS:   I1,I2    NUMBERS                           (RETURNED)
C                  NOT_USED ANY MORE                              (SENT)
C                  PROMPT   SOLICITATION MESSAGE                  (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*******************************************************************

      SUBROUTINE RDPRMI(I1,I2,NOT_USED,PROMPT)

      INCLUDE      'CMBLOCK.INC' 

      CHARACTER *(*)  PROMPT

      VAL1 = 0
      VAL2 = 0

      CALL RDPRINC(PROMPT,2,.TRUE.,NOT_USED,VAL1,VAL2,VAL3,IRTFLG)
      IF (IRTFLG .EQ. -1) RETURN

      I1 = VAL1
      I2 = VAL2

      RETURN
      END

