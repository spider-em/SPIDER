
C++*********************************************************************
C
C RDPRM.F                 USED RDPRINC             FEB 99  Ardean Leith 
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
C    RDPRM(REALNO,NOT_USED,PROMPT)
C
C    PURPOSE:      READ A FLOATING POINT NUMBER
C
C    PARAMETERS:   F1       REAL NUMBER                       (RETURNED)
C                  NOT_USED                                       (SENT)
C                  PROMPT   SOLICITATION MESSAGE                  (SENT)
C
C--*******************************************************************

      SUBROUTINE RDPRM(F1,NOT_USED,PROMPT)

      INCLUDE         'CMBLOCK.INC' 
       
      CHARACTER *(*) PROMPT

      F1 = 0.0

      CALL RDPRINC(PROMPT,1,.FALSE.,NOT_USED,F1,F2,F3,IRTFLG)
      IF (IRTFLG .EQ. -1) RETURN

      RETURN
      END

