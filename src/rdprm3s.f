
C++*********************************************************************
C
C RDPRM3S.F               USED RDPRINC             FEB 99  ARDEAN LEITH 
C                         ADAPTED FROM RDPRM2S     AUG 99  ARDEAN LEITH
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
C RDPRM3S(F1,F2,F3,NOT_USED,PROMPT,IRTFLG)
C
C PURPOSE:        READ A TRIPLET OF FLOATING POINT NUMBERS
C                 CAN DO REGISTER SUBSTITUTION
C
C   PARAMETERS:
C        F1       FIRST NUMBER ENTERED                        (RETURNED)
C        F2       SECOND NUMBER ENTERED                       (RETURNED)
C        F3       THIRD NUMBER ENTERED                        (RETURNED)
C        NOT_USED                                                 (SENT)
C        STRING   SOLICITATION MESSAGE                            (SENT)
C        IRTFLG   RETURN FLAG (0 IS NORMAL,                   (RETURNED)
C                             -1 IS GOTO PREVIOUS QUESTION
C                             -3 IS ACCEPTED NULL RETURN
C
C NOTE: DOES NOT ALTER RECEIVED VALUES FOR F1, F2, F3 (UNLIKE SPIDER)
C
C--*******************************************************************

      SUBROUTINE RDPRM3S(F1,F2,F3,NOT_USED,PROMPT,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      CHARACTER *(*) PROMPT

      CALL RDPRINC(PROMPT,3,.FALSE.,NOT_USED,F1,F2,F3,IRTFLG)
      IF (IRTFLG .EQ. -1) RETURN

      RETURN
      END

