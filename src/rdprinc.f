
C++************************************************************************ 
C
C   RDPRINC.F             REMOVED FROM RDPRI       FEB 99  ArDean Leith
C                         LENGTHENED ANSW        11/17/00  ArDean Leith
C                         EXPRESSQ SHOULD BE SUB.  MAY 01  ArDean Leith
C                         DID NOT PRINT IF X       JAN 02  ArDean Leith
C                         ~PROMPT                  FEB 03  ArDean Leith
C                         NLOG                   11/26/03  ArDean Leith
C                         RDPR PARAMETERS        04/14/05  ArDean Leith
C                         ?..? LEVELS            11/28/05  ArDean Leith
C                        [] REWRITE              12/02/05  ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C   RDPRINC(PROMPT,NVAL,INTS,NOT_USED,VAL1,VAL2,VAL3,IRTFLG)
C
C   PURPOSE:    EVALUATE INPUT FOR RDPR** SUBROUTINES
C
C   PARAMETERS:
C        PROMPT      PROMPT                                       (SENT)
C        NVAL        NUMBER OF VALUES TO RETURN                   (SENT)
C        INTS        LOGICAL FLAG FOR INTEGER RETURN              (SENT)
C        NOT_USED                                                 (SENT)
C        VAL1..      VALUES (ALWAYS AS FLOATS!)               (RETURNED)
C        IRTFLG      RETURN FLAG (0 IS NORMAL,1 IS ERROR      (RETURNED)
C                             -1 IS GOTO PREVIOUS QUESTION
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE RDPRINC(PROMPT,NVAL,INTS,NOT_USED,
     &                   VAL1,VAL2,VAL3,IRTFLG)

      CHARACTER(LEN=*)  :: PROMPT
      LOGICAL           :: INTS
      REAL              :: VALUES(3)

      VALUES(1) = VAL1
      VALUES(2) = VAL2
      VALUES(3) = VAL3

      IRTFLG = 654321            ! FLAG TO ACCEPT <CR> or *

      CALL RDPRA(PROMPT,NVAL,0,INTS,VALUES,NGOT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (NGOT > 0) THEN
C        COPY THE RETURNED VALUES 
         IF (NGOT .GE. 1) VAL1 = VALUES(1)
         IF (NGOT .GE. 2) VAL2 = VALUES(2)
         IF (NGOT .GE. 3) VAL3 = VALUES(3)
      ENDIF
      END

