C++*************************************************** 12/11/79 1/11/81 VAX
C
C EVALNQ.F                    ADAPTED FROM EVALN.FOR FOR CHAR AUG 89 AL
C                             REWRITTEN FOR POLISH EVAN CHEN  MAR 98 EC
C                             REWRITTEN                       MAY 98 AL
C                             IRTFLG RETURNED IF NRPN=0       NOV 09 AL
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
C EVALNQ(EXPR,NLET,RETVALUE,IRTFLG)
C
C PURPOSE:   EVALUATES EXPRESSIONS, RETURNS VALUE
C
C       EXPR       CHARACTER STRING CONTAINING EXPRESSION     (SENT)
C       NLET       LENGTH OF EXPR                             (SENT)
C       RETVALUE   RETURNED VALUE OF EXPRESSION               (RETURNED)
C       IRTFLG     ERROR FLAG                                 (RETURNED)
C 
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE EVALNQ(ISTACK,EXPR,NLET,RETVALUE,IRTFLG)

        CHARACTER(LEN=*)      :: EXPR

        INTEGER,PARAMETER     :: IVALEN  = 40  ! RPN LENGTH LIMIT
        INTEGER,PARAMETER     :: IRPNLEN = 80  ! RPN LENGTH LIMIT
        INTEGER,PARAMETER     :: LENEXP  = 80  ! EXPR LENGTH LIMIT

        INTEGER               :: IRPN(IRPNLEN)
        REAL                  :: VAL(IVALEN)

C       CONVERT EXPR TO REVERSE POLISH NOTATION ARRAY 
	CALL POLISH(ISTACK,EXPR(1:NLET),NLET,
     &               IRPN,NRPN,VAL,NVAL,IRTFLG)

        IF (IRTFLG .NE. 0 .OR. NRPN .EQ. 0) THEN
C          EMPTY OR UNDECIPHERABLE EXPRESSION
           IRTFLG = 1
           RETURN
        ENDIF

C       EVALUATE RPN EXPRESSION, PIXVAL IS NOT USED HERE
	CALL CALC(IRPN,NRPN,VAL,PIXVAL,RETVALUE,IRTFLG)

	END
