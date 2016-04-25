
C++*********************************************************************
C
C EXPRQ.F           ADAPTED FROM EXPR.FOR FOR CHAR.  AUG 89 ArDean Leith
C                   INDEXTOREG CHANGED TO SUBROUTINE MAR 01 ArDean Leith
C                   REG_GET                          JAN 06 ArDean Leith
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
C    EXPRQ(OLDSTR,NCHAR,VALUE,IRTFLG)
C
C    PURPOSE:    SUBROUTINE TO EVALUATE GENERAL ARITHMETIC EXPRESSION
C
C    PARAMETERS:  OLDSTR     INPUT STRING                   (SENT)
C                 NCHAR      LENGTH OF INPUT STRING         (SENT)
C                 VALUE      VALUE OF EXPRESSION            (RETURNED)
C                 IRTFLG     ERROR FLAG                     (RETURNED)
C
C **********************************************************************

        SUBROUTINE EXPRQ(OLDSTR,NCHAR,VALUE,IRTFLG)

        CHARACTER(LEN=*)    :: OLDSTR
        CHARACTER(LEN=80)   :: NEWSTR
        CHARACTER(LEN=4)    :: CREG
        LOGICAL             :: ISCHAR

C       SET DEFAULT RETURN VALUES
        IRTFLG = 0
        VALUE = 0.0

C       REMOVE ALL BLANKS FROM OLDSTR
        CALL SHRINKQ(OLDSTR(1:NCHAR),NCHAR,NEWSTR,L2)

        IF (NCHAR .EQ. 1 .AND. ISCHAR(NEWSTR(1:1))) THEN

C          EXPRESSION HAS LENGTH 1, AND IT'S NOT A NUMBER, INTERPRET
C          AS A OLD-STYLE DO-LOOP INDEX - EVALUATE IT!
           CALL REG_GET(0,0,CXREG,VALUE,LOOPREG,IRTFLG)
CC         CALL REG_GET_VAR(0,NEWSTR(1:1),.TRUE.,VALUE,LOOPREG,IENDVAR,IERR)
           IF (VALUE .EQ. 0) IRTFLG = -1

        ELSEIF ((NEWSTR(:1) .EQ. 'P' .OR. NEWSTR(:1) .EQ. 'p')
     &         .AND. NCHAR .LE. 3) THEN

C          EXPRESSION STARTS WITH 'P', AND HAS LENGTH 2 OR 3, IT IS A
C          OBSOLETE SYMBOLIC REFERENCE (P1)!
           IRTFLG = -1

        ELSE
C          IF NONE OF ABOVE APPLIES, HAVE TRUE ARITHMETIC EXPRESSION!
           ILEVEL = 0
           CALL EVALNQ(ILEVEL,NEWSTR,NCHAR,VALUE,IRTFLG)
        ENDIF

        RETURN
        END
