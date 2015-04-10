
C++*********************************************************************
C
C ARASQ.F       ADAPTED FROM ARAS.FOR FOR CHAR. AUG 89  al
C               CHANGED ERROR MESG.             SEP 97  al
C               USED FILNAMSUB                  DEC 99  ARDEAN LEITH
C               REMOVED FILNAMSUB               APR 01  ARDEAN LEITH
C               [] REG. SUPPORT                 NOV 05  ARDEAN LEITH
C               GLOBAL                          JUN 09  ARDEAN LEITH
C               ERROR MSG. FORMAT               NOV 09  ARDEAN LEITH
C               SKIP GLOBAL BUG                 NOV 09  ARDEAN LEITH
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
C ARASQ(STRING,NLET,GLOBAL,IRTFLG)
C
C PURPOSE:       SUBROUTINE TO EVALUATE ARITHMETIC EXPRESSION AND SET
C                REGISTERS FROM OPERATION LINE
C
C PARAMETERS:    STRING        INPUT STRING TO BE EVALUEATED        SENT
C                NLET          LENGTH OF INPUT STRING               SENT
C                GLOBAL        VARIABLE IS TO BE GLOBAL             SENT
C                IRTFLG        ERROR RETURN FLAG (0 IS NORMAL)      RET.
C
C--*********************************************************************

	SUBROUTINE ARASQ(OLDSTR,NLET,MAKEGLOBAL,IRTFLG)

	COMMON/UNITS/LUN,NIN,NOUT

        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*)      :: OLDSTR
        CHARACTER(LEN=160)    :: NEWSTR
        LOGICAL               :: MAKEGLOBAL,ISGLOBAL

	IRTFLG = 1

C       REMOVE ALL BLANKS
	CALL SHRINKQ(OLDSTR,NLET,NEWSTR,NLET2)
	IF (NLET2 .LE. 0) GOTO 9999
        
        NEQ  = INDEX(NEWSTR(:NLET2),'=')
        IF (NEQ .EQ. 0) THEN
C          NO EQUALS SIGN FOUND, EVALUATE 
           CALL EXPRQ(NEWSTR(:NLET2),NLET2,VALUE,IRTFLG)

           WRITE(NOUT,14) VALUE
14         FORMAT('  ',1pG15.8)

           RETURN
        ENDIF

C       NEQ NOW POINTS TO POSITION OF EQUAL SIGN

C       SKIP ANY GLOBAL FLAG
        IGO = VERIFY (NEWSTR(1:NEQ-1),'GLOBALglobalCc')
        IGO = MAX(1,IGO)

C       EVALUATE LEFT HAND NEQ-1 CHARACTERS FOR REGISTER NUMBER
	IF (NEWSTR(IGO:IGO) .NE. '[' ) GOTO 9999

        IF (MAKEGLOBAL) THEN
           CALL REG_FIND_IREG('GLO',NEWSTR(IGO:NEQ-1),
     &                        ISGLOBAL,IREG,IRTFLG)
!       ELSEIF(NEWSTR(1:3) .EQ. 'LOC' .OR. 
!     &         NEWSTR(1:3) .EQ. 'loc') THEN
!           CALL REG_FIND_IREG('LOC',NEWSTR(IGO:NEQ-1),
!     &                        ISGLOBAL,IREG,IRTFLG)
        ELSE
           CALL REG_FIND_IREG('LOC',NEWSTR(IGO:NEQ-1),
     &                        ISGLOBAL,IREG,IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 9999

C       NOW EVALUATE RIGHT HAND NLET2-NEQ CHARS
	CALL EXPRQ(NEWSTR(NEQ+1:NLET2),NLET2-NEQ,VALUE,IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 9999

C       NOW SET REGISTER FROM LEFT HAND NEQ-1 CHARACTERS
        CALL REG_SET_BYNUM(IREG,VALUE,ISGLOBAL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

	RETURN

9999	WRITE(NOUT,*)'*** ERROR: ARITHMETIC ASSIGNMENT FAILED'
        RETURN

	END

