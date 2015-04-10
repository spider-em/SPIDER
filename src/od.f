
C++*********************************************************************
C
C    OD.FOR
C                   LONG FILENAMES                 JAN 89 al
C                   MAXNAM                         JUL 14 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C    OD(LUN1,LUN2,LUN3,NSAM,NROW,MAXDIM)
C        
C    OD:   OD-CONVERSION ACCORDING TO LOOKUP-TABLE.  APPEARS TO BE
C          ANCIENT AND DOES NOT USE SPIDER AUX FILE OPENING!!
C
C    CALLER:  UTIL3
C
C--*********************************************************************

	SUBROUTINE OD(LUN1,LUN2,LUN3,NSAM,NROW,MAXDIM)

	COMMON YINT(801),B(1)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FILNAM
       
	CALL FILERD(FILNAM,NLET,DATEXC,'TABLE',IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

	OPEN(UNIT=LUN3,FILE=FILNAM,STATUS='OLD')

	READ(LUN3,101) YINT
101	FORMAT(1X,10F12.6)
	CLOSE(LUN3)

	DO  I=1,NROW
          CALL REDLIN(LUN1,B,NSAM,I)
          DO  K=1,NSAM
            J = INT(200.*B(K)+1.1)
            IF (J .LT. 1) J=1
            IF (J .GT. 801) THEN
               B(NSAM+K)=YINT(801)+(YINT(801)-YINT(800))* FLOAT(I-801)
            ELSE
               B(NSAM+K)=YINT(J)
            ENDIF
	  ENDDO
	  CALL WRTLIN(LUN2,B(NSAM+1),NSAM,I)
	ENDDO

        END
