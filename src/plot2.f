
C++*********************************6/23/80 *** VAX 11/20/81 *** al 17/9/87
C
C PLOT2.F                                  LONG FILE NAMES JAN 25 89 al
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
C	PLOT2
C
C       COMMANDS:    TP
C
C--*******************************************************************

        SUBROUTINE PLOT2(MAXDIM)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM)   ::   FILNAM
        COMMON /COMMUN/ FILNAM

	PARAMETER (NFUNC=1)

        CHARACTER *2   FUNC(NFUNC)

	DATA FUNC/'TP'/

	LUN1   = 12
        MAXIM  = 0

 	DO 10 I=1,NFUNC
           IF (FCHAR(1:2) .NE. FUNC(I)) GOTO 10
           IFUNC=I
           GOTO 12
10	CONTINUE
        RETURN

C	FIND AND OPEN INPUT FILE
12	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &		   MAXIM,'IMAGE INPUT',.FALSE.,IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

C       3-D POSTSCRIPT PLOT
        CALL D3DPLT(LUN1,NSAM,NROW,NSLICE,MAXDIM)

9999    CLOSE(LUN1)

        RETURN
	END

