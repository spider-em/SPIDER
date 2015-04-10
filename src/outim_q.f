C++*********************************************************************
C
C OUTIM_Q.F             USED OPFILE                          NOV 99 AL
C                       USED WRITEV                          SEP 00 BR
C                       LEN=MAXNAM                           JUL 14 AL
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014 Health Research Inc.,                         *
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
C  OUTIM_Q.F
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE  OUTIM_Q(X,LSD,NSAM,NROW,IT)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC' 
 
         CHARACTER(LEN=MAXNAM) :: FINPAT,FINPIC,DOCFIL,OUTIMA,OUTDOC
         COMMON /FISPEC/
     &     FINPAT,FINPIC,DOCFIL,OUTIMA,OUTDOC,NLET,NLETI,NLIMA,NLDOC

         DIMENSION     X(LSD,NROW)

         DATA INPIC/69/

         CALL  FILGET(OUTIMA,FINPIC,NLIMA,IT,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         MAXIM  = 0
         NSLICE = 1
         IFORM  = 1
         CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'U',IFORM,NSAM,NROW,NSLICE,
     &               MAXIM,' ',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CALL WRITEV(INPIC,X,LSD,NROW,NSAM,NROW,NSLICE)
 
         CLOSE(INPIC)
         END
