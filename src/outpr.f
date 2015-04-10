C++*********************************************************************
C
C  OUTPR.F              USED SAVDN1                NOV 99 ARDEAN LEITH
C                       LEN=MAXNAM                 JUL 14 ARDEAN LEITH
C
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
C OUTPR(PARA,NIMA,IT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE OUTPR(PARA,NIMA,IT)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FINPAT,FINPIC,DOCFIL,OUTIMA,OUTDOC
        COMMON  /FISPEC/
     &     FINPAT,FINPIC,DOCFIL,OUTIMA,OUTDOC,NLET,NLETI,NLIMA,NLDOC

        DIMENSION     PARA(3,NIMA),DLIST(4)

        DATA          NDOC/55/

        CALL  FILGET(OUTDOC,FINPIC,NLDOC,IT,INTFLAG)

        NRUN = 0
C       REWIND AND OVERWRITE EXISTING DOC FILE
        IAP  = 0

        DO  I=1,NIMA
            DLIST(1) = I
            DLIST(2) = PARA(1,I)
            DLIST(3) = PARA(2,I)
            DLIST(4) = PARA(3,I)

            CALL SAVDN1(NDOC,FINPIC,DLIST,4,NRUN,IAP)
            NRUN     = 1
        ENDDO

        CLOSE(NDOC)

        END
