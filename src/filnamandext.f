
C++*********************************************************************
C
C  FILNAMANDEXT.F -- SIMPLIFIED FROM FILCAN.F  JAN 99 --  ArDean Leith
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
C   FILNAMANDEXT(FILIN,EXTEN,FILOUT,NLET,CALLERRT,IRTFLG)
C
C   PURPOSE: CONSTRUCTS FILE NAME FROM FILENAME & FILE EXTENSION 
C            STRINGS.
C            THIS SUBROUTINE IS SOMEWHAT OS DEPENDENT. IT IS WRITTEN 
C            FOR THE FILE-NAMING CONVENTIONS OF UNIX OS.
C
C            THE FILE NAME HAS THE FOLLOWING FORMAT: <FILEN>.<EXTEN>
C
C  PARAMETERS:  FILIN    INPUT FILE NAME              (SENT)
C               EXTEN    EXTENSION NAME               (SENT)
C               FILOUT   OUTPUT FILE NAME             (RETURNED)
C               NLET     LENGTH OF FILOUT             (RETURNED)
C               CALLERRT LOGICAL FLAG TO CALL ERRT    (SENT)
C               IRTFLG   ERROR FLAG (0 IS NORMAL)     (RETURNED)
C
C **********************************************************************

        SUBROUTINE FILNAMANDEXT(FILIN,EXTEN,FILOUT,NLET,CALLERRT,IRTFLG)

        INCLUDE 'CMBLOCK.INC'

        LOGICAL         CALLERRT
        CHARACTER *(*)  FILIN,FILOUT,EXTEN

C       SET DEFAULT IRTFLG
        IRTFLG = 1

C       FIND NUMBER OF USED CHARACTERS IN INPUT FILANEM
        NLETF = LNBLNKN(FILIN)
        IF (NLETF .LE. 0) THEN
C          NO FILIN SENT TO THIS ROUTINE
           WRITE(NOUT,*) '*** PGM ERROR, MISSING FILENAME'
           IF (CALLERRT) CALL ERRT(100,'FILNAMANDEXT',IDUM)
           RETURN
        ENDIF

C       FIND NUMBER OF USED CHARACTERS IN INPUT EXTENSION
        NLETE = LNBLNKN(EXTEN)

C       FIND MAX. NUMBER OF CHARACTERS IN RETURNED FILOUT
        LENO  = LEN(FILOUT)

C       FIND NUMBER OF CHARACTERS WANTED IN RETURNED FILOUT
        NLET = NLETF + NLETE + 1

        IF (NLET .GT. LENO .OR. NLET .LE. 1) THEN
C          OVERFLOW OR UNDERFLOW OF FILOUT

           WRITE(NOUT,90) FILIN(1:NLETF),EXTEN 
90         FORMAT(' *** ERROR, MAKING FILE NAME FROM: ',A,' & ',A)
           IF (CALLERRT) CALL ERRT(100,'FILNAMANDEXT',IDUM)
           RETURN
        ENDIF

C       PUT BASE FILENAME IN FILOUT
        FILOUT(1:NLETF) = FILIN(1:NLETF)
        NLET = NLETF

        IF (NLETE .GT. 0) THEN
C          ADD DOT, THEN ADD EXTENSION
           FILOUT(NLETF+1:) = '.'

C          ADD FILENAME EXTENSION TO FILOUT
           FILOUT(2+NLETF:) = EXTEN(1:NLETE)
           NLET = NLETF + 1 + NLETE
        ENDIF

        IRTFLG = 0

        RETURN
        END

