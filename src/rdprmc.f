
C++*********************************************************************
C
C    RDPRMC.F -- CREATED FROM RDPRMB.F          08/15/89 ArDean Leith
C                REWRITTEN                       03/4/98 ArDean Leith
C                LENGTHENED ANST                11/17/00 ArDean Leith
C                SKIP BLANK LINE ECHO VERBOSE    03/7/02 ArDean Leith
C                IQ P SPECIAL LINES             06/26/02 ArDean Leith
C                NLOG                           11/26/03 ArDean Leith
C                NO INITIAL BLANKS ON ECHO      03/30/05 ArDean Leith
C                .OPERATION.... BUG             05/25/05 ArDean Leith
C                RDPR PARAMETERS                04/14/05 ArDean Leith
C                TO NOUT ALWAYS                 02/21/06 ArDean Leith
C                IF (FCHAR(1:2) == 'FR'         07/31/06 ArDean Leith
C                DEBRAKREG                      09/05/06 ArDean Leith
C                IF (FCHAR(1:4) == 'FR N'       08/28/09 ArDean Leith
C                NECHO                          09/20/12 ArDean Leith
C                KLUDGE TO ACCEPT <CR> IN DISP  10/20/25 ArDean Leith 
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email:                                                             *
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
C  RDPRMC(ANS,NCHAR,STRIP,PROMPT,CDUM,IRTFLG)
C
C  PURPOSE: READ AN ALPHANUMERIC STRING, CHECK FOR ANY SPECIAL CHARS
C           IN AN OPERATION,
C           RETURN STRING MINUS ANY LEADING OR TRAILING BLANKS, AND
C           NUMBER OF CHARACTERS IN STRING AND A ERROR FLAG.  
C           NORMALLY CONVERTS INPUT TO UPPER CASE.  
C           STRIPS OFF ANY TRAILING SPIDER COMMENT
C
C  PARAMETERS:  ANS       ANSWER                                   RET.
C               NCHAR     NUMBER OF CHARACTERS IN THE ANSWER       RET.
C               STRIP     LOGICAL FLAG TO STRIP BLANKS FROM ANS   SENT
C               PROMPT    SOLICITATION MESSAGE                    SENT
C               CDUM      (UNUSED)
C               IRTFLG    RETURN FLAG                          SENT/RET
C                           (0 IS NORMAL, -1 IS GOTO PREVIOUS
C                           QUESTION, 1 IS END-OF_FILE) 
C                           IRTFLG: -999 -654321 ON INPUT
C                           DOES NOT CONVERT INPUT TO UPPERCASE
C
C  CALLED BY:   VERY MANY SPIDER ROUTINES
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE RDPRMC(ANS,NCHAR,STRIP,PROMPT,CDUM,IRTFLG)

        INCLUDE 'CMBLOCK.INC' 

        CHARACTER(LEN=*)     :: ANS,PROMPT,CDUM
        CHARACTER(LEN=161)   :: ANST
        LOGICAL              :: STRIP,GETANS
        LOGICAL              :: UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       SEE IF WANT TO KEEP LOWER CASE INPUT UNALTERED
        UPPER     = (IRTFLG .NE. -999)

C       KLUDGE TO ACCEPT <CR> IN SUBROUTINE DISP BUT NOT TO UPPERCASE 
        IF (IRTFLG .EQ. 654321) THEN
           UPPER = .FALSE.  ! BUT RETAIN IRTFLG FOR RDPR USE

        ELSE IF (IRTFLG .EQ. -999) THEN
           UPPER  = .FALSE.
           IRTFLG = 0

        ELSE
           UPPER  = .TRUE.
           IRTFLG = 0
        ENDIF

        GETANS    = .TRUE.

C       WANT TO SUBSTITUTE FOR VARIABLES NOW?
        WANTSUB   = .TRUE. 
        IF ((FCHAR(1:2) == 'FR' .AND. NALPH <= 2) .OR.
     &      (FCHAR(1:4) == 'FR T')  .OR. 
     &      (FCHAR(1:4) == 'FR N')  .OR.
     &      (FCHAR(1:2) == 'AR')) WANTSUB = .FALSE.
 
        SAYPRMT   = .NOT. SILENT
        SAYANS    = .FALSE.
        ENDATSEMI = .TRUE.

C       MOVE BLANKS TO THE ANSWER STRING, NECESSARY FOR SOME SPIDER CODE
        LENA = LEN(ANS)
        ANS(1:LENA) = ' '
    
C       PRINT PROMPT, READ ANSWER STRING, SKIP ANY INPUT WHICH HAS
C       COMMENT IN FIRST COL. AND READ ANOTHER INPUT LINE
        CALL RDPR(PROMPT,NCHAR,ANST,GETANS,
     &            UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NCHAR <= 0) THEN
C          IF NULL ANSWER STRING, MUST RETURN ZERO LENGTH ANSWER
           IF (MYPID <= 0) THEN
              WRITE(NOUT,*) ' '
              IF (NLOG .NE. 0) THEN
                 WRITE(NLOG,*) ' '
                 NECHO = NECHO + 1
              ENDIF
           ENDIF
           RETURN
        ENDIF

C       SET RETURNED ANSWER, TRUNCATE TO FIT LENGTH OF ANS IN CALL
        IF (NCHAR > LENA) NCHAR       = LENA
        IF (NCHAR > 0)   ANS(1:NCHAR) = ANST(1:NCHAR)

        NLET = NCHAR
        IF (.NOT. SILENT .AND. MYPID <= 0) THEN
           IF (NCHAR > 0) THEN
C             CONVERT [_x**] BACK TO X** FOR ECHO
              CALL DEBRAKXREG(ANST,NLET)

              IF (COPT == 'I') THEN 
                 WRITE(NOUT,92) ANST(1:NLET)
 92              FORMAT(5X,A)
              ELSE
                 WRITE(NOUT,90) ANST(1:NLET)
 90              FORMAT('  ',A)
              ENDIF
           ENDIF
        ENDIF
        IF (NLOG .NE. 0 .AND. MYPID <= 0 .AND. NLET > 0) THEN 
            WRITE(NLOG,*) ANST(1:NLET)
            NECHO = NECHO + 1
         ENDIF

        IF (ANS(1:1) == '^' .AND. NCHAR == 1) IRTFLG = -1

        END
