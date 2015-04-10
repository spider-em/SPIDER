
C++*********************************************************************
C
C FILERD.F   ADAPTED FROM FILRD AND SPFILE         OCT 88 ARDEAN LEITH
C            REWRITTEN FOR NEW SUBSTITUTION        JUN 97 ARDEAN LEITH
C            ADDS EXTEN NOW IF SENT                FEB 99 ARDEAN LEITH
C            CAN USE PROMPT FOR INPUT --           AUG 99 ARDEAN LEITH
C            REMOVED AVGX11 & AVG0I RULES --       SEP 99 ARDEAN LEITH
C            TRAILING ~ ON PROMPT OMITS "FILE"     NOV 00 ARDEAN LEITH
C            $ EFFECTIVE NOW IF ONLY 1 CHAR.       JAN 01 ARDEAN LEITH
C            ALTERED COMMENT HANDLING              JAN 01 ARDEAN LEITH
C            ALTERED COMMENT HANDLING              MAR 02 ARDEAN LEITH
C            NLOG                                  NOV 03 ARDEAN LEITH
C            REMOVED IRTFLG INPUT                  APR 04 ARDEAN LEITH
C            TEMPLATTING CAPABILITY                DEC 04 ARDEAN LEITH
C            RDPR PARAMETERS                     04/14/05 ARDEAN LEITH
C            RDPR IRTFLG FOR LEGACY REGS         06/24/09 ARDEAN LEITH
C            PROMPT(LENP-2:LENP) .EQ. '~9~'      12/07/10 ARDEAN LEITH
C            ACCEPT EXTENSION BUG                  JUN 11 ARDEAN LEITH 
C            ! COMMENT DELIMITER                   DEC 11 ARDEAN LEITH
C            NECHO                                 SEP 12 ARDEAN LEITH
C
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
C  FILERD(FILN,NLET,EXTEN,PROMPT,IRTFLG)
C
C  PARAMETERS:  FILN    CHARACTER ARRAY FOR FILE NAME           (RET.)
C               NLET    NO. OF CHARACTERS IN FILE NAME          (RET.)
C               EXTENT  OPTIONAL FILE NAME EXTENSION            (SENT)
C                       AN INPUT EXTENSION SUPERCEDES THIS
C                       PARAMETER
C               PROMPT  FILE NAME SOLICITATION PROMPT           (SENT)
C                       A ~ (TILDE) IN FIRST CHAR. SAYS USE
C                       PROMPT INSTEAD OF INPUT
C                       A ~ (TILDE) IN LAST CHAR. SAYS SKIP
C                         "FILE" AT END OF PROMPT
C                       A ~9  IN NEXT TO LAST OR 
C                          NEXT-TO-NEXT-TO LAST
C                          ACCEPTS AN EXTENSION 
C                          (OTHERWISE DISCARDED!)
C               IRTFLG  ERROR FLAG, ZERO IS NORMAL          (SENT/RET.)
C
C  PURPOSE:     READS IN FILE NAME AND TITLE.  CAN
C               SUBSTITUTE REGISTERS AND DO LOOP INDICES INTO
C               CERTAIN FILENAME
C               
C  DESCRIPTION: PRINTS OUT <PROMPT> MESSAGE:,
C               READS IN FILE NAME FROM INPUT DEVICE. 
C               SUBSTITUTES FOR REGISTER AND LOOP INDEX. 
C               EXAMPLES IF (X0 = 1, X11=11, X21=1111, INDEX:A=44)
C
C                   avg{***A }  ---> avg044
C                   avg{**A }   ---> avg44
C                   avg{***x0}  ---> avg001
C                   avg{****x0} ---> avg0001
C                   avg{****}   ---> avg****
c                   avg{***x21} ---> ERROR (DAMAGES INVARIANT PART OF 
C                                    FILE NAME)
C
C               IF THE FIRST CHARACTER IS A '$', FILNAM FROM PREVIOUS 
C               CALL REMAINS UNCHANGED.  
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE FILERD(FILN,NLET,EXTENT,PROMPT,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*)      :: FILN,EXTENT,PROMPT
        INTEGER               :: NLET,IRTFLG

        CHARACTER(LEN=MAXNAM) :: EXTEN,COMMENT
        CHARACTER(LEN=160)    :: FILNAM,RESPONSE
        CHARACTER(LEN=100)    :: PROMPTP
        CHARACTER(LEN=1)      :: JCHAR
        LOGICAL               :: NOX,EXTENOK,GOTEX
        LOGICAL               :: ISDIGI,ISCHAR,GETANS,STRIP,LEGACYREGS
        LOGICAL               :: UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI
 
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        INTEGER               :: ICOMM,MYPID,MPIERR,LENP

        CALL SET_MPI(ICOMM,MYPID,MPIERR)

        LEGACYREGS = (IRTFLG == -999)  ! CONVERT x**
  
        LENP = LNBLNKN(PROMPT)

C       '~9' IN LAST 3 CHAR. OF PROMPT INDICATES THAT EXTENSION CAN 
C       BE READ IN ON NAME (USED TO BE IRTFLG == 9)
        EXTENOK = .FALSE.
        IF (PROMPT(LENP-1:LENP-1) == '~' .OR.
     &      PROMPT(LENP-2:LENP-2) == '~') THEN
           EXTENOK = .TRUE.
           LENP    = LENP - 2
        ENDIF
      
C       SET ERROR RETURN
        IRTFLG = 1

        IF (PROMPT(:1) == '~') THEN
C          DO NOT SOLICIT NAME, INSTEAD USE PROMPT FOR: FILN 
           FILNAM   = PROMPT(2:LENP)
           NLET     = LENP - 1
           IRTFLGT  = 0

        ELSE
C          READ THE FILE NAME FROM THE INPUT, DO NOT UPPERCASE IT!
C          THIS WILL FAIL IF PROMPT > 96 CHAR!!!        
           LENU = MIN(96,LENP)
           IF (PROMPT(LENU:LENU) == '~') THEN
C             DO NOT ADD "FILE" TO PROMPT
              LENU    = LENU - 1
              PROMPTP = PROMPT(1:LENU) 
           ELSE
              PROMPTP = PROMPT(1:LENU) // ' FILE' 
              LENU    = LENU + 5
           ENDIF

           GETANS    = .TRUE.
           UPPER     = .FALSE.
           WANTSUB   = .TRUE.
           SAYPRMT   = .NOT. SILENT
           SAYANS    = .FALSE.
           ENDATSEMI = .TRUE.
           STRIP     = .TRUE.
           IRTFLGT   = 0
           IF (.NOT. LEGACYREGS) IRTFLGT   = -999  ! DO NOT CONVERT x**
           CALL RDPR(PROMPTP(1:LENU),NLET,FILNAM,GETANS,
     &             UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLGT)
        ENDIF

        IF (NLET <= 0 .OR. IRTFLGT .NE. 0) THEN
C          NO FILE NAME IN INPUT LINE
           IRTFLG = -1
           RETURN
        ENDIF

C       REPLACE ANY TAB CHARACTERS WITH BLANKS
        ITAB = INDEX(FILNAM,CHAR(9))
        DO WHILE (ITAB > 0)
           FILNAM(ITAB:ITAB) = ' '
           ITAB = INDEX(FILNAM,CHAR(9))
        ENDDO
   
C       REMOVE ANY INITIAL BLANKS
        FILNAM = ADJUSTL(FILNAM)

C       SAVE ANY TRAILING COMMENT FROM THE INPUT LINE
        ISEMI = SCAN(FILNAM,';!')

        IF (ISEMI > 1) THEN
            NLET    = LNBLNKN(FILNAM(1:ISEMI-1)) 
            FILNAM  = FILNAM(1:NLET)
            COMMENT = FILNAM(NLET+1:)
            NLETC   = LNBLNKN(COMMENT)
        ELSE
            NLETC   = 0
            NLET    = LNBLNKN(FILNAM)
        ENDIF

        IF (FILNAM(:1) == '?') THEN
C          ?PROMPT? [var] LINE RETURNED, READ RESPONSE FROM CALLER
C          SUBSTITUTE FOR SYMBOLS & REG. FROM CALLING PROCEDURE
           CALL FRSYMPAR(FILNAM(1:NLET),FILNAM,NLET,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

	IF (FILNAM(1:NLET) == '$') THEN
C          DO NOT ALTER FILN  IF INPUT IS "$"
	   CALL ECHONAME(FILNAM,NLET,COMMENT,NLETC,MYPID)
           IRTFLG = 0
           RETURN 
   
        ELSEIF (FILNAM(1:1) == '^' .OR. FILNAM(1:1) == '*') THEN
C          RETURN IF  INPUT IS "^" OR "*"
           FILN(1:) = FILNAM(1:1)
	   CALL ECHONAME(FILN,NLET,COMMENT,NLETC,MYPID)
           IRTFLG = -1
           RETURN

        ENDIF

C       REMOVE ANY EXTENSION FROM THE FILENAME (UNLESS IT IS WANTED)
        J      = NLET
        GOTEX  = .FALSE.
        DO WHILE (J >= 1 )
           JCHAR = FILNAM(J:J)

           IF (JCHAR == ' ') THEN
C             DISCARD THIS TRAILING BLANK (or embedded blank!)
              NLET = J - 1

           ELSEIF (JCHAR == ']' .OR. 
     &             JCHAR == ':'  .OR. 
     &             JCHAR == '/' ) THEN
C             BEGINNING OF DEVICE OR DIRECTORY, HALT EXTENSION SEARCH
              EXIT

           ELSEIF (JCHAR == '.') THEN
              IF (EXTENOK) THEN
C                WANT TO USE INPUTTED FILE EXTENSION
                 IEXT  = J + 1    ! MAY HAVE > 1 EXTENSION
                 GOTEX = .TRUE.
              ENDIF
           ENDIF
           J = J - 1
        ENDDO

        IF (GOTEX) THEN
C           RECEIVED EXTENSION ON INPUT LINE
            EXTEN = FILNAM(IEXT:NLET) // NULL
            NLET  = IEXT - 2  
        ENDIF

        IF (NLET <= 0) THEN 
           CALL ERRT(101,'*** ABNORMAL FILENAME',NE)
           IRTFLG = 1
           RETURN
        ENDIF

        ILEN = LEN(FILN)
        IF (NLET >= ILEN) THEN
           WRITE(NOUT,*) 'FILENAME: ',FILNAM(1:NLET)
           CALL ERRT(102,'FILENAME TOO LONG, MUST BE < ',ILEN)
           IRTFLG = 1
           RETURN
        ENDIF
        FILN = FILNAM(1:NLET) // NULL

C       USE SENT EXTENSION IF NOT NULL AND NO EXTENSION INPUT AND OK
        LENET = LNBLNKN(EXTENT)
        IF (.NOT. GOTEX) EXTEN = NULL
        IF (LENET > 0 .AND. .NOT. GOTEX) EXTEN = EXTENT

        LENE = LNBLNKN(EXTEN)    ! LENGTH OF EXTENSION TO BE APPENDED 
        
        !WRITE(6,*) ' exten1:',lene,':',exten(1:4)

        IF (LENE <= 0) THEN
C          DO NOT PUT ANY EXTENSION AT END OF FILE NAME, ECHO NAME
           CALL ECHONAME(FILN,NLET,COMMENT,NLETC,MYPID)

        ELSE
C          PUT EXTENSION (SENT OR INPUT) AT END OF FILE NAME
           IF (LENE + NLET > LEN(FILN)) THEN
              NLET = NLET + LENE + 1
              CALL ERRT(102,'FILENAME TOO LONG',NLET)
              IRTFLG = 1
              RETURN
           ENDIF

           IF (EXTEN(1:LENE) == DATEXC(1:3)) THEN
C             USES STANDARD SPIDER EXTENSION, DO NOT ECHO IT
              CALL ECHONAME(FILN,NLET,COMMENT,NLETC,MYPID)
           ENDIF

           FILN = FILN(1:NLET) // '.' // EXTEN(1:LENE) // NULL
           NLET = NLET + LENE + 1
           IF (EXTEN(1:LENE) /= DATEXC(1:3)) THEN
C             USES NON-STANDARD SPIDER EXTENSION, ECHO IT
              CALL ECHONAME(FILN,NLET,COMMENT,NLETC,MYPID)
           ENDIF
        ENDIF

C       SET NORMAL ERROR RETURN
        IRTFLG = 0

        END



	SUBROUTINE ECHONAME(FILN,NLET,COMMENT,NLETC,MYPID)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'

        CHARACTER(LEN=*)      :: FILN
        INTEGER               :: NLET
        CHARACTER(LEN=*)      :: COMMENT
        INTEGER               :: NLETC,MYPID
 
        IF (NLETC > 0 .AND. MYPID <= 0) THEN

          IF (.NOT. SILENT) WRITE(NOUT,93) FILN(1:NLET),COMMENT(1:NLETC)
          IF (NLOG .NE. 0) THEN
             WRITE(NLOG,93) FILN(1:NLET),COMMENT(1:NLETC)
             NECHO = NECHO + 1
         ENDIF

93        FORMAT(' ',A,' ',A)

        ELSEIF (MYPID <= 0) THEN
          IF (.NOT. SILENT) WRITE(NOUT,93) FILN(1:NLET)
          IF (NLOG .NE. 0) THEN
            WRITE(NLOG,93) FILN(1:NLET)
            NECHO = NECHO + 1
         ENDIF
        ENDIF

        END
