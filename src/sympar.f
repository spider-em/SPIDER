
C++*********************************************************************
C
C    SYMPAR    REWRITTEN FROM SETSYMPAR           JUN 2002 ARDEAN LEITH
C              VARIABLES NOT PARAMETERS           OCT 2005 ARDEAN LEITH
C              SYMPAR TEXT ADDED                  OCT 2006 ARDEAN LEITH
C              'FR GS' & 'FR LS'                  JAN 2009 ARDEAN LEITH
C              PROMPT                             SEP 2009 ARDEAN LEITH
C              TEXT FILE PROMPT                   NOV 2009 ARDEAN LEITH
C              ! COMMENT DELIMITER                DEC 2011 ARDEAN LEITH
C              'FR N' / INPUT BUG                 APR 2014 ARDEAN LEITH
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
C   SYMPAR(LUNDOC)
C
C   PURPOSE:   CALLS ROUTINES REMOVED FROM DRIVER IN MAR 93
C
C   PARAMETERS: LUNDOC     IO UNIT FOR 'FR F' OPERATION
C
C   CALL TREE:  DRIV1 
C                v          
C               SYMPAR       
C                v 
C       'FR F'   |-> FILESYMPAR -> ----PARSESYMPAR
C                |                     SETSYMPAR
C                |                                     
C       'FR L'   |-> LOCALSYMPAR   
C                |     v            
C          'FR'  |-> RDPRMC ----->  RDPR -> FRSYMPAR  -> PARSESYMPAR                    ^
C                     ^                                  EVALSYMPAR
C                     ^                                  SETSYMPAR
C      ?..? [ID] -----' (FILERD)                              
C                                                   
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE SYMPAR(LUNDOC)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=1)   :: NULL,FCVAL
        CHARACTER (LEN=160) :: PROMPTNID,SYMPARVAL
        LOGICAL             :: CLOSEIT,WANTSUB,BIND

        NULL = CHAR (0)

        SELECT CASE (FCHAR(4:4))

        CASE('F') 
C          GLOBAL VARIABLE & VALUE FROM  TEXT FILE --------------- FR F
           CALL FILESYMPAR(LUNDOC,IRTFLG)

        CASE('N') 
C          NEXT GLOBAL VARIABLE & VALUE FROM  TEXT FILE ---------- FR N
C          READS LINE_BY_LINE
           CLOSEIT = (FCHAR(5:5) .EQ. 'E')
           LUNTEXT = 103
           CALL SYMPARTEXT(CLOSEIT,LUNTEXT,IRTFLG)

        CASE('G')
C          GLOBAL VARIABLE & VALUE FROM  INPUT ------------------- FR G
C          GET GLOBAL  PARAMETER & ASSOCIATED VALUE FROM INPUT 
           BIND = (FCHAR(5:5) .EQ. 'S')
           CALL LOCALSYMPAR(.FALSE.,BIND,SYMPARVAL,IRTFLG)

        CASE('L')
C          LOCAL VARIABLE & VALUE FROM  INPUT -------------------- FR L
C          GET LOCAL  PARAMETER & ASSOCIATED VALUE FROM INPUT 
           BIND = (FCHAR(5:5) .EQ. 'S')
           CALL LOCALSYMPAR(.TRUE.,BIND,SYMPARVAL,IRTFLG)

        CASE DEFAULT
C          FILE READ ----------------------------------------------- FR 
C          GET "?-----? PROMPT,PARAMETER NUMBER, AND ASSOCIATED VALUE
C          FROM CALLER (CALLER CAN BE PROCEDURE OR INTERACTIVE RUN)

C          KEEP LOWERCASE INPUT BY SETTING IRTFLG = -999  
           IRTFLG = -999
           CALL RDPRMC(PROMPTNID,NCHAR,.TRUE.,
     &        'PROMPT (?PROMPT?) & VARIABLE NAME ([NAME])',NULL,IRTFLG)

C          READ AND SET A SYMBOL
           CALL FRSYMPAR(PROMPTNID(1:NCHAR),SYMPARVAL,NCHAR,IRTFLG)

         END SELECT

       RETURN
        END


C++*********************************************************************
C
C  FRSYMPAR.F -- CREATED 6/8/02 ARDEAN LEITH 
C
C **********************************************************************
C
C FRSYMPAR(PROMPTNID,SYMPARVAL,NCHAR,IRTFLG)
C
C PURPOSE: 
C     MEANT TO BE USED INSIDE A PROCEDURE!
C     TAKES IN PROMPT & PARAMETER LABEL, THEN QUERIES CALLING
C     PROCEDURE OR TERMINAL FOR ASSOCIATED VALUE USING THIS PROMPT. 
C     ASSOCIATED VALUE IS USED A LOCAL PARAMETER.
C             
C PARAMETERS:     PROMPTNID    PROMPT AND ID                      SENT
C                 SYMPAROUT    VARIABLE PARAMETER VALUE        RETURNED
C                 NCHARV       LENGTH OF SYMPARVAL             RETURNED
C                 IRTFLG       ERROR FLAG (0 IS NORMAL)        RETURNED
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE FRSYMPAR(PROMPTNID,SYMPAROUT,NCHARV,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 
 
      CHARACTER (LEN=*) ::         SYMPAROUT,PROMPTNID
      CHARACTER (LEN=2*MAXNAM) ::  SYMPARID,PROMPT,SYMPARIN,COMMENTSTR
      CHARACTER (LEN=1) ::         NULL,CDUM
      LOGICAL   ::                 CALLERRT

C     FOR VARIABLE  HANDLING 
      INTEGER, DIMENSION(MAXPRC) :: IPSTACK,IPNUMSTACK,IPARNUM
      COMMON /QSTR_STUFF1/ ISTOP,ITI,ITIN,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM
#ifdef USE_MPI
      include 'mpif.h'
      icomm = MPI_COMM_WORLD
      call MPI_COMM_RANK(icomm, mypid, ierr)
#else
      MYPID = -1
#endif

      NULL = CHAR(0)

C     EXTRACT PROMPT & ID  FROM PROMPTNID INPUT STRING
      CALLERRT = .TRUE.
      CALL PARSESYMPAR(PROMPTNID,NULL,PROMPT,NCHARP,
     &                 SYMPARID,NCHARI,CDUM,NDUM,CALLERRT,IRTFLG)
      IF (PROMPT .EQ. NULL) RETURN

      IF (SYMPARID .EQ. NULL .AND. CALLERRT) THEN
C        MUST CREATE A NUMERICAL LABEL
         INUM           = IPARNUM(ISTOP) + 1
         IPARNUM(ISTOP) = INUM
         SYMPARID(1:1)  = '<'
         CALL INTTOCHAR(INUM,SYMPARID(2:),NCHARI,1)
         SYMPARID(NCHARI+2:NCHARI+2) = '>'
         NCHARI  = NCHARI + 2
c        write(6,*)'symparid(1:',nchari,'): ',symparid
      ENDIF

C     INPUT ASSOCIATED VALUE FOR THIS VARIABLE  


      IF (FROMBATCH) THEN
C        FROM BATCH MODE, NOT FROM INTERACTIVE MODE
C        SO GET SYMPARIN FROM CALLING PROCEDURE FILE

C        INCREMENT BATCH LINE POINTER FOR FURTHER READS
         IPSTACK(ISTOP) = IPSTACK(ISTOP) + 1
         CALL PROC_GETPLINE(IPSTACK(ISTOP),IPNUMSTACK(ISTOP-1),SYMPARIN,
     &                         NCHAR,IRTFLG)

      ELSE
C        '?...?' FROM BATCH TO INTERACTIVE MODE

C        WRITE  ?---? PROMPT TO TERMINAL 
         IF (MYPID .LE. 0) THEN
            WRITE(ITI,991,ADVANCE='NO') PROMPT(1:NCHARP)
         ENDIF
991      FORMAT( ' .',A,': ')

C        GET SYMPARIN FROM CALLING TERMINAL
         READ(ITIN,80) SYMPARIN
 80      FORMAT(A)

      ENDIF

C     STRIP LEADING & TRAILING BLANKS IN SYMPARIN BEFORE COMMENT
      NCHARR = lnblnkn(SYMPARIN)
      CALL PARSE_RESPONSE(SYMPARIN,NCHARR,.TRUE.,.TRUE.,
     &                    SYMPARIN,NCHAR,COMMENTSTR,NCHARC,IRTFLG)

      NLENBRAK = 1
C     LOOP TO CHECK FOR ALL [] PAIRS
      DO WHILE (NLENBRAK .GT. 0)
         CALL CHARINSIDE(SYMPARIN(1:NCHAR),'[',']',.FALSE.,.FALSE.,
     &                   IGOBRAK,IENDBRAK,NLENBRAK)

         IF (NLENBRAK .GT. 0) THEN      
C           CONVERT  [] VARIABLE DELIMITERS TO QSTRQ <> FORMAT
            SYMPARIN(IGOBRAK:IGOBRAK)   = '<'
            SYMPARIN(IENDBRAK:IENDBRAK) = '>'
         ENDIF
      ENDDO

C     SUBSTITUTE FOR VARIABLES & REGISTERS IN HIGHER LEVEL PROCEDURES
      CALL EVALSYMPAR(SYMPARIN(1:NCHAR),SYMPAROUT,NCHARV,IRTFLG)

      IF (CALLERRT) THEN
C        SET VARIABLE AT THIS LEVEL
         CALL SETSYMPAR(SYMPARID(1:NCHARI),SYMPAROUT(1:NCHARV),
     &               .TRUE.,IRTFLG)
      ENDIF

C     WRITE TO  RESULTS FILE
      IF (MYPID .LE. 0) THEN
         WRITE(NOUT,*) ' ',SYMPAROUT(1:NCHARV)
      ENDIF 

#ifdef USE_MPI
      call MPI_BARRIER(icomm,ierr)
#endif
      RETURN
      END



C++*********************************************************************
C
C  FILESYMPAR.F -- CREATED 6/8/02 ARDEAN LEITH 
C
C **********************************************************************
C
C FILESYMPAR(PROMPT,NCHAR,ANS,UPPER,SAYIT,IRTFLG)
C
C PURPOSE: 
C             
C PARAMETERS:     LUNT          UNIT FOR DOC FILE               SENT
C                 IRTFLG        ERROR FLAG (0 IS NORMAL)        RETURNED
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE FILESYMPAR(LUNT,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 
 
      CHARACTER (LEN=MAXNAM) ::  FILNAM,RECLIN,SYMPARID,SYMPARVAL
      CHARACTER (LEN=MAXNAM) ::  COMMENTSTR
      CHARACTER (LEN=1)      ::  CDUM

      
      LENREC = 0
      CALL OPAUXFILE(.TRUE.,FILNAM,DATEXC,LUNT,LENREC,
     &                 'O','TEXT',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      IRTFLG = 1

C     ---------------------- INPUT LOOP 
 
10    READ(LUNT,80,END=998,ERR=999) RECLIN
80    FORMAT(A)

C     IGNORE COMMENT LINES 
      NCHARR = lnblnkn(RECLIN)
      CALL PARSE_RESPONSE(RECLIN,NCHARR,.TRUE.,.TRUE.,
     &                    RECLIN,NCHAROUT,COMMENTSTR,NCHARC,IRTFLG)
      IF (NCHAROUT .LE. 0) GOTO 10

C     CONVERT OLD <> VARIABLE DELIMITER TO NEW: []
      IGOANG = INDEX(RECLIN(1:NCHAROUT),'<')
      IF (IGOANG .GE. 1) THEN
         RECLIN(IGOANG:IGOANG) = '['
         IENDANG = INDEX(RECLIN(1:NCHAROUT),'>')
         IF (IENDANG .GE. 1) THEN
            RECLIN(IENDANG:IENDANG) = ']'
         ENDIF
      ENDIF
 
      IEND = INDEX(RECLIN(1:NCHAROUT),']')
      IF (IEND .LE. 1) THEN
         WRITE(NDAT,*) '*** UNDECIPHERABLE LINE: ',RECLIN(1:NCHAROUT)
         CALL ERRT(101,'FILESYMPAR',NE)
         GOTO 10
      ENDIF

C     EXTRACT VARIABLE ID & VALUES FROM RECLIN
      CALL PARSESYMPAR(CHAR(0),RECLIN(1:NCHAROUT),CDUM,NDUM,
     &                 SYMPARID,NCHARI,
     &                 SYMPARVAL,NCHARV,.TRUE.,IRTFLG)
      IF (SYMPARID .EQ. CHAR(0)  .OR. IRTFLG .NE. 0) GOTO 999 

C     SET GLOBAL VARIABLE ID & VALUE
      CALL SETSYMPAR(SYMPARID(:NCHARI),SYMPARVAL(:NCHARV),
     &                .FALSE.,IRTFLG)

      GOTO 10

C     ------------------ END INPUT LOOP

998   IRTFLG = 0

999   CLOSE (LUNT)
      RETURN
      END


C++*********************************************************************
C
C  LOCALSYMPAR.F -- CREATED 6/8/02 ARDEAN LEITH 
C
C **********************************************************************
C
C LOCALSYMPAR(LOCAL,BIND,SYMPARVAL,IRTFLG)
C
C PURPOSE:  
C             
C PARAMETERS:      
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************


      SUBROUTINE LOCALSYMPAR(LOCAL,BIND,SYMPARVAL,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 
 
      CHARACTER (LEN=*)      :: SYMPARVAL
      CHARACTER (LEN=160)    :: RESPONSE
      CHARACTER (LEN=MAXNAM) :: SYMPARID
      LOGICAL                :: LOCAL,GETANS,STRIP,BIND
      CHARACTER (LEN=1)      :: NULL,CDUM
      LOGICAL                :: UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI

      NULL = CHAR(0)

C     DO NOT UPPERCASE THE INPUT LINE, DO NOT SUBSTITUTE FOR REGS
      GETANS    = .TRUE.
      UPPER     = .FALSE.
      SAYPRMT   = .TRUE.
      SAYANS    = .TRUE.
      ENDATSEMI = .TRUE.
      STRIP     = .TRUE.
      WANTSUB   = .FALSE.

      CALL RDPR('VARIABLE NAME & ASSOCIATED VALUE',NCHAR,RESPONSE,
     &       GETANS,UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
 
C     EXTRACT PROMPT (OLD STYLE) & ID FROM RESPONSE
      CALL PARSESYMPAR(RESPONSE(1:NCHAR),NULL,CDUM,NDUM,SYMPARID,NCHARI,
     &                 CDUM,NDUM,.FALSE.,IRTFLG)
      IF (SYMPARID .EQ. CHAR(0)  .OR. IRTFLG .NE. 0) THEN
         CALL ERRT(101,'SYMPAR',NE)
         RETURN
      ENDIF

C     EXTRACT  SYMBOL VALUE FROM RESPONSE(S)
      CALL PARSESYMPAR(NULL,RESPONSE(1:NCHAR),CDUM,NDUM,CDUM,NDUM,
     &                 SYMPARVAL,NCHARV,.FALSE.,IRTFLG)

      IF (NCHARV .LE. 0) THEN
C        NO SYMBOL VALUE IN RESPONSE,  MUST GET SYMBOL VALUE NOW
         CALL RDPR(SYMPARID(1:NCHARI),NCHARV,SYMPARVAL,
     &       GETANS,UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF

      IF (BIND) THEN
C        SUBSTITUTE FOR STRING VAR. NOW
         CALL SUBSYMPAR(SYMPARVAL(1:NCHARV),SYMPARVAL,NCHARV,0,IRTFLG)
      ENDIF

C     SET LOCAL SYMBOL NAME & VALUE
      CALL SETSYMPAR(SYMPARID(1:NCHARI),SYMPARVAL(1:NCHARV),
     &               LOCAL,IRTFLG)
      
      END


C      *********************** PARSE_RESPONSE ********************************

       SUBROUTINE PARSE_RESPONSE(RESPONSE,NCHARR,ENDATSEMI,STRIP,
     &                           ANSW,NCHAR,COMMENTSTR,NCHARC,IRTFLG)

C      FINDS LOCATION OF COMMENT AND ANY TRAILING BLANKS BEFORE COMMENT

       CHARACTER(LEN=*) :: RESPONSE,ANSW,COMMENTSTR
       LOGICAL          :: ENDATSEMI,STRIP
       LOGICAL          :: KEEPGO
       CHARACTER(LEN=1) :: CTEMP,CTEMPJ

       NCHAR  = 0
       NCHARC = 0

       DO I = 1,NCHARR
          CTEMP = RESPONSE(I:I)

          IF ((CTEMP == ';' .OR. CTEMP == '!') .AND. ENDATSEMI) THEN
             COMMENTSTR = RESPONSE(I:)
             NCHARC     = NCHARR - I + 1
             EXIT

          ELSEIF ((CTEMP .LT. '!' .OR. CTEMP .GT. '~') .AND.
     &             .NOT. STRIP) THEN
C            GOT NON PRINTING CHAR LIKE A BLANK
             NCHAR = NCHAR + 1
             ANSW(NCHAR:NCHAR) = CTEMP   ! DO NOT REPLACE WITH BLANK

          ELSEIF ((CTEMP .GE. '!' .AND. CTEMP .LE. '~')) THEN
C            GOT PRINTING CHAR
             NCHAR             = NCHAR + 1
             ANSW(NCHAR:NCHAR) = CTEMP

          ELSEIF ((CTEMP .LT. '!' .OR. CTEMP .GT. '~') .AND.
     &             NCHAR .GT. 0 .AND. I .LT. NCHARR) THEN
C            GOT NON PRINTING CHAR LIKE A BLANK AFTER A PRINTING CHAR
             KEEPGO = .FALSE.
             DO J = I+1,NCHARR
                CTEMPJ = RESPONSE(J:J)
                IF ((CTEMP .GE. '!' .AND. CTEMP .LE. '~')) THEN
C                  GOT PRINTING CHAR
                   KEEPGO = .TRUE.
                   EXIT
                ENDIF
             ENDDO
             IF (.NOT. KEEPGO) EXIT
             NCHAR             = NCHAR + 1
             ANSW(NCHAR:NCHAR) = CTEMP
           ENDIF
       ENDDO
       IRTFLG = 0
       END
  

C++*********************************************************************
C
C SYMPARTEXT                  NEW                 OCT 2006 ARDEAN LEITH
C
C **********************************************************************
C
C  SYMPARTEXT(CLOSEIT,LUNT,IRTFLG)
C
C  PURPOSE:  SUPPORTS OPERATION TO RETRIEVE A SYMBOLIC VARIABLE FROM
C            A TEXT FILE. ALWAYS SOLICITS FILENAME,  OPENS FILE IF NAME
C            NOT SAME AS PREVIOUS TEXT FILE USED BY THIS OPERATION.
C             
C  TYPICAL USAGE: 
C            FR N
C            Filename            (RETRIEVED)
C            FR NE               (CLOSES FILE OPENED WITH FR N)
C
C  PARAMETERS:    CLOSEIT   CLOSE CURRENT FILE                 (SENT)
C                 LUNT      LUN NUMBER OF FILE                 (SENT)
C                 IRTFLG    ERROR RETURN FLAG                  (RET.)
C
C--*********************************************************************

	SUBROUTINE SYMPARTEXT(CLOSEIT,LUNT,IRTFLG)

        INCLUDE 'CMLIMIT.INC' 
        INCLUDE 'CMBLOCK.INC'

        LOGICAL               :: CLOSEIT
        CHARACTER(LEN=MAXNAM) :: FRNAME,FRNAMET
        CHARACTER(LEN=MAXNAM) :: OLDNAM = '-'
        CHARACTER(LEN=MAXNAM) :: SYMPARID,SYMPARVAL
        INTEGER               :: LUNOLD = 0
        INTEGER               :: NLINE  = 0
        CHARACTER(LEN=1)      :: NULL

        NULL = CHAR(0)

        CALL REG_GET_USED(NSEL_USED)

C       GET VARIABLE LIST NAME
        IRTFLG = -999    ! CONVERT LEGACY REGISTERS x**
	CALL FILERD(FRNAMET,NLET,NULL,'VARIABLE LIST',IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            FIRTFLG = IRTFLG
            CALL REG_SET_NSEL(1,1,FIRTFLG,0.0,0.0,0.0,0.0,IRTFLG)
            RETURN
        ENDIF

        IF (CLOSEIT) THEN
C          END USE OF CURRENT FILE
           CLOSE(LUNOLD)
           LUNOLD = 0
           NLINE  = 0
           OLDNAM = NULL
           IRTFLG = 0
           RETURN
        ENDIF

        IF (LUNT .NE. LUNOLD .OR. FRNAMET .NE. OLDNAM) THEN
C         NEW FILE OR DIFFERENT LUN FOR FILE
          IF (LUNOLD .NE. 0) CLOSE(LUNOLD)
          LUNOLD = 0
          OLDNAM = NULL
          NLINE  = 0

C         OPEN THE DOC FILE
          FRNAME = FRNAMET
          LENREC = 0             !SEQUENTIAL ACCESS, FORMATTED
          CALL OPAUXFILE(.FALSE.,FRNAME,DATEXC,LUNT,LENREC,
     &                       'O', NULL,.TRUE.,IRTFLG)
          IF (IRTFLG .NE. 0) THEN
             FIRTFLG = IRTFLG
             CALL REG_SET_NSEL(1,1,FIRTFLG,0.0,0.0,0.0,0.0,IRTFLG)
             RETURN
          ENDIF

          OLDNAM = FRNAMET
          LUNOLD = LUNT
        ENDIF

C       SPECIAL CODE IN RDPRMC FOR NO SUBSTITUTION OF VARIABLE
        IRTFLG = -999
        CALL RDPRMC(SYMPARID,NCHAR,.TRUE.,
     &         'VARIABLE NAME (ENCLOSED WITH[])',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           FIRTFLG = IRTFLG
           CALL REG_SET_NSEL(1,1,FIRTFLG,0.0,0.0,0.0,0.0,IRTFLG)
           RETURN
        ENDIF

        NLINE = NLINE + 1
        READ(LUNT,'(A)',IOSTAT=IERR) SYMPARVAL

        IF (IERR .NE. 0) THEN
C          ERROR ON READ, PROBABLY END OF FILE
           NLETT   = lnblnkn(OLDNAM)

           IF (NSEL_USED > 0) THEN
              FIRTFLG = IERR
              CALL REG_SET_NSEL(1,1,FIRTFLG,0.0,0.0,0.0,0.0,IRTFLG)
           ELSE
              WRITE(NOUT,90) SYMPARID(1:NCHAR),OLDNAM(1:NLETT),NLINE
90            FORMAT(' *** UNABLE TO RETRIEVE: ',A,
     &               '  FROM: ',A,
     &               '  LINE: ',I6)
              CALLERRT(102,'UNABLE TO READ INPUT LINE',NLINE)
           ENDIF
           IRTFLG  = 1
           RETURN
        ENDIF

C       SET THE VARIABLE
        LENID = lnblnk(SYMPARID)
        SYMPARID(1:1) = '<'
        SYMPARID(LENID:LENID) = '>'
        LENVAR = lnblnk(SYMPARVAL)

c        write(6,*) ' FOR ID: ',SYMPARID(1:LENID),
c     &             '  VALUE: ',SYMPARVAL(:LENVAR)

        CALL SETSYMPAR(SYMPARID(1:LENID),SYMPARVAL(:LENVAR),
     &                 .FALSE.,IRTFLG)

C       DO NOT CLOSE FILE UNTIL 'FR NE' IS GIVEN!

	END




