
C++*********************************************************************
C
C  RDPR.F -- CREATED 2/8/90 ARDEAN LEITH 
C	  -- ADD ON-LINE HELP                3/29/93 MAHIEDDINE LADJADJ  
C         -- CONVERTED FROM READCH                  DEC 96 ARDEAN LEITH
C         -- F90 CHANGES                            OCT 97 ARDEAN LEITH
C         -- STRIPS COMMENT                         AUG 99 ARDEAN LEITH
C         -- LUNDONOW ADDED                         OCT 99 ARDEAN LEITH
C         -- TRAILING BLANKS IN COMMENT REMOVED     NOV 99 ARDEAN LEITH
C         -- PUT IN <1> VARIABLE HANDLING           SEP 00 ARDEAN LEITH
C         -- MULTIPLE VARIABLE SUBSTITUTION         JAN 01 ARDEAN LEITH
C         -- USED PROC_GETLINE                      JAN 01 ARDEAN LEITH
C         -- FLAG FOR ; OK                          MAR 01 ARDEAN LEITH
C         -- ADDED FILNAMSUB                        APR 01 ARDEAN LEITH
C         -- ADDED VERBOSE FOR ;                    APR 01 ARDEAN LEITH
C         -- DELAYED PROMPT FOR .NOT. VERBOSE       JUN 01 ARDEAN LEITH
C         -- MOVED SSUPCASE LATER                   SEP 01 ARDEAN LEITH
C         -- NO PROMPT FOR .OP COMMENT LINES        MAR 02 ARDEAN LEITH
C         -- SYMPAR REWRITTEN                       JUN 02 ARDEAN LEITH
C         -- NO SYMPAR FOR 'RR'                     AUG 02 ARDEAN LEITH
C         -- '[]' --> '<>'                          SEP 02 ARDEAN LEITH
C         -- PARAMETERS CHANGED                     APR 05 ARDEAN LEITH
C         -- [] DEFAULT FOR VARIABLES               OCT 05 ARDEAN LEITH
C         -- NDOLINE                                MAY 07 ARDEAN LEITH
C         -- ?prompt?[  FR BUG                      JUN 07 ARDEAN LEITH
C         -- REMOVED IMCx33 OBSOLETE SYNTAX         JUN 09 ARDEAN LEITH
C	  -- REMOVED ON-LINE HELP                   AUG 09 ARDEAN LEITH  
C         -- $DATEXT x11 BUG                        AUG 09 ARDEAN LEITH
C         -- '@@' SUPPORT                           NOV 09 ARDEAN LEITH
C         -- NDOLINE                                NOV 09 ARDEAN LEITH
C         -- VMS COMMAND DOES NOT <> --> []         SEP 10 ARDEAN LEITH
C         -- ! COMMENT DELIMITER                    DEC 11 ARDEAN LEITH
C         -- RECURSIVE FILNAMSUB                    MAR 12 ARDEAN LEITH
C         -- IRTFLG 654321                          MAR 12 ARDEAN LEITH
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
C  RDPR(PROMPT,NCHAR,ANS,UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,IRTFLG)
C
C  PURPOSE: OUTPUTS PROMPT
C           READS AN ALPHANUMERIC STRING FROM STORED PROC. LINE, TERMINAL,
C              OR PROMPT.
C           ECHO & SKIP LINES WHICH ONLY CONTAIN A COMMENT AT START OF LINE
C           HANDLES INTERACTIVE HELP
C           CAN ECHO LINE TO CURRENT INTERACTIVE DO-LOOP IFLE.
C           CONVERTS OLD @B01[X11] PROC. ARG. FORMAT TO TO NEW: () ARG.  
C           CONVERTS OLD <> VARIABLE FORMAT TO NEW [] VARIABLE FORMAT
C           CAN INVOKE VARIABLE SUBSTITUTION FOR [string]. 
C           CAN CONVERT TO UPPERCASE
C           CONVERTS OLD X REGISTER TO  [] VARIABLE FORMAT
C           SUBSTITUTES FOR {***[]} AND ${ENV} STRINGS
C               
C           RETURNS NCHAR=LENGTH OF STRING WITHOUT TRAILING BLANKS OR COMMENT. 
C           COMMENT IS LIMITED TO 80 CHAR.
C           VARIABLE VALUE RESPONSE IS LIMITED TO 160 CHAR.
C             
C           REGISTER SUBSTITUTION OCCURS IN RDPRINC
C            
C  PARAMETERS:  PROMPT    INPUT PROMPT                     (SENT)
C               NCHAR     LAST NON_BLANK CHAR IN           (RETURNED)
C                            ANS RESPONSE BEFORE COMMENT
C               ANS       USER RESPONSE                    (RETURNED)
C               GETANS    READ ANSWER (NOT PROMPT)         (SENT)
C               UPPER     CONVERT TO UPPERCASE             (SENT)
C               WANTSUB   WANT SYM. PARAMETER SUBSTITUTION (SENT)
C                             HERE NOW (USUAL)
C               SAYPRMT   ECHO PROMPT TO OUTPUT            (SENT)
C               SAYANS    ECHO RAW ANSWER TO OUTPUT        (SENT)
C               ENDATSEMI IGNORE SEMICOLON COMMENT         (SENT) 
C                            (FOR vms.f)
C               IRTFLG    RETURN FLAG (0 IS NORMAL)        (RETURNED)
C
C  CALLED BY:   RDPRMC -> RDPR -> SUBSYMPAR &  SSUPCAS & FILNAMSUB
C
C               RDPRM2 -> RDPRINC -> RDPRA -> RDPR -> SUBSYMPAR & 
C                                                     SSUPCAS & 
C                                                     FILNAMSUB
C                               ---> EXPRESS3Q
C                               ---> CHKSTR
C
C               SPIDER -> RDPRMC
C
C               FILERD --> RDPR & FILNAMSUB
C               INQUIREREG
C               VMS,VMS_CD
C               SYMPAR, UTIL4
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE RDPR(PROMPT,NCHAR,ANS,
     &       GETANS,UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      COMMON /LUNDOECHO/ LUNDONOW,NDOLINE

      CHARACTER(LEN=*)   :: PROMPT, ANS
      CHARACTER(LEN=1600):: ANS_CYCLE
      CHARACTER(LEN=80)  :: COMMENTSTR
      CHARACTER(LEN=1 )  :: CTEMP
      LOGICAL            :: LDUM
      LOGICAL            :: GETANS,UPPER,WANTSUB,SAYPRMT,SAYANS
      LOGICAL            :: ENDATSEMI,STRIP,LEGACYREGS,ACCEPTCR 

      CALL SET_MPI(ICOMM,MYPID,MPIERR)

      LEGACYREGS = (IRTFLG .NE. -999)  ! DO NOT CONVERT x**
      ACCEPTCR   = (IRTFLG == 654321)  ! ACCEPT <CR> FOR RDPRI & M
  
      IDOL = INDEX(PROMPT,'$') - 1
      IF (IDOL .LE. 0) IDOL = LEN(PROMPT)

      IRTFLG = 0

10    CONTINUE
C     PROMPT OUTPUT IS DELAYED IN BATCH TO IGNORE COMMENT / BLANK LINES
      IF (SAYPRMT .AND. COPT == 'I' .AND. MYPID <= 0) THEN
         WRITE(NOUT,90,ADVANCE='NO') PROMPT(1:IDOL)
 90      FORMAT(' .',A,': ')
      ENDIF

      IF (GETANS) THEN
C        READ ANSWER STRING

         IF (NIN == 1) THEN
C           UPDATE THE BATCH COUNTER FOR CURRENT PROCEDURE LINE
            IBCNT = IBCNT + 1

C           READ FROM CURRENT STORED PROCEDURE LINE IBCNT
            CALL PROC_GETPLINE(IBCNT,0,ANS,NCHAR,IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               CALL ERRT(101,'PROCEDURE LACKS: RE',NE)
               ANS    = 'RE'
               NCHAR  = 2
               RETURN
            ENDIF

         ELSE
C           UPDATE THE BATCH COUNTER FOR CURRENT PROCEDURE LINE
            IBCNT = IBCNT + 1

C           READ FROM TERMINAL
            READ(NIN,FMT='(A)',IOSTAT=IERR) ANS
            IF (IERR .NE. 0) THEN
               CALL ERRT(101,'RAN OUT OF INPUT',NE)
               ANS    = 'EN'
               NCHAR  = 2
               RETURN
            ENDIF
            NCHAR = lnblnk(ANS)
 
         ENDIF
      ELSE
C        READ FROM PROMPT INSTEAD OF FROM INPUT
         ANS   = PROMPT
         NCHAR = IDOL
      ENDIF
      IF (NCHAR <= 0) RETURN

C     SEE IF THIS IS A COMMENT ONLY LINE WITH ; IN FIRST POSITION
C     (IF ; IS PROCEEDED BY SPACE MAYBE THE USER INPUT A BLANK??)
      LOCSEMI = SCAN(ANS(1:NCHAR),';!')
      NCHARA  = NCHAR           ! NOTHING BEFORE ; FLAG
      IF (LOCSEMI > 0) NCHARA = lnblnk(ANS(1:LOCSEMI-1))

      !write(6,*) 'rdpr, nchara,locsemi:',nchara,locsemi,acceptcr
      !write(6,*) 'rdpr, ANS(NCHARA:NCHARA):',ANS(NCHARA:NCHARA),':'

      IF (ACCEPTCR .AND. 
     &   (NCHARA == 0 .OR. ANS(NCHARA:NCHARA) == '*') ) THEN

C        ACCEPT <CR> or *  RESPONSE BEFORE COMMENT
         IF (MYPID <= 0 .AND. SAYPRMT .AND. COPT == 'B') THEN
            WRITE(NOUT,95)  PROMPT(1:IDOL),': ',ANS(1:NCHAR)
95          FORMAT('  ',A,A,A)
         ENDIF
         NCHAR  = 0
         IRTFLG = 0
         !write(6,*) 'rdpr, nchar:',nchar,irtflg
         RETURN

      ELSEIF (LOCSEMI == 1 .AND. ENDATSEMI) THEN
C        NOTHING BEFORE COMMENT
         IF (VERBOSE .AND. MYPID <= 0) THEN
C           ECHO COMMENT
            IF (SAYPRMT .AND. COPT == 'B') 
     &          WRITE(NOUT,90,ADVANCE='NO')  PROMPT(1:IDOL)
	    IF (NOUT .NE. 0)  WRITE(NOUT,91) ANS(1:NCHAR)
91          FORMAT(' ',A)
         ENDIF
C        READ ANOTHER INPUT LINE
         GOTO 10 
      ENDIF

      IF (STRIP) THEN
C        REMOVE LEADING AND TRAILING NON-PRINTING CHAR. FROM ANSWER
         I      = 1
         J      = 0
         DO WHILE (I <= NCHAR)
            CTEMP = ANS(I:I)
            IF (CTEMP == ';' .OR. CTEMP == '!') THEN
               COMMENTSTR = ANS(I:)
               EXIT
            ELSEIF (J .GT. 0 .OR. 
     &              (CTEMP >= '!' .AND. CTEMP <= '~')) THEN
               J        = J + 1
               ANS(J:J) = ANS(I:I)
            ENDIF
            I = I + 1
         ENDDO
         NCHAR         = lnblnkn(ANS(1:J))
         ANS(NCHAR+1:) = ' '
         NCHARCOM      = lnblnkn(COMMENTSTR)
      ENDIF
                      
      IQUES = INDEX(ANS(:NCHAR),'?')
C     PROMPT OUTPUT IS DELAYED UNTIL HERE IN BATCH TO IGNORE BLANK LINES 
      IF (SAYPRMT .AND. COPT == 'B' .AND. IQUES <= 0) THEN
         IF (MYPID <= 0) THEN
            WRITE(NOUT,94,ADVANCE='NO') PROMPT(1:IDOL)
         ENDIF
94       FORMAT(' .',A,': ')
      ENDIF

C     HANDLE OBSOLETE INTERACTIVE HELP
      IF ( COPT .EQ. 'I'  .AND. 
     &     (IQUES > 0    .OR. 
     &     (INDEX(ANS(:NCHAR),'HELP')     > 0)  .OR.
     &     (INDEX(ANS(:NCHAR),'help')     > 0)) .AND.
     &     (INDEX(ANS(IQUES+1:NCHAR),'>') == 0) .AND.
     &     (INDEX(ANS(IQUES+1:NCHAR),'[') == 0)) THEN
         IF (LOCSEMI <= 0 .OR. LOCSEMI > IQUES) THEN
            WRITE(NOUT,*)' USE YOUR WEB BROWSER FOR SPIDER MANUAL PAGES'
C           READ ANOTHER INPUT LINE
            GOTO 10
         ENDIF            
      ENDIF

      IF (LUNDONOW > 0 .AND. MYPID <= 0) THEN
C        MUST COPY INPUT LINE TO CURRENT INTERACTIVE DO-LOOP FILE
         WRITE(LUNDONOW,*) ANS(1:NCHAR)
         NDOLINE = NDOLINE + 1
         !write(6,*) ' rdpr lundo: ',ANS(1:NCHAR),':',NDOLINE !!!!
      ELSEIF ((COPT == 'I' .AND. NDOLINE > 0)) THEN
         NDOLINE = NDOLINE + 1
         !write(6,*) ' rdpr ndoline: ',ANS(:NCHAR),':',NDOLINE !!!!
      ENDIF

C     CHECK IF JUST BLANKS BEFOR ; & STRIP OFF ANY TRAILING BLANKS
      IF (LOCSEMI > 0 .AND. ENDATSEMI) THEN 
C        PRESERVE COMMENT FOR LATER USE
         COMMENTSTR = ANS(LOCSEMI:)
         NCHAR      = LNBLNKN(ANS(1:LOCSEMI-1))
      ELSEIF(.NOT. ENDATSEMI) THEN 
         LOCSEMI = 0
      ENDIF

      IF (SAYANS) THEN
C        ECHO ANSWER IN RAW FORMAT
         WRITE(NOUT,*) ' ',ANS(1:NCHAR)
      ENDIF
       
C     CONVERT OLD <> VARIABLE FORMAT TO NEW [] VARIABLE FORMAT
      NLENANG = 1
      DO WHILE (NLENANG > 0 .AND. ENDATSEMI)
         CALL CHARINSIDE(ANS(1:NCHAR),'<','>',.FALSE.,.FALSE.,
     &                   IGOANG,IENDANG,NLENANG)

         IF (NLENANG > 0) THEN      
C           CONVERT OLD <> VARIABLE FORMAT TO NEW [] VARIABLE FORMAT
C           write(6,*) 'CONVERT OLD <> VAR. DELIMIT. TO NEW:',ans 
            ANS(IGOANG:IGOANG)   = '['
            ANS(IENDANG:IENDANG) = ']'      !  MAY BE MORE VARIABLES
         ENDIF
      ENDDO

C     SEE IF '[' AND ']' NEED SYMBOL SUBSTITUTION
      IGOBRAK = INDEX(ANS(1:NCHAR), '[') 
      IF (IGOBRAK > 0 .AND. WANTSUB) THEN
C         '[' AND ']' NEED SYMBOL SUBSTITUTION E.G. [str]
          CALL SUBSYMPAR(ANS(1:NCHAR),ANS,NCHAR,0,IRTFLG)
      ENDIF
         
C     SEE IF NEED TO CONVERT OLD x11 REGISTER FORMAT 
      IX = SCAN(ANS(1:NCHAR),'xX')
c!!!  IF (LEGACYREGS .AND. IX > 0) THEN

      IF (IX .GT. 0 .AND. NCHAR > IX) THEN
C        POSSIBLE Xdd OR Xd IN INPUT
C        CONVERT OLD x11 REGISTER FORMAT TO TO NEW: [name] FORMAT
         !write(6,*) ' calling dexreg:',ans(1:nchar),':'
         CALL DEXREG(ANS,NCHAR)
      ENDIF
 
      IF (WANTSUB) THEN
         IF (NCHAR > 1600) THEN
            CALL ERRT(102,'RESPONSE OVERFLOW',NCHAR)
            RETURN
         ENDIF
         
         DO 
            ISUB   = SCAN(ANS(:NCHAR), '{[*$')
            IF (ISUB <= 0) EXIT             ! NO SUB NEEDED

            ANS_CYCLE   = ANS(1:NCHAR)
            NCHAR_CYCLE = NCHAR

C           SUBSTITUTE FOR: {***[]}   {---[]}    ***[]   ${ENV}  .1[] 
            CALL FILNAMSUB(ANS,NCHAR,0,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

C           EXIT IF THERE WERE NO SUBSTITUTIONS
            IF (NCHAR       == NCHAR_CYCLE .AND.
     &          ANS(:NCHAR) == ANS_CYCLE(:NCHAR)) EXIT

         ENDDO
      ENDIF
 
      IF (UPPER) THEN
C        CONVERT INPUT STRING TO ALL UPPER CASE 
         CALL SSUPCAS(ANS(1:NCHAR))
      ENDIF

      IF (LOCSEMI .GT. 0) THEN
C        PUT COMMENT STRING BACK AT END OF INPUT STRING
         ANS = ANS(1:NCHAR) // COMMENTSTR  
      ENDIF
      IRTFLG = 0

      END


C      *********************** DEXREG ********************************

       SUBROUTINE DEXREG(CINPUT,NCHAR)

      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*)   :: CINPUT
      CHARACTER(LEN=161) :: CSUB
      CHARACTER(LEN=1 )  :: CTEMP
      LOGICAL            :: INSUB

C     CONVERT OLD x11 REGISTER FORMAT TO TO NEW: [name] FORMAT
      I     = 1
      J     = 0
      INSUB = .TRUE.
      DO WHILE (I < NCHAR)
         CTEMP = CINPUT(I:I)
         IF (INSUB .AND. (CTEMP == 'X' .OR. CTEMP == 'x')) THEN
C           PROBABLE REGISTER START x or X
            NDIG = VERIFY(CINPUT(I+1:NCHAR),'0123456789')
            IF (NDIG > 0) THEN
               NDIG = NDIG - 1
            ELSE
               NDIG = NCHAR - I
            ENDIF

            IF (NDIG > 0) THEN
               CSUB = '[_' // CINPUT(I+1:I+NDIG) // ']' // CHAR(0)
               CALL SUBCHAR(CSUB(1:NDIG+3),CINPUT,I,I+NDIG,
     &                          NCHAR,IRTFLG)
               I = I + NDIG
            ENDIF 

         ELSEIF (INSUB .AND. CTEMP == '[') THEN
            INSUB = .FALSE.

         ELSEIF (.NOT. INSUB .AND. CTEMP == ']') THEN
            INSUB = .TRUE.

         ENDIF
         I = I + 1
      ENDDO

      END





C      *********************** DECOMMENT ********************************

       SUBROUTINE DECOMMENT(CINPUT,NCHAROUT,LOCSEMI)

C      FINDS LOCATION OF COMMENT AND ANY TRAILING BLANKS BEFORE COMMENT

       CHARACTER *(*) CINPUT

C      IGNORE SEMICOLON DENOTED COMMENT AT END OF CINPUT STRING
       LOCSEMI = SCAN(CINPUT,';!')

       IF (LOCSEMI <= 0) THEN
          NCHAROUT = LNBLNKN(CINPUT)

       ELSEIF (LOCSEMI == 1) THEN
          NCHAROUT = 0

       ELSEIF (LOCSEMI > 1) THEN
C         STRIP COMMENT & TRAILING BLANKS
          NCHAROUT = LNBLNKN(CINPUT(1:LOCSEMI-1))
       ENDIF
       RETURN
       END

