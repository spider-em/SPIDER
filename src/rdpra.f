
C++************************************************************************ 
C
C   RDPRA.F      REMOVED FROM RDPRI                FEB 99  ArDean Leith
C                LENGTHENED ANSW                 11/17/00  ArDean Leith
C                EXPRESSQ SHOULD BE SUB.           MAY 01  ArDean Leith
C                DID NOT PRINT IF X                JAN 02  ArDean Leith
C                ~PROMPT                           FEB 03  ArDean Leith
C                NLOG                            11/26/03  ArDean Leith
C                RDPR PARAMETERS                 04/14/05  ArDean Leith
C                ?..? LEVELS                     11/28/05  ArDean Leith
C                [] REWRITE                      12/02/05  ArDean Leith
C                LEGACY () INPUT                 06/25/06  ArDean Leith
C                SUBSYMPAR BUG                   07/20/06  ArDean Leith
C                FORMAT(..1PG12.3, MPISET        11/04/10  ArDean Leith
C                NECHO                           09/20/12  ArDean Leith
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
C   RDPRA(PROMPT,NVAL,INTS,VALUES,NGOT,IRTFLG)
C
C   PURPOSE:    EVALUATE INPUT FOR RDPR** SUBROUTINES
C
C   PARAMETERS:
C        PROMPT      PROMPT                                       (SENT)
C        NVAL        NUMBER OF VALUES TO RETURN                   (SENT)
C        ILEVEL      STACK LEVEL FOR REGISTER VARIABLES           (SENT)
C        INTS        LOGICAL FLAG FOR INTEGER RETURN              (SENT)
C        VALUES      VALUES (ALWAYS AS FLOATS)                (RETURNED)
C        NVAL        NUMBER OF VALUES RETURNED                (RETURNED)
C        IRTFLG      SENT TO RDPRA AS FLAG!!                      (SENT)
C                    FLAG = 0 IS NORMAL                       (RETURNED)
C                    FLAG = 1 IS ERROR                        (RETURNED)
C                    FLAG = -1 IS GOTO PREVIOUS QUESTION      (RETURNED)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE RDPRA(PROMPT,NVAL,ILEVEL,INTS,VALUES,NGOT,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 

      CHARACTER(LEN=*)         :: PROMPT
      REAL, DIMENSION(NVAL)    :: VALUES
      LOGICAL                  :: INTS

C     MAXANS IS LENGTH OF ANSW, RESPONSE
      INTEGER, PARAMETER       :: MAXANS = 600
      CHARACTER(LEN=MAXANS)    :: ANSW,RESPONSE

      LOGICAL                  :: GETANS,UPPER,WANTSUB
      LOGICAL                  :: SAYPRMT,SAYANS,ENDATSEMI,STRIP
      LOGICAL                  :: INPARLOOP

      INTEGER, PARAMETER       :: MAXB = 200
      REAL, DIMENSION(MAXB)    :: FBUF

      SAVE            FBUF

      CALL SET_MPI(ICOMM,MYPID,MPIERR)

C     NOFF IS OFFSET FOR NUMBER OF INPUTTED MULTILINE INPUT VALUE
      NOFF = 1

      GETANS    = .TRUE.
      UPPER     = .FALSE.
      WANTSUB   = .FALSE.
      SAYPRMT   = .TRUE.
      SAYANS    = .FALSE.
      ENDATSEMI = .TRUE.
      STRIP     = .TRUE.

10    IF (PROMPT(:1) == '~') THEN
C         USE PROMPT FOR INPUT LINE
          ANSW    = PROMPT(2:)
          NCHAR   = lnblnkn(ANSW)
      ELSE
          CALL RDPR(PROMPT,NCHAR,ANSW,GETANS,
     &              UPPER,WANTSUB,SAYPRMT,SAYANS,
     &              ENDATSEMI,STRIP,IRTFLG)
         IF (IRTFLG == -1) RETURN
      ENDIF
        
      !write(6,*) ' in rdpra, got:',nchar,':',answ(1:16)

      IF (NCHAR <= 0) THEN
         NGOT   = 0
         IRTFLG = 0
         !write(6,*) 'rdpra, returning:',ngot,irtflg
         RETURN
      ENDIF

      DO WHILE (ANSW(NCHAR:NCHAR) == ',') 
C        INPUT CONTINUATION LINE 
         CALL RDPR('NEXT LINE OF INPUT',NCHAR2,
     &             RESPONSE,GETANS,
     &             UPPER,WANTSUB,SAYPRMT,SAYANS,
     &             ENDATSEMI,STRIP,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF ((NCHAR + NCHAR2) > MAXANS) THEN
            CALL ERRT(102,'ANSWER TOO LONG',NCHAR + NCHAR2)
            RETURN
         ENDIF
         ANSW(NCHAR+1:) = RESPONSE(1:NCHAR2)
         NCHAR          = NCHAR + NCHAR2
      ENDDO

      IF (ANSW(:1) == '?') THEN
C         ?PROMPT? [var] LINE RETURNED, FIND NEXT RESPONSE
          CALL STACK_RESPONSE(ANSW(1:NCHAR),RESPONSE,IRTFLG)

C         STRIP SEMICOLON DENOTED COMMENT & TRAILING BLANKS
          CALL DECOMMENT(RESPONSE,NCHAR,LOCSEMI)

C         SEE IF NEED TO CONVERT OLD x11 REGISTER FORMAT 
          IX = SCAN(RESPONSE(1:NCHAR),'xX')
          IF (IX > 0) THEN
             IF (NCHAR .LT. MAXANS) RESPONSE(NCHAR+1:) = ' '
C            CONVERT OLD x11 REGISTER FORMAT TO TO NEW: [name] FORMAT
             CALL DEXREG(RESPONSE,NCHAR)
          ENDIF
          ANSW   = RESPONSE(1:NCHAR)
          ISTACK = -1
      ELSE
          ISTACK = ILEVEL
      ENDIF

C     THIS SHOULD BE DONE IN RDPR BUT IT FOULS UP BY
C     SUBSTITUTING FOR: {***[]}   {---[]}    ***[]   ${ENV}   .1[] 
      IGOBRAK = INDEX(ANSW(1:NCHAR), '[') 
      IF (IGOBRAK > 0) THEN
C        HAS '[...]' NEED SYMBOL SUBSTITUTION E.G. [str]
         CALL SUBSYMPAR(ANSW(1:NCHAR),ANSW,NCHAR,0,IRTFLG)
      ENDIF

      !write(6,*) ' in rdpra, got1:',nchar,':',answ(1:16)

C     EVALUATE EXPRESSION(S) OR LIST OF VALUES FROM ANSW
      CALL EXPRESS3Q(ANSW(1:NCHAR),ISTACK,MAXB,
     &                  FBUF,INUM,INPARLOOP,IRTFLG)

      IF (IRTFLG .NE. 0) THEN
C        GOT BAD INPUT PARAMETERS, RE-ENTER
         CALL ERRT(16,'CAN NOT INTERPRET INPUT',NE)
         GOTO 10
      ENDIF
      !write(6,*) ' in rdpra, got2:',nchar,':',answ(1:16)

C     PREVIOUSLY ANY  EXPRESSION WITH REGISTERS DID NOT NEED ()
      IBRAK = SCAN(ANSW(1:NCHAR),'[')
      IF (IBRAK > 0) INPARLOOP = .TRUE.

      IF (.NOT. LEGACYPAR .OR. INPARLOOP .OR. NLOOP .LE. 1) THEN
C        INPUT HAS () AROUND IT OR IS NOT IN A MULTIPLE VALUE LOOP.
C        IF IN A LOOP, USES SAME INPUT FOR ALL LOOP INDICES
         ILOC = 1
C        NGOT IS NUMBER OF VALUES LEFT IN FBUF
         NGOT = INUM
      ELSE
C        INPUT HAS NO () AROUND IT, USES DIFFERENT SET OF INPUTS
C        FOR EACH INDEX OF THE CURRENT LOOP
 
         NTOT   = NOFF + INUM - 1 
         NEEDED = NVAL * NLOOP 
         IF (NTOT < NEEDED)  THEN
C           NEEDS MORE INPUT TO GET "NVAL" INPUTS, READ ANOTHER LINE.
C           INCREMENT CURRENT INPUT VALUE INDEX FIRST
            NOFF = NOFF + INUM
            GOTO 10
         ENDIF

C        ILOC IS POINTER TO CURRENT LOCATION IN FBUF
         ILOC = (ILOOP - 1) * NVAL + 1
C        NGOT IS NUMBER OF VALUES LEFT IN FBUF
         NGOT = NTOT - ILOC + 1
      ENDIF
 
      NGOT = MIN(NGOT,NVAL)
      IF (INTS) THEN
C        CONVERT VALUES TO INTEGERS
         DO I = 1,NGOT
            VALUES(I) = INT(FBUF(ILOC + I -1))
         ENDDO
      ELSE
C        FLOATING POINT VALUES WANTED
         DO I = 1,NGOT
            VALUES(I) = FBUF(ILOC + I -1)
         ENDDO
      ENDIF
  
      IF (INTS .AND. MYPID <= 0) THEN
C        INTEGER VALUES WANTED
         IF (NOUT .NE .0) WRITE(NOUT,90) (INT(VALUES(I)),I=1,NGOT)
         IF (NLOG .NE .0) THEN
            WRITE(NLOG,90) (INT(VALUES(I)),I=1,NGOT)
            NECHO = NECHO + 1
         ENDIF

 90      FORMAT(2X,10(1X,I7))
         
      ELSEIF (MYPID <= 0) THEN
C        FLOATING POINT VALUES WANTED
         IF (NOUT .NE .0) WRITE(NOUT,91) (VALUES(I),I=1,NGOT)
         IF (NLOG .NE .0) THEN
            WRITE(NLOG,91) (VALUES(I),I=1,NGOT)
            NECHO = NECHO + 1
         ENDIF
  91     FORMAT(2X,10(1PG12.3,1X))
      ENDIF
 
      IRTFLG = 0

      RETURN
      END


C      *********************** STACK_RESPONSE *************************

      SUBROUTINE STACK_RESPONSE(PROMPTNID,RESPONSE,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 
 
      CHARACTER (LEN=*)        :: PROMPTNID,RESPONSE

      CHARACTER (LEN=2*MAXNAM) :: PROMPT,SYMPARID
      CHARACTER (LEN=1)        :: NULL,CDUM

C     FOR LOCAL VARIABLE HANDLING 
      INTEGER, DIMENSION(MAXPRC) :: IPSTACK,IPNUMSTACK,IPARNUM
      COMMON /QSTR_STUFF1/ ISTOP,ITI,ITIN,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      CALL SET_MPI(ICOMM,MYPID,MPIERR)

      NULL = CHAR(0)

C     EXTRACT PROMPT FROM PROMPTNID INPUT STRING
      CALL PARSESYMPAR(PROMPTNID,NULL,PROMPT,NCHARP,
     &                 SYMPARID,NCHARI,CDUM,NDUM,CALLERRT,IRTFLG)
      IF (PROMPT == NULL) RETURN

C     WRITE PROMPT TO  RESULTS FILE
      IF (MYPID .LE. 0) THEN
         WRITE(NOUT,*) ' ',PROMPT(1:NCHARP)
      ENDIF 

      IF (FROMBATCH) THEN
C        FROM BATCH MODE, NOT FROM INTERACTIVE MODE
C        SO GET RESPONSE FROM CALLING PROCEDURE FILE

C        INCREMENT BATCH LINE POINTER FOR FURTHER READS
         IPSTACK(ISTOP) = IPSTACK(ISTOP) + 1
         CALL PROC_GETPLINE(IPSTACK(ISTOP),IPNUMSTACK(ISTOP-1),
     &                      RESPONSE, NCHAR,IRTFLG)

      ELSE
C        '?...?' FROM BATCH TO INTERACTIVE MODE

C        WRITE  ?---? PROMPT TO TERMINAL 
         IF (MYPID .LE. 0) THEN
            WRITE(ITI,991,ADVANCE='NO') PROMPT(1:NCHARP)
         ENDIF
991      FORMAT( ' .',A,': ')

C        GET RESPONSE FROM CALLING TERMINAL
         READ(ITIN,80) RESPONSE
80       FORMAT(A)
      ENDIF

C     STRIP SEMICOLON DENOTED COMMENT & TRAILING BLANKS
      CALL DECOMMENT(RESPONSE,NCHAR,LOCSEMI)

#ifdef USE_MPI
      call MPI_BARRIER(icomm,ierr)
#endif
      RETURN
      END

