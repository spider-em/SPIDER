C++*********************************************************************
C
C ERRT.F               ALTERED OUTPUT UNITS           SEPT  97 A. Leith
C                      ADDED ERROR 44                 JULY  98 A. Leith
C                      ADDED ERROR 102                JULY  00 A. Leith
C                      ADDED ERROR CLEANUP EXECUTION  APRIL 02 A. Leith
C                      ADDED SPIREOUT                 JUNE  05 A. Leith
C                      ADDED ERR 8                    FEB   06 A. Leith
C                      MPI                            OCT   07 A. Leith
C                      MPI NOW PRINTS FROM ALL PID    OCT   08 A. Leith
C                      DOUBLE OUTPUT BUG FIXED        JUN   09 A. Leith
C                      IEEE                           MAR   11 A. Leith
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
C  ERRT(IERRT,PROG,NE)
C
C  PURPOSE:      WRITE OUT AN ERROR MESSAGE, HALTS IF IN BATCH MODE
C 
C  PARAMETERS:   IERRT   ERROR CODE NUMBER                     (SENT)
C                PROG    STRING ALPHANUMERIC CHARACTERS CONTAINING
C                        PROGRAM NAME OR ERROR MESSAGE (IF IERRT IS
C                        EQUAL TO 101 or 102)                  (SENT)
C                NE      ONLY USED FOR ERROR 46 & 102          (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE ERRT(IERRT,PROG,NE)

#ifndef SP_GFORTRAN
#ifndef USE_MPI
#ifndef __APPLE__
        USE, intrinsic :: ieee_exceptions
#endif
#endif
#endif

      INCLUDE 'CMBLOCK.INC' 

      CHARACTER (LEN=*)   :: PROG
      CHARACTER (LEN=78)  :: MESG
      CHARACTER (LEN=180) :: CSTRING
      CHARACTER (LEN=1)   :: NUL
C     INTEGER*4           :: NE

C     FOLLOWING COMMON BLOCKS ARE USED BY THIS MODULE
C     COMMON /UNITS/, /FUNCTION/, /PROC/

      CALL SET_MPI(ICOMM,MYPID,MPIERR)

      IF (IERRT .EQ. 16) THEN

C         WRITES:  '.BAD INPUT PARAMETER(S). RE-ENTER' TO NOUT, THEN
C                 WILL ALLOW RE-INPUT

          WRITE(NOUT,116)
116       FORMAT(' .BAD INPUT PARAMETER(S). RE-ENTER: ')
          GOTO 2201

      ELSEIF (IERRT .EQ. 8) THEN

C         WRITES:  'ERROR WITH REGISTER VARIABLE --' THEN WRITES
C                 PROG (ERROR MESSAGE)  
C         TO TERMINAL AND RESULTS FILE

          WRITE(NDAT,904) PROG
          IF (NDAT .NE. 6) WRITE(6,904) PROG
904       FORMAT(' *** ERROR UNDEFINED REGISTER VARIABLE: ',A,/)

C         SKIP USUAL MESSAGE OUTPUT
          GOTO 2201

      ELSEIF (IERRT .EQ. 46 .AND. NE .GT. 0) THEN

C         WRITES:  'ERROR MEMORY ALLOCATION FAILED' THEN WRITES
C                 PROG (ERROR MESSAGE) AND NE (ALLOCATION VALUE) 
C         TO TERMINAL AND RESULTS FILE

          WRITE(NDAT,903) PROG,NE
          IF (NDAT .NE. 6) WRITE(6,903) PROG,NE
903       FORMAT(' *** ERROR MEMORY ALLOCATION FAILED -- ',A,': ',I12/)

C         SKIP USUAL MESSAGE OUTPUT
          GOTO 2201


      ELSEIF (IERRT .EQ. 101) THEN
C         WRITES  PROG (ERROR MESSAGE) TO TERMINAL AND RESULTS FILE

          WRITE(NDAT,901) PROG                       ! 3 = NDAT
          IF (NOUT .NE. NDAT) WRITE(NOUT,901) PROG   ! 6 = NOUT
          IF (NOUT .EQ. NDAT) WRITE(0,901) PROG      ! 0 = TERMINAL
901       FORMAT(' *** ERROR: ',A/)

C         SKIP USUAL MESSAGE OUTPUT
          GOTO 2201


      ELSEIF (IERRT .EQ. 102) THEN
C         WRITES  PROG (ERROR MESSAGE) AND NE (INTEGER VALUE) 
C         TO TERMINAL AND RESULTS FILE

          WRITE(NDAT,902) PROG,NE                       ! 3 = NDAT
          IF (NOUT .NE. NDAT) WRITE(NOUT,902) PROG,NE   ! 6 = NOUT
          IF (NOUT .EQ. NDAT) WRITE(0,902) PROG,NE      ! 0 = TERMINAL
902       FORMAT(' *** ERROR -- ',A,': ',I8/)

C         SKIP USUAL MESSAGE OUTPUT
          GOTO 2201

      ELSEIF (IERRT .EQ. 100) THEN
C         SKIP ANY MESSAGE OUTPUT
          GOTO 2201

       ENDIF

C      PROCESS ALL OTHER ERROR CODES (1 - 47)

       SELECT CASE(IERRT)
          CASE (1)
          MESG = 'INCONSISTENT PICTURE DIMENSIONS'

          CASE (2)
          MESG = 'OPERATION NOT CONSISTENT WITH DATA FORMAT'

          CASE (3)
          MESG = 'ROWLENGTH TOO LARGE FOR BUFFER'

          CASE (4)
          MESG = 'OPENING FILE'

          CASE (5)
          MESG = 'NORMALIZING DATA'

          CASE (6)
          MESG = 'INSUFFICIENT BUFFER SPACE'

          CASE (7)
          MESG = 'I/O RECORD NUMBER OUT OF LIMITS'

          CASE (8)
C         THIS ERROR ALREADY PROCESSED BEFORE
          CONTINUE

          CASE (9)
          MESG = 'LABEL SPACE INSUFFICIENT'

          CASE (10)
          MESG = 'DIMENSIONS NOT POWER OF TWO'

          CASE (11)
          MESG = 'PICTURE DIMENSIONS EXCEED FRAME'

          CASE (12)
          MESG = 'READING FILE'

          CASE (13)
          MESG = 'OPERATION NOT ALLOWED IN INTERACTIVE MODE'

          CASE (14)
          MESG = 'INCONSISTENT INPUT PARAMETERS'

          CASE (15)
          MESG = 'FILE IS WRITE PROTECTED'

          CASE (16)
C         THIS ERROR ALREADY PROCESSED BEFORE
          CONTINUE

          CASE (17)
          MESG = 'END-OF-FILE ON INPUT'

          CASE (18)
          MESG = 'FILE DOES NOT EXIST'

          CASE (19)
          MESG = 'TOO MANY FILE NUMBERS ENTERED'

          CASE (20)
          MESG = 'NOT DEFINED...'

          CASE (21)
          MESG = 'NUMBER OF PROCEEDURES EXCEEDS MAXPRC'

          CASE (22)
          MESG = 'TOO MANY PROCEDURES'

          CASE (23)
          MESG = 'UNKNOWN OPTION'

          CASE (24)
          MESG = 'SLICE NUMBER OUT OF BOUNDS'

          CASE (25)
          MESG = 'NOT CONTAINED IN TABLE'

          CASE (26)
          MESG = 'UNDEFINED ERROR'

          CASE (27)
          MESG = 'UNABLE TO CONSTRUCT FILE NAME'

          CASE (28)
          MESG = 'NUMBER OF PROJECTIONS TOO LARGE'

          CASE (29)
          MESG = 'POOR PHASES'

          CASE (30)
           MESG = 'TMPARY FULL, SHORTEN YOUR FILENAMES?'

          CASE (31)
           MESG = 'PARAMETER VALUE OUT OF LEGAL RANGE'

          CASE (32)
           MESG = 'PROGRAM TRAP'

          CASE (33)
           MESG = 'ACCURACY NOT ACHIEVED'

          CASE (34)
           MESG = 'OVERFLOW PROTECTION'

          CASE (35)
           MESG = 'EXIT ERROR'

          CASE (36)
           MESG = 'ILLEGAL ARGUMENT RANGE'

          CASE (37)
           MESG = 'ONLY FIXED IMAGE SIZE ALLOWED'

          CASE (38)
          MESG = 'DIMENSION IS A PRODUCT OF A PRIME NUMBER > 23'

          CASE (39)
          MESG = 'COMMAND NOT SUPPORTED FOR THIS FILE FORMAT'

          CASE (40)
          MESG = 'CONFLICT OF DATA FORMATS'

          CASE (41)
          MESG = 'COMMAND TEMPORARILY DISABLED'

          CASE (42)
          MESG = 'CAN NOT USE THIS OPERATION ON VOLUMES'

          CASE (43)
          MESG = 'INVALID ARITHMETIC EXPRESSION'

          CASE (44)
          MESG = 'INVALID COLUMN, ROW, OR SLICE NUMBER'

          CASE (45)
          MESG = 'INSUFFICIENT DYNAMIC MEMORY'

          CASE (46)
          MESG = 'MEMORY ALLOCATION FAILED'

          CASE (47)
          MESG = 'CANNOT PERFORM FFT'

          CASE DEFAULT
C         PROCESS ERROR CODE > 47
          MESG = 'UNIDENTIFIED,'

        END SELECT


C         WRITE TO TERMINAL AND RESULTS FILE ------------------------

          WRITE(NDAT,900) MESG,IERRT,PROG                      ! 3 = NDAT
          IF (NOUT .NE. NDAT) WRITE(NOUT,900) MESG,IERRT,PROG  ! 6 = NOUT
          IF (NOUT .EQ. NDAT) WRITE(0,900) MESG,IERRT,PROG     ! 0 = TERMINAL
900       FORMAT(' *** ERROR: ',A,'#',I3,' SUBROUTINE: ',A,' ***'/)

          CSTRING = ' *** ERROR: ' // MESG // ' SUBROUTINE: ' // 
     &              PROG // ' ***' 
          CALL SPIREOUT(CSTRING,IRTFLG)


2201      IF (FROMBATCH .OR. COPT .EQ. 'B') THEN
C            IN BATCH MODE, MUST STOP IMMEDIATELY
             NUL = CHAR(0)
             CALL ENDIT('TERMINATED ON ERROR IN BATCH MODE',.FALSE.,NUL)

#ifdef  HAS_IEEE
#ifndef SP_GFORTRAN
#ifndef USE_MPI
#ifndef __APPLE__
C            DO NOT REPORT IEEE INEXACT ....
             call ieee_set_flag( ieee_inexact,.FALSE. )
             call ieee_set_flag( ieee_invalid,.FALSE. )
             call ieee_set_flag( ieee_denorm,.FALSE. )
#endif
#endif
#endif 
#endif 
             STOP '*** FATAL ERROR ENCOUNTERED IN BATCH MODE'
          ENDIF

C         SET REGISTER 9 EQUAL TO 1. SO THAT YOU CAN TEST AND ABORT 
C         THE PROGRAM IF NECESSARY WITHIN A SUBROUTINE (UNUSED??)
C         PARAM(10) = 1.0
          CALL REG_SET(9,1.0, CHAR(0), IRTFLG) 
     
      RETURN
      END
