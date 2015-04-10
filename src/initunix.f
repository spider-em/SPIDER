
C ++********************************************************************
C
C INITUNIX.F   ABILITY TO REGISTERS ON START         JUN 00 ARDEAN LEITH
C              SET MEM REMOVED                       JAN 01 ARDEAN LEITH
C              SET MP REMOVED                        JUN 02 ARDEAN LEITH
C              SPIREOUT                              JUN 05 ARDEAN LEITH
C              REG SET REMOVED                       JAN 06 ARDEAN LEITH
C              MPI CHANGES                           OCT 08 ARDEAN LEITH
C              SET_MPI                               DEC 10 ARDEAN LEITH
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
C  INITUNIX(NUMARG,FCHART,NALPHT,CXNUM,COMMANDLINE)
C
C  PURPOSE:  RECOVERS & PROCESSES COMMAND LINE ARGUMENTS                                    *
C            SETS USE_SPIRE IN CMBLOCK.INC
C            
C  PARAMETERS:  NUMARG      NUMBER OF COMMAND LINE ARGUEMENTS (RETURNED)
C               FCHART      OPERATION ON COMMAND LINE         (RETURNED)
C               NALPHT      NO. OF CHARS IN OPERATION         (RETURNED)
C               CXNUM       RESULTS FILE NUMBER               (RETURNED)
C               COMMANDLINE SPIDER NAME BEING RUN             (RETURNED)
C                                                     
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE INITUNIX(NUMARG,FCHART,NALPHT,CXNUM,COMMANDLINE)

#ifdef SP_NT
	use dflib
#endif
        INCLUDE 'CMBLOCK.INC'

        CHARACTER(LEN=*)   :: FCHART,CXNUM,COMMANDLINE

        INTEGER FUNCTION  iargc

        CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SETS ICOMM AND MYPID

        NUMARG = 0

C       GET COMMAND LINE USING iargc
        IF (MYPID .LE. 0) NUMARG = iargc() 
C        write(6,*) mypid,'  numarg:',numarg

C       GET COMMAND BEING RUN FOR RETURNING TO CALLER --------------
        CALL GETARG(0,COMMANDLINE)

C       CHECK FOR DATA EXTENSION ON COMMAND LINE --------------------          
        FCHART = ''
        IF (NUMARG .GE. 1) THEN
C          HAS PROJECT/DATA EXTENSION ON COMMAND LINE (ARG. 1)
           CALL GETARG(1,FCHART)
           PRJEXC(1:3) = FCHART(1:3)

           IF (FCHART(4:4) .NE. '/') THEN
              DATEXC(1:3) = PRJEXC(1:3)
           ELSE
              DATEXC(1:3) = FCHART(5:7)
           ENDIF
           IF (DATEXC(1:1) .EQ. ' ' .OR. DATEXC(2:2) .EQ. ' ' .OR.
     &         DATEXC(3:3) .EQ. ' ') 
     &        STOP '*** INVALID DATA EXTENSION ON COMMAND LINE!'
        ENDIF

C       CHECK FOR FIRST OPERATION (USUALLY BATCH) ON COMMAND LINE
        IF (NUMARG .GE. 2) THEN
C          HAS FIRST OPERATION GIVEN ON COMMAND LINE (ARG. 2)
           CALL GETARG(2,FCHART)
           NALPHT = lnblnk(FCHART)

C          CHECK FOR SPIRE OPERATION  ON COMMAND LINE 
           USE_SPIRE = (FCHART(1:5) .EQ. 'SPIRE')
           IF (USE_SPIRE) THEN
              FCHART = FCHART(6:)
              NALPHT = MAX(0,NALPHT-5)
           ENDIF
        ENDIF

C       CHECK FOR RESULTS FILE VERSIONING ON COMMAND LINE (ARG. 3) ---
        CXNUM = CHAR(0)
        IF (NUMARG .GE. 3) THEN
C          RESULTS FILE VERSIONING GIVEN ON COMMAND LINE (ARG. 3)
           CALL GETARG(3,CXNUM)
           IF (CXNUM(1:1) .EQ. '-') THEN
              NUMARG = 2
              CXNUM  = CHAR(0)
           ENDIF
	ENDIF

C       REGISTER SETTING FROM COMMAND LINE (ARG. 4...) DONE IN SPIDER
        IF (NALPHT .LE. 0) NUMARG = MIN(NUMARG,1)

#ifdef USE_MPI

        CALL BCAST_MPI('INITUNIX','DATEXC',DATEXC,3,'C',ICOMM)
        CALL BCAST_MPI('INITUNIX','PRJEXC',PRJEXC,3,'C',ICOMM)

        ILEN = LEN(FCHART)
        CALL BCAST_MPI('INITUNIX','FCHART',FCHART,ILEN,'C',ICOMM)

        ILEN = LEN(CXNUM)
        CALL BCAST_MPI('INITUNIX','CXNUM',CXNUM,ILEN,'C',ICOMM)
C       write(6,*) ' cxnum:',cxnum

        ILEN = LEN(COMMANDLINE)
        CALL BCAST_MPI('INITUNIX','COMMANDLINE',COMMANDLINE,ILEN,
     &                 'C',ICOMM)
C       write(6,*) ' commandline:',commandline

        CALL BCAST_MPI('INITUNIX','NUMARG',NUMARG,1,'I',ICOMM)
        CALL BCAST_MPI('INITUNIX','NALPHT',NALPHT,1,'I',ICOMM)
#endif
        END



