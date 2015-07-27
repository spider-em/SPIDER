
C++*************************************************************************
C
C VMS.F                         
C                             REMOVED FILNAMSUB   APR 2001 ARDEAN LEITH
C                             MULTILINE           SEP 2003 ARDEAN LEITH
C                             RDPR PARAMETERS     04/14/05 ARDEAN LEITH
C                             [_d] --> xd         03/24/06 ARDEAN LEITH
C                             NULL AT END BUG     01/06/10 ARDEAN LEITH
C                             NULL AT END BUG     08/29/12 ARDEAN LEITH
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
C   VMS  SOLICTS SYSTEM COMMAND AND THEN RUNS THAT COMMAND 
C
C--*******************************************************************

	SUBROUTINE VMS(MULTILINE)

        IMPLICIT NONE
	INCLUDE 'CMBLOCK.INC' 

        LOGICAL             :: MULTILINE

	CHARACTER(LEN=160)  :: COMLIN
	CHARACTER(LEN=1600) :: COMMAN
        INTEGER             :: system,lnblnkn
        LOGICAL             :: GETANS,STRIP
        LOGICAL             :: UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI

        INTEGER             :: NC,IRTFLG,NCM,NDUM,I,NCT,IRET

        INTEGER             :: ICOMM,MYPID,MPIERR

        CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SETS ICOMM AND MYPID

        COMMAN    = CHAR(0)

C       DO NOT UPPERCASE THE INPUT LINE, DO NOT STRIP AFTER ;
        UPPER     = .FALSE.
        WANTSUB   = .TRUE.
        SAYPRMT   = .TRUE.
        SAYANS    = .FALSE.
        ENDATSEMI = .FALSE.
        GETANS    = .TRUE.
        STRIP     = .FALSE.

        CALL RDPR('SYSTEM COMMAND',NC,COMMAN,GETANS,
     &             UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
        IF (MULTILINE) THEN
C          DO NOT UPPERCASE THE INPUT LINE, DO NOT STRIP AFTER ;
10         SAYPRMT   = .FALSE.
           CALL RDPR('',NCM,COMLIN,GETANS,
     &             UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (NCM > 1 .OR. COMLIN(1:1) .NE. '.') THEN
C              CONCATENATE WITH COMMAN
               IF ((NC + NCM) > 1600) THEN
                  CALL ERRT(101,'COMMAND LIMITED TO 1600 CHAR.',NDUM)
                  RETURN
               ENDIF
               COMMAN = COMMAN(1:NC) // COMLIN(1:NCM)
               NC     = NC + NCM
               GOTO 10
           ENDIF
        ENDIF                

cc      write(6,*) comman(:nc)

        IF (NC < 160) COMMAN(NC+1:NC+1) = CHAR(0)
        NC = LNBLNKN(COMMAN)

        CALL DEBRAKXREG(COMMAN,NC)
        NC = LNBLNKN(COMMAN)

        IF (MYPID <= 0) THEN
           IF (NC <= 80) THEN
              WRITE(NOUT,90)COMMAN(1:NC)
              IF (NDAT .NE. NOUT) WRITE(NDAT,90)COMMAN(1:NC)
   90         FORMAT(/,2X,A)
           ELSE
              DO I=1,NC,80
                 NCT = lnblnkn(COMMAN(I:I+79)) -1
                 WRITE(NOUT,91) COMMAN(I:I+NCT-1)
                 IF (NDAT .NE. NOUT) WRITE(NDAT,91)COMMAN(I:I++NCT-1)
   91            FORMAT(' ',A)
              ENDDO
           ENDIF
           IF (VERBOSE) WRITE(NOUT,*) ' '

           IRET = system(COMMAN(1:NC))
        ENDIF

	END

C      *********************** DEBRAKXREG ********************************

       SUBROUTINE DEBRAKXREG(CINPUT,NCHAR)

       IMPLICIT NONE

       CHARACTER(LEN=*)   :: CINPUT
       INTEGER            :: NCHAR

       CHARACTER(LEN=161) :: CSUB
       LOGICAL            :: ISDIGI
       INTEGER            :: I,J,NDIG,IRTFLG

C        CONVERT NEW: [name] FORMAT to OLD x11 REGISTER FORMAT
         I     = 1
         DO WHILE (I < (NCHAR - 3))
            IF (CINPUT(I:I+1) == '[_') THEN
C              START OF [_
               J = I + 2
               DO WHILE (J <= NCHAR)
                  IF (ISDIGI(CINPUT(J:J))) THEN
                     J = J + 1
                     CYCLE 
                  ENDIF
                  EXIT
               ENDDO

               IF (CINPUT(J:J) .NE. ']') THEN
                  I = I + 1
                  CYCLE                 ! NOT A [_ddd]
               ENDIF

               NDIG = J - I -2 
               IF (NDIG < 1) THEN
                  I = I + 1
                  CYCLE                 ! NO dd in:  [_]
               ENDIF

               CINPUT(I:I) = 'X'
               CSUB = 'X' // CINPUT(I+2:J)
               CALL SUBCHAR(CSUB(1:NDIG+1),CINPUT,I,J,NCHAR,IRTFLG)
               I = J -1
            ENDIF
            I = I + 1
         ENDDO

         END
