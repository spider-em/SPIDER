
C++***************************************************6/12/80 1/5/81 VAX
C
C  ENDIT.F                                          REWRITTEN SEPT 97 al
C                                    DELETE CLOSED RESULT FILE JAN 03 al
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
C  ENDIT(MESG,DELETIT,RESULT)
C
C--*********************************************************************

	SUBROUTINE ENDIT(MESG,DELETIT,RESULT)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*)       :: MESG,RESULT
        LOGICAL                :: DELETIT

        LOGICAL                :: ISOPEN 
        INTEGER                :: IER 
        INTEGER                :: N,ISIZE,IRET,system,NLET
        CHARACTER (LEN=MAXNAM) :: FILOPEND,LOG
        CHARACTER (LEN=3)      :: LOGM = 'LOG'

C       PRINT COMPLETION TIME, -1 IS NOT STARTING NEW PAGE IN FILE
C       IN RESULTS FILE
	CALL PDATES(MESG,-1)

        INQUIRE(UNIT=NDAT,OPENED=ISOPEN,NAME=FILOPEND)

C       CLOSE RESULTS FILE (3) (IF NOT DIVERTED TO null)
        IF (INDEX(FILOPEND,'dev/null') == 0) CLOSE(NDAT)

	IF (DELETIT) THEN
C          DELETE RESULTS FILE (NDAT) 
           OPEN(NDAT,FILE=RESULT,STATUS='OLD',IOSTAT=IER)
           IF (IER == 0) CLOSE(NDAT,STATUS='DELETE',IOSTAT=IER)
       ENDIF

C       CLOSE LOG FILE (1), (MAY ALREADY BE CLOSED, IF IN PROCEDURE)
        INQUIRE(UNIT=1,OPENED=ISOPEN,SIZE=ISIZE,NEXTREC=N,NAME=FILOPEND)

C       SIZE AND NEXTREC ARE ALWAYS ZERO!!! COMPILER PROBLEM??
        !write(6,'(A,i5,a,L)') ' NEXTREC: ',N,    ' ISOPEN:',ISOPEN
        !write(6,'(A,i5,a,a)') ' SIZE:',    ISIZE,' NAME:',FILOPEND
        !write(6,'(A,i5)')     ' NECHO:',   NECHO
        !write(6,'(A,a)')      ' LOG:',     LOG

        IF (NECHO <= 0) THEN
C          LOG FILE IS EMPTY. DELETE IT
           IF (.NOT. ISOPEN) THEN
              CALL FILNAMANDEXT(LOGM,PRJEXC,LOG,NLET,.FALSE.,IER)
              OPEN(UNIT=1,FILE=LOG,STATUS='OLD',IOSTAT=IER)
           ENDIF
           CLOSE(UNIT=1,STATUS='DELETE',IOSTAT=IER)
        ELSE
C          LOG FILE IS NOT EMPTY. CLOSE IT
	   IF (ISOPEN) CLOSE(UNIT=1)
        ENDIF

C       STOP IS CALLED AFTER RETURN IN CALLER

	END
