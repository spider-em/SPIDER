
C++*************************************************************************
C
C VMS_CD.F                         
C                             NEW                 JAN 2009 ARDEAN LEITH
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
C   VMS  SOLICTS SYSTEM COMMAND AND THEN RUNS THAT COMMAND 
C
C--*******************************************************************

	SUBROUTINE VMS_CD()

	INCLUDE 'CMBLOCK.INC' 

	CHARACTER(LEN=160) :: COMMAN
        INTEGER            :: chdir
        LOGICAL            :: GETANS,STRIP
        LOGICAL            :: UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI

#ifdef USE_MPI
        INCLUDE 'mpif.h'
        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID, IERR)
#else
        MYPID = -1
#endif

        COMMAN = CHAR(0)

C       DO NOT UPPERCASE THE INPUT LINE, STRIP AFTER ;
        UPPER     = .FALSE.
        WANTSUB   = .TRUE.
        SAYPRMT   = .TRUE.
        SAYANS    = .FALSE.
        ENDATSEMI = .TRUE.
        GETANS    = .TRUE.
        STRIP     = .TRUE.

        CALL RDPR('DIRECTORY',NC,COMMAN,GETANS,
     &             UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
        IF (MYPID .LE. 0) THEN

cc         write(6,*) comman(:nc)
           IF (VERBOSE) WRITE(NOUT,*) ' '

           IRET = chdir(COMMAN(1:NC))
        ENDIF

	END
