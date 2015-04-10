
C++************************************************************* 5/15/97
C
C MOD1.F                  USER SUPPLIED PROGRAMS.
C                         ADDED FLUSH RESULTS      NOV 00 ARDEAN LEITH
C                         ADDED 'MY MP'            JUN 01 ARDEAN LEITH
C                         OPFILEC                  FEB  03 ARDEAN LEITH
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
C  MOD1(MAXDIM)
C
C  PURPOSE:      A PLACE FOR USERS TO ADD SUBROUTINES TO SPIDER
C
C  PARAMETERS:   MAXDIM   LENGTH OF COMMON BLOCK MEMORY          (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE MOD1(MAXDIM)

	INCLUDE 'CMLIMIT.INC'
	INCLUDE 'CMBLOCK.INC'

        CHARACTER *2   FUNC

	DATA FUNC/'MY'/

C       CHECK TO BE SURE THAT OPERATION IS HANDLED BY MOD1
        IFUNC = 0
        IF (INDEX(FUNC,FCHAR(1:2)) .GT. 0) 
     &      IFUNC = (INDEX(FUNC,FCHAR(1:2)) / 3) + 1
        IF (IFUNC .LE. 0) RETURN

        IF (FCHAR(4:6) == 'ONE') THEN

C          USER SUPPLIED OPERATION  'MY ONE'
	   CALL  MYMODS(MAXDIM)

        ELSE IF (FCHAR(4:5) == 'FL') THEN

C          USER SUPPLIED OPERATION  'MY FL' FLUSHES RESULTS FILE
           CALL  FLUSHRESULTS

        ELSE IF (FCHAR(4:6) == 'MP') THEN

C          USER SUPPLIED OPERATION  "MY MP" TESTS PARALLEL OPS
5          IDELAY = 1
           IPROCS = 16
           CALL RDPRIS(IPROCS,IDELAY,NOT_USED,
     &              'NUMBER OF PROCESSORS AND DELAY INTERVAL',IRTFLG)
           IF (IRTFLG == -1) GOTO 5
           IF (IRTFLG .NE. 0) RETURN
           CALL  MYMPCHECK(IPROCS,IDELAY)

        ELSE IF (FCHAR(4:6) == 'QQ') THEN

C          USER SUPPLIED OPERATION  "TESTS QUADRI"
           CALL  QUADCHECK()

	ELSE
C          ERROR  UNKNOWN SUBOPTION TO "MY"
	   CALL ERRT(101,'UNKNOWN OPTION',NE)
	ENDIF

        RETURN
	END



        SUBROUTINE MYMPCHECK(IPROCS,IDELAY)

        CHARACTER(LEN=5) :: DATE

        DATE = 'date'//CHAR(0)

c$omp   parallel do private(i)
        DO I=1,IPROCS
            IRET   = system(DATE(1:5))
#if defined(SP_IBMSP3)
            CALL sleep_(IDELAY)
#else
            CALL SLEEP(IDELAY)
#endif
         ENDDO

         END


        SUBROUTINE QUADCHECK()

	INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FILNAM
        REAL, ALLOCATABLE     :: FDATA(:)

        LUNIN = 10
        MAXIM = 0

C       OPEN INPUT FILE, NO FOURIER INPUT ALLOWED 
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIN,'O',IFORM,NSAM,NROW,NSLICE,
     &              MAXIM,'INPUT',.FALSE.,IRTFLG)

        MEMTOT = NSAM  * NROW
        ALLOCATE(FDATA(MEMTOT))

        ILOC = 1
        DO IREC = 1, NROW
           CALL REDLIN(LUNIN,FDATA(ILOC),NSAM,IREC)
           ILOC = ILOC + NSAM
        ENDDO

        XX = NSAM + 1
        YY = NROW + 1
        DO I = 1,  MEMTOT
           IF (XX .GT. NSAM) XX = 1
           IF (YY .GT. NROW) YY = 1
           SUM = SUM + QUADRI(XX, YY, NSAM, NROW, FDATA)
           XX = XX + 0.03
           YY = YY + 0.02
         ENDDO

         WRITE (6,*) 'SUM: ',SUM
         IF (ALLOCATED(FDATA))  DEALLOCATE(FDATA)

         END
