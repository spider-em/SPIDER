
C ++*******************************************************************
C
C   CONF1.F    LONG FILE NAMES              JAN 89   ARDEAN LEITH
C                                          4/30/93   MAHIEDDINE LADJADJ
C                                           9/1/93   JING SU
C              REMOVED UNUSED 'RS'        10/27/08   ARDEAN LEITH
C              TYPET                        5/7/12   ARDEAN LEITH
C              'EP TM' PROMPT              1/28/14   ARDEAN LEITH
C
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
C                                                                      *
C **********************************************************************

        SUBROUTINE CONF1(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

C       ACTUAL MAX. LENGTH OF A0 IS SET IN CMLIMIT.INC
        COMMON /IOBUF/ A0(NBUFSIZ)

        PARAMETER (NFUNC=3)
        CHARACTER(LEN=2)      :: CFUNC(NFUNC)
        CHARACTER(LEN=MAXNAM) :: FILNAM,FILTST
        CHARACTER             :: DISP
        LOGICAL               :: EX,ISOPEN
        CHARACTER(LEN=3)      :: TYPET
        CHARACTER(LEN=1)      :: NULL = CHAR(0)

C       COMMAND OPTIONS:
        DATA CFUNC /'PT', 'CI', 'EP'/

        LUN1   = 10
        LUN2   = 7
        LUN3   = 8
        LUN4   = 9
        LUN5   = 12
        MAXIM  = 0
        MAXIM1 = 0
        MAXIM2 = 0
        MAXIM3 = 0
        MAXIM4 = 0

        DO IFUNC = 1, NFUNC
          IF (FCHAR(1:2) == CFUNC(IFUNC)) THEN
            GOTO ( 2, 3, 4), IFUNC
          ENDIF
        ENDDO
        RETURN


C       OPERATION ----------------------------------------------- 'PT'
2       CONTINUE

        IF (FCHAR(4:4)  ==  'U') THEN
           CALL ERRT(101,'UNKNOWN OPERATION: PT U',NE)
           RETURN
        ENDIF

C       GET NAME FOR EXISTING OR NEW IMAGE FILE
10      CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG  ==  -1) RETURN

C       MERGE FILNAM WITH DATEXC IF NECESARY
        IF (FILNAM(1:1) .NE. '_' .AND. INDEX(FILNAM,'@')  ==  0) THEN
           CALL FILCAN(FILTST,NLET,NULL,NULL,FILNAM,DATEXC,IER)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
           FILTST = FILNAM
        ENDIF

C       FIND IF FILE EXISTS
        TYPET = 'FI'
        CALL INQUIREIF1(LUN1,FILTST,TYPET,EX,ISOPEN,LUNOP,
     &                  INLNED,IMGNUM,IRTFLG)

        DISP = 'U'
        IF (EX) DISP = 'O'

C       OPEN IMAGE FILE
	IF (FCHAR(4:4)  ==  '3')  THEN
	   IFORM  = 3
           NSLICE = 0
	ELSE
           IFORM  = 1
           NSLICE = 1
	ENDIF
        NSAM   = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUN1,DISP,IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'IMAGE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 10
        IF (.NOT. EX) THEN

C          FILL THE OUTPUT IMAGE WITH ZEROS
           A0(1:NSAM) = 0.0

           DO  I = 1, NROW*NSLICE
              CALL WRTLIN(LUN1,A0,NSAM,I)
           ENDDO

           FMIN = 0.0
           FMAX = 0.0
        ELSE
           IF (IMAMI.NE.1) 
     &        CALL NORM3(LUN1,NSAM,NROW,NSLICE,FMAXD,FMIND,AVD)
        ENDIF

C       CREATE PATTERN
        CALL PTTERN(LUN1,NSAM,NROW,NSLICE,FMAX,FMIN)
        GOTO 9000

C      OPERATION ------------------------------------------------- 'CI'
C      INPUT FILE (MEAN FILE)
3      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM1,NROW1,NSLIC1,
     &  	   MAXIM,'AVERAGE',.FALSE.,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9000

       IF (NSLIC1 .NE. 1) THEN
          CALL ERRT(101,'DOES NOT WORK ON VOLUMES',NE)
          GOTO 9000
       ENDIF

C      SECOND INPUT FILE CONTAINING THE VARIANCE
       
       CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM,NSAM2,NROW2,NSLIC2,
     &  	   MAXIM2,'VARIANCE',.FALSE.,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9000

       IF (NSLIC2 .NE. 1) THEN
          CALL ERRT(101,'DOES NOT WORK ON VOLUMES',NE)
          GOTO 9000

       ELSEIF (NSAM1 .NE. NSAM2 .OR. NROW1 .NE. NROW2) THEN
          CALL ERRT(101,'INCONSISTENT DIMENSIONS',NE)
          GOTO 9000
       ENDIF

C      OUTPUT FILE OR FILE TO RECEIVE UPPER LIMIT OF CONFIDENCE INTERVAL
       NROW3 = NROW2
       NSAM3 = NSAM2
       NSAM  = NSAM2
       NROW  = NROW2
       IFORM = 1
       CALL OPFILEC(0,.TRUE.,FILNAM,LUN3,'U',IFORM,NSAM3,NROW3,1,
     &  	   MAXIM3,'UPPER LIMIT OUTPUT',.FALSE.,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9000

C      FILE TO RECEIVE LOWER LIMIT OF CONFIDENCE INTERVAL
       NROW4 = NROW2
       NSAM4 = NSAM2
       IFORM = 1
       CALL OPFILEC(0,.TRUE.,FILNAM,LUN4,'U',IFORM,NSAM4,NROW4,1,
     &  	   MAXIM,'LOWER LIMIT OUTPUT',.FALSE.,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9000

C      CONFIDENCE INTERVAL
       CALL RDPRMI(N,IDUM,NOT_USED,'NUMBER OF FILES ADDED')

       CALL RDPRM(ALPHA,NOT_USED,'ERROR PROBABILITY IN %')

       CALL CONF(N,ALPHA,LUN1,LUN2,LUN3,LUN4,NSAM,NROW)
       GOTO 9000

C      OPERATION	------------------------------------------- 'EP'
C      SECOND INPUT FILE CONTAINING THE VARIANCE
4      IF(FCHAR(4:5) == 'TT')  THEN
C
	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &  	     MAXIM,'FIRST AVERAGE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000

	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM,NSAM,NROW,NSLICE,
     &  	     MAXIM1,'FIRST VARIANCE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000

          CALL RDPRMI(N1,IDUMP,NOT_USED,'NUMBER OF IMAGES AVERAGED')

	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN3,'O',IFORM,NSAM,NROW,NSLICE,
     &  	     MAXIM2,'SECOND AVERAGE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000

	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN4,'O',IFORM,NSAM,NROW,NSLICE,
     &  	     MAXIM3,'SECOND VARIANCE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000

          CALL RDPRMI(N2,IDUMP,NOT_USED,'NUMBER OF IMAGES AVERAGED')
									
	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN5,'U',IFORM,NSAM,NROW,NSLICE,
     &  	     MAXIM4,'SIGNIFICANCE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000
	  IF (IFORM < 1)  THEN
             CALL ERRT(101,'DOES NOT WORK ON FOURIER FILES',NE)
             GOTO 9000
          ENDIF

          CALL  TTEST(LUN1,LUN2,LUN3,LUN4,LUN5,N1,N2,
     &                NSAM,NROW,NSLICE)
          CLOSE(LUN5)

	ELSEIF(FCHAR(4:5) == 'TM')  THEN
C
	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &  	     MAXIM,'AVERAGE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000

	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM,NSAM,NROW,NSLICE,
     &  	       MAXIM1,'VARIANCE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000

          CALL RDPRMI(N1,IDUMP,NOT_USED,'NUMBER OF IMAGES AVERAGED')

	  CALL RDPRM(UM,NOT_USED,'TESTED POPULATION MEAN')

	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN5,'U',IFORM,NSAM,NROW,NSLICE,
     &  	     MAXIM2,'SIGNIFICANCE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000
	  IF (IFORM < 1)  THEN
             CALL ERRT(101,'DOES NOT WORK ON FOURIER FILES',NE)
             GOTO 9000
          ENDIF

          CALL TTEST1(LUN1,LUN2,LUN5,N1,UM,NSAM,NROW,NSLICE)
          CLOSE(LUN5)

	ELSEIF(FCHAR(4:5) == 'TP')  THEN
          CALL ERRT(101,'UNDOCUMENTED OPERATION REMOVED',NE)
          !CALL TPAIRED  July 2014 al

	ELSEIF(FCHAR(4:5) == 'MM')  THEN

          CALL ERRT(101,'UNDOCUMENTED OPERATION REMOVED',NE)
          !CALL HOTM July 2014 al

	ELSEIF(FCHAR(4:5) == 'T2')  THEN

          CALL ERRT(101,'UNDOCUMENTED OPERATION REMOVED',NE)
          !CALL HOTELLING July 2014 al

        ELSE
	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM,
     &               NSAM2,NROW2,NSLIC2,
     &  	     MAXIM,'VARIANCE',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000

          IF (NSLIC2 .NE. 1) THEN
             CALL ERRT(101,'DOES NOT WORK ON VOLUMES',NE)
             GOTO 9000
          ENDIF

          NROW3 = NROW2
          NSAM3 = NSAM2
          NSAM  = NSAM2
          NROW  = NROW2
          IFORM = 1
	  CALL OPFILEC(0,.TRUE.,FILNAM,LUN3,'U',IFORM,NSAM3,NROW3,1,
     &  	     MAXIM2,'OUTPUT',.FALSE.,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000

C         ERROR PROBABILITY
          CALL RDPRMI(N,IDUM,NOT_USED,'NUMBER OF FILES ADDED')

          CALL RDPRM(A,NOT_USED,'MAGNITUDE OF CONFIDENCE INTERVAL')

          CALL PROBER(A,N,LUN2,LUN3,NSAM,NROW,NZERO)

	  IF (NZERO .GT. 0)  WRITE(NOUT,401)  NZERO
401	  FORMAT('  WARNING ! ',I0,' NEGATIVE VALUES ENCOUNTERED')
        ENDIF
        GOTO 9000

9000    CLOSE(LUN1)
        CLOSE(LUN2)
        CLOSE(LUN3)
        CLOSE(LUN4)

        END

