C++*********************************************************************
C                                                                      *
C  INQUIREIF.F                NEW ROUTINE          SEP 97 ARDEAN LEITH *
C                             READPQ ADDED         AUG 99 ARDEAN LEITH *
C                             USED REG_            AUG 00 ARDEAN LEITH *
C                             ACCEPTS *            JAN 05 ARDEAN LEITH *
C                             'IQ DI'              JAN 10 ARDEAN LEITH *
C                             ACCEPTS EXTENSION    JUN 11 ARDEAN LEITH *
C                             EXTENSION ../NAM BUG AUG 11 ARDEAN LEITH *
C                             LESS VERBOSE         MAY 14 ARDEAN LEITH *
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
C   INQUIREIF    DETERMINES IF A FILE/DIR EXISTS
C
C--*********************************************************************

        SUBROUTINE INQUIREIF(TYPET)

        INCLUDE 'CMLIMIT.INC' 
        INCLUDE 'CMBLOCK.INC' 
 
        CHARACTER (LEN=*)      :: TYPET
        CHARACTER (LEN=MAXNAM) :: FILNAM
        CHARACTER (LEN=1)      :: NULL
        LOGICAL                :: EX,ISOPEN
        LOGICAL, PARAMETER     :: FROMBACK = .TRUE.

        CALL SET_MPI(ICOMM,MYPID,MPIERR)

        NULL   = CHAR(0)
        EX     = .FALSE.
        ISOPEN = .FALSE.

        IF (TYPET .EQ. 'DIR') THEN
           CALL FILERD(FILNAM,NLET,NULL,
     &         'QUERY EXISTENCE OF DIRECTORY~',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 1000

        ELSE
           CALL FILERD(FILNAM,NLET,NULL,'QUERY EXISTENCE OF~9',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 1000

C          MERGE FILNAM WITH DATEXC IF NECESARY
           ILOCAT  = INDEX(FILNAM,'@')
           ILOCDOT = INDEX(FILNAM,'.')
           !write(6,*) 'filnam:',filnam(1:20)

           IF ((FILNAM(1:1) .NE. '_') .AND. 
     &         (ILOCAT  == 0)      ) THEN

C              SEE IF HAS EXTENSION ALREADY
               ISLASH = INDEX(FILNAM(1:NLET),'/',BACK=FROMBACK)
               IDOT   = INDEX(FILNAM(1:NLET),'.',BACK=FROMBACK)

               !write(6,*) 'islash,idot:',islash,idot

               IF ((ISLASH <= 0 .AND. 
     &              IDOT   <= 0) .OR.
     &             (IDOT   < ISLASH)) THEN

C                MERGE FILNAM WITH DATEXC 
                 CALL FILCAN(FILNAM,NLET,NULL,NULL,FILNAM,DATEXC,IRTFLG)
                 IF (IRTFLG .NE. 0) RETURN
              ENDIF
           ENDIF
        ENDIF

C       FIND IF FILE/DIR EXISTS
        CALL INQUIREIF1(33,FILNAM,TYPET,EX,ISOPEN,LUNOP,
     &                  INLNED,IMGNUM,IRTFLG)

C       SET REGISTERS (IF NECESSARY)
1000    VAL1 = 0.0
        IF (EX) VAL1 = 1.0

        VAL2 = 0.0
        IF (ISOPEN) VAL2 = 1.0

        CALL REG_SET_NSEL(1,2,VAL1,VAL2, 0.0, 0.0, 0.0,IRTFLG)

        IF (MYPID <= 0) THEN
           IF (TYPET == 'DIR') THEN
              IF (EX) THEN
                 WRITE(NOUT,*)' DIRECTORY EXISTS: ',FILNAM(1:NLET)
              ELSE  
                 WRITE(NOUT,*)' DIRECTORY DOES NOT EXIST: ',
     &                        FILNAM(1:NLET)
              ENDIF
           ELSE
              IF (EX) THEN
                 WRITE(NOUT,*)' FILE EXISTS: ',FILNAM(1:NLET)
              ELSE  
                 WRITE(NOUT,*)' FILE DOES NOT EXIST: ',FILNAM(1:NLET)
              ENDIF
           ENDIF
           IF (VERBOSE) WRITE(NOUT,*) ' '
        ENDIF

        END





