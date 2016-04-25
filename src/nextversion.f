 
C++*********************************************************************
C
C NEXTVERSION.F                                                  3/2/94
C                       ADDED CXNUMB                AUG 00 ARDEAN LEITH
C                       NEXTRESULTS IN BIN DIR      FEB 09 ARDEAN LEITH
C                       OUTPUT FORMATTING           AUG 14 ARDEAN LEITH
C                       OUTPUT FORMATTING,DELETE    JAN 16 ARDEAN LEITH
C
C **********************************************************************
C *  AUTHOR: MAHIEDDINE LADJADJ                                            *
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
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
C   NEXTVERSION(FILIN,FILOUT,LUNT,CXNUMB)
C
C   PARAMETERS:     FILIN  FILE NAME FOR WHICH THE HIGHEST
C                           VERSION NUMBER HAS TO BE FOUND.      (SENT)
C                   FILOUT  FILE NAME WILL HAVE A VERSION
C                                   NUMBER APPENDED TO IT.)      (RET.)
C                   LUNT    READ UNIT                            (SENT)
C                   CXNUMB  VERSION NUMBER  FOR RES. FILE    (SENT/RET)                                (RET.)
C                                       
C       CALLS A SCRIPT CALLED Nextresults THAT WILL
C       DO A CHECK OF THE DIRECTORY FILES AND FIND THE NEXT VERSION
C       NUMBER FOR A FILE. Nextresults WILL WRITE THE NAME OF THE FILE
C       WITH THE CORRECT VERSION TO A FILE CALLED SPIDER_JUNK.TMP.
C       THIS ROUTINE READS SPIDER_JUNK.TMP AND DELETES IT.
C
C--*********************************************************************

        SUBROUTINE NEXTVERSION(FILIN,FILOUT,LUNT,CXNUMB,SAYIT)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
  
        CHARACTER (LEN=*)      :: FILIN, FILOUT
        INTEGER                :: LUNT
        CHARACTER (LEN=4)      :: CXNUMB
        LOGICAL                :: SAYIT

        CHARACTER (LEN=1)      :: NULL = CHAR(0)
        CHARACTER (LEN=280)    :: SCRIPT
        CHARACTER (LEN=9)      :: FIND_DIR = 'SPBIN_DIR'
	LOGICAL                :: EX
        CHARACTER (LEN=MAXNAM) :: DIR_NAME, FILNAM
        INTEGER                :: NCHAR,NLET,INUMB,IDUM,IRTFLG,NLEN
        INTEGER                :: IERR,ILOC
        
        INTEGER                :: lnblnkn


        FILOUT = FILIN
        NCHAR  = lnblnkn(FILOUT)

        IF (CXNUMB .NE. NULL) THEN
C          RESULTS EXTENSION SET BY CALLER
           NLET   = lnblnkn(CXNUMB)
           FILOUT = FILIN(1:NCHAR) // '.' // CXNUMB(1:NLET) // ' '
           IF (SAYIT) THEN
              NCHAR  = lnblnkn(FILOUT)
              WRITE(NOUT,97) FILOUT(1:NCHAR)
           ENDIF

           RETURN
        ENDIF
 


C       SET DEFAULT VALUE USING LAST 3 DIGITS OF SYSTEM CLOCK
        CALL SYSTEM_CLOCK(INUMB,IDUM,IDUM)
        INUMB                   = MOD(INUMB,1000)
        FILOUT(NCHAR+1:NCHAR+1) = '.' 
        CALL INTTOCHAR(INUMB,FILOUT(NCHAR+2:NCHAR+4),NLET,3)

#ifndef SP_IBMSP3
C       OMIT THIS FOR NOW ON PARALLEL ARCH. al Oct 00

        CALL MYGETENV(FIND_DIR,DIR_NAME,NLEN,
     &                'DIRECTORY_OF_BINARY_FILES', IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

	FILNAM = DIR_NAME(1:NLEN) // 'Nextresults'// NULL
	INQUIRE(FILE=FILNAM,EXIST=EX,IOSTAT=IERR)
        IF (IERR .NE. 0) GOTO 999

        IF (.NOT. EX) THEN
	  WRITE(NOUT,93) FILNAM
93	  FORMAT(//,'*** FILE NOT FOUND: ',A)
	  GOTO 999
        END IF
        
C       CALL THE SCRIPT Nextresults TO FIND THE NEXT VERSION NUMBER OF
C       FILE FILENAME.

        SCRIPT = '$' // FIND_DIR // '/Nextresults ' // FILIN //
     &               ' > SPIDER_JUNK.TMP' // NULL

        CALL system(SCRIPT)

C       READ THE FILE FILENAME.NBR' FOUND AND STORED IN SPIDER_JUNK.TMP
        OPEN(LUNT, FILE='SPIDER_JUNK.TMP',STATUS='OLD',IOSTAT=IERR)
        IF (IERR .NE. 0) THEN
           WRITE(NOUT,*) 
     &     '*** Could not open SPIDER_JUNK.TMP for results file version')
           GOTO 999
        ENDIF

        READ(LUNT,12,IOSTAT=IERR) FILOUT
12      FORMAT(A)
        IF (IERR .NE. 0) THEN
           WRITE(NOUT,*) 
     &     '*** Could not read SPIDER_JUNK.TMP for results file version')
           GOTO 999
        ENDIF

C       DELETE THE TEMP FILE.
        CLOSE(LUNT,STATUS='DELETE')

C       DELETE THE TEMP FILE.
cc      SCRIPT = 'rm -f SPIDER_JUNK.TMP' // NULL
cc      CALL system(SCRIPT)  NO LONGER NEEDED AFTER STATUS=DELETE

#endif

C       GET CXNUMB FROM FILOUT
999     ILOC  = INDEX(FILOUT,'.',.TRUE.)
        READ(FILOUT(ILOC+1:),80,IOSTAT=IERR) INUMB
80      FORMAT(I8)

        IF (IERR == 0) THEN
            CALL INTTOCHAR(INUMB,CXNUMB(1:3),NLET,3)
        ENDIF

        IF (SAYIT) THEN
           NCHAR  = lnblnkn(FILOUT)
           WRITE(NOUT,97) FILOUT(1:NCHAR)
        ENDIF
97      FORMAT('  Results file: ',A)

        END        

