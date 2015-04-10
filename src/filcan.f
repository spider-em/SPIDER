
C++*********************************************************************
C
C    FILCAN.FOR   ADAPTED FROM FILCON.FOR          OCT 88 al
C                 LONG FILNAMES
C                 MAXNAM                           JUL 14 ARDEAN LEITH
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
C
C   FILCAN(FILNAM,NF,DEV,DIR,MASTER,EXTEN,IER)
C
C   PURPOSE: CONSTRUCTS FILE NAME FROM CHARACTER STRINGS THAT CONTAIN
C            1) DEVICE ,2) DIRECTORY NAME, 3) MASTER FILE NAME AND 
C            4) FILE EXTENSION.
C            THIS SUBROUTINE IS INSTALLATION-DEPENDENT. IT IS WRITTEN 
C            FOR THE FILE-NAMING CONVENTIONS OF THE VAX/VMS OPERATING
C            SYSTEM.
C
C            THE FILE NAME HAS THE FOLLOWING FORMAT:
C	         <dev>:[<dir>]<master>.<exten>
C
C  PARAMETERS:  FILNAM   COMPLETE FILE NAME           (RETURNED)
C               NLET     LENGTH OF COMPLETE FILE NAME (RETURNED)
C               DEV      DEVICE NAME                  (RECEIVED)
C               DIR      DIRECTORY NAME               (RECEIVED)
C               MASTER   MASTER FILE NAME             (RECEIVED)
C               EXTEN    EXTENSION NAME               (RECEIVED)
C               IER      ERROR FLAG (0 IS NORMAL)
C
C **********************************************************************

	SUBROUTINE FILCAN(FILNAM,NLET,DEV,DIR,MASTER,EXTEN,IER)

 
        INCLUDE 'CMLIMIT.INC'

	COMMON/UNITS/LUN,NIN,NOUT

        CHARACTER *(*)  FILNAM,DEV,DIR,MASTER,EXTEN
        CHARACTER * 1   NULL

        CHARACTER(LEN=MAXNAM) :: MASTMP

        NULL = CHAR(0)
    
C       SAVE MASTER IN CASE FILNAM IS SAME AS MASTER
C       NULL OUT THE COMPLETE FILE NAME
#ifdef __osf__
        NP = INDEX(MASTER,NULL) - 1
        IF (NP .LT. 0) NP = LEN(MASTER)
	MASTMP = MASTER(1:NP) // NULL
	NLET = LEN(FILNAM)
	FILNAM(1:1) = NULL
#else
	MASTMP = MASTER // NULL
	FILNAM = NULL
#endif

	IER  = 0
C       ZERO THE RETURNED FILE NAME LENGTH
	NLET = 0      

        INDNUL = INDEX(DEV,NULL)
        IF (INDNUL .GT. 1 .OR. 
     &     (INDNUL .EQ. 0  .AND. LEN(DEV) .GT. 1)) THEN
C          DEVICE PROVIDED
           CALL LONGER(FILNAM,NLET,DEV,IRTFLG)
           ILOC = 1
           IF (IRTFLG .EQ. 1) GOTO 799
           CALL LONGER(FILNAM,NLET,':',IRTFLG)
           ILOC = 2
           IF (IRTFLG .EQ. 1) GOTO 799
        ENDIF

        INDNUL = INDEX(DIR,NULL)
        IF (INDNUL .GT. 1 .OR. 
     &     (INDNUL .EQ. 0  .AND. LEN(DIR) .GT. 1)) THEN
C          DIRECTORY PROVIDED
           CALL LONGER(FILNAM,NLET,'[',IRTFLG)
           ILOC = 3
           IF (IRTFLG .EQ. 1) GOTO 799
           CALL LONGER(FILNAM,NLET,DIR,IRTFLG)
           ILOC = 4
           IF (IRTFLG .EQ. 1) GOTO 799
           CALL LONGER(FILNAM,NLET,']',IRTFLG)
           ILOC = 5
           IF (IRTFLG .EQ. 1) GOTO 799
        ENDIF

        CALL LONGER(FILNAM,NLET,MASTMP,IRTFLG)
        ILOC = 6
        IF (IRTFLG .EQ. 1) GOTO 799

        INDNUL = INDEX(EXTEN,NULL)
        IF (INDNUL .GT. 1 .OR. 
     &     (INDNUL .EQ. 0  .AND. LEN(EXTEN) .GT. 1)) THEN
C          EXTENSION PROVIDED IN CALL

C          PUT DOT FOR EXTENSION IN FILENAME
           CALL LONGER(FILNAM,NLET,'.',IRTFLG)
           ILOC = 7
           IF (IRTFLG .EQ. 1) GOTO 799

C          PUT EXTENSION IN FILENAME
           CALL LONGER(FILNAM,NLET,EXTEN,IRTFLG)
           ILOC = 8
           IF (IRTFLG .EQ. 1) GOTO 799
        ENDIF

C       NORMAL RETURN
	RETURN  

799	WRITE(NOUT,51) ILOC,MASTMP
51	FORMAT(' *** ',I1,' CAN NOT CONSTRUCT FILENAME FROM: ',A)
	IER = 1
	RETURN

	END

