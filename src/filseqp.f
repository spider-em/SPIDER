
C++*********************************************************************
C
C    FILSEQP.FOR         GETS INPUT WITH PROMPT TO CREATE SEQUENTIAL 
C                        LONG FILENAMES
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
C     FILSEQP(FILPAT,NLET,ILIST,NMAX,NUM,PROMPT,IRTFLG)
C
C     PURPOSE:          INPUTS FILE NAME TEMPLATE AND NUMBERS FOR FILE
C                       NAME LOOP.  USUALLY USED WITH FILGET.FOR
C                       I.E.  CALL FILGET(FILPAT,FILNAM,NLET,INUM,IRTFLG)
C                         
C     PARAMETERS:       FILNAM    FILE NAME PATTERN         (RETURNED)
C                       NLET      LENGTH OF FILNAM          (RETURNED)
C                       ILIST     ARRAY FOR NUMBERS         (RETURNED)
C                       NMAX      MAX. LENGTH OF ILIST      (NEEDED)
C                                 IF ZERO ONLY GETS FILPAT NOT ILIST
C
C                       NUM       NUMBER OF VALUES IN ILIST (RETURNED)
C                       PROMPT    PROMPT                    (NEEDED)
C                       IRTFLG    ERROR FLAG; 0 IS NORMAL   (RETURNED)
C
C--*********************************************************************

        SUBROUTINE FILSEQP(FILPAT,NLET,ILIST,NMAX,NUM,PROMPT,IRTFLG)

 

	INCLUDE 'CMBLOCK.INC'

   	CHARACTER *(*)FILPAT,PROMPT
	CHARACTER *1  NULL

C       ILIST IS DIMENSIONED AS (1) HERE SO NMAX=0 IS ACCEPTED
C**	INTEGER*4     ILIST(NMAX)      ! ACTUAL SIZE
	INTEGER*4     ILIST(1)

        NULL=CHAR(0)

C       SET NORMAL ERROR RETURN FLAG

C       DO NOT CHANGE CASE OF THE LETTERS
5       IRTFLG = -999

C       READ IN FILE NAME TEMPLATE
        CALL RDPRMC(FILPAT,NLET,.TRUE.,PROMPT,NULL,IRTFLG)
        IF (IRTFLG .EQ. -1) RETURN

        IF (NLET .EQ. 3 .AND. FILPAT(NLET:NLET) .NE. '*') THEN
C           MAKE NEW STYLE TEMPLATE
            FILPAT(4:7) = '***' // NULL
            NLET = 6

        ELSE
            FILPAT(NLET+1:NLET+1) = NULL
        ENDIF

        IF (NMAX .GT. 0) THEN
C          FILL THE NUMBERS ARRAY ALSO
C          SET NUM TO NMAX FOR NUMBER OF FILES ALLOWED
           NUM = NMAX
10         CALL RDPRAI(ILIST,NMAX,NUM,0,9999999,'FILE NUMBERS',
     &         NULL,IRTFLG)
           IF (IRTFLG .EQ. -1) GOTO 5
        ENDIF

        RETURN
        END


