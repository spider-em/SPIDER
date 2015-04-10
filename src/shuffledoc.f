C++*********************************************************************
C
C SHUFFLEDOC     M.RADERMACHER 5/95
C                REWRITTEN                     MAY 97   ARDEAN LEITH                        
C                OPENDOC PARAMETERS CHANGED    DEC 2000 ARDEAN LEITH
C                INCORE OPENDOC                JUL 03   ARDEAN LEITH
C                PROMPTS                       JAN 13   ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C SHUFFLEDOC(MAXDIM)
C
C PURPOSE:              SHUFFLES A DOCUMENT FILE
C                  
C PARAMETERS:           MAXDIM   LENGTH OF COMMON BUFFER          (SENT)           
C
C--*********************************************************************

        SUBROUTINE SHUFFLEDOC(MAXDIM)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	COMMON        PLIST(16),BTUF(1)	

        CHARACTER(LEN=MAXNAM) :: DOCF1
	INTEGER               :: P1,P2,P3

        DATA          NDOCT/23/
	
C       INPUT DOCUMENT FILENAME, ~9 ALLOWS DIFFERENT EXTENSION
        CALL FILERD(DOCF1,NLET,DATEXC,'SOURCE DOCUMENT~9',IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

        CALL OPENDOC(DOCF1,.FALSE.,NLET,NDOCT,NDOC,.FALSE.,' ',
     &             .TRUE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       ECHO FIRST HEADER FROM FILE
        CALL LUNDOCSAYHDR(NDOC,NOUT,IRTFLG)

C       FIND MAXKEYT & MAXREGT BY READING DOC FILE
        CALL LUNDOCINFO(NDOC,MAXKEYT,MAXREGT,KEYUSED,.TRUE.,IRTFLG)
        CLOSE(NDOCT)
        IF (IRTFLG .NE. 0) RETURN

        CALL RDPRMI(KEY1,KEY2,NOT_USED,'FIRST & LAST KEY NUMBERS')

        WRITE(NOUT,100)
100     FORMAT(' THE OUTPUT DOCFILE WILL HAVE THE FORMAT:'/
     &         ' NEWKEY, 4 OLDKEY, + 3 COLUMNS')

       	CALL RDPRMI(ICOL0,IDUM,NOT_USED,
     &              'INPUT REGISTER #  FOR OUTPUT REG. 2')

	CALL RDPRMI(ICOL1,ICOL2,NOT_USED,
     &              'INPUT REGISTER #s FOR OUTPUT REG. 3 AND 4')

C	THE SIZE OF MAXKEY IS IN REALITY MMKEY (ADD 10 FOR GOOD MEASURE)
	MMKEY    = MAX0(KEY1,KEY2) + 10
        IREGSRED = MAXREGT + 1
	P1       = 1
	P2       = P1 + IREGSRED * MAXKEYT 
	P3       = P2 + MMKEY 
	MEMNEED  = P3 + 4 * MMKEY
	IF (MEMNEED .GE. (MAXDIM -20)) THEN
	   WRITE(NOUT,414) MEMNEED, MAXDIM
414	   FORMAT('*** COMMON BLOCK MEM NEEDED: ',I7,
     &                 ' MEM AVAILABLE ONLY: ',I7)
           CALL ERRT(100,'SHUFFLEDOC',NE)
           RETURN
	ENDIF
	
        NREG  = ICOL2 
        IF (ICOL1 .GT. NREG)  NREG = ICOL1
        IF (ICOL0 .GT. NREG)  NREG  = ICOL0
        IRN   = KEY2 - KEY1 + 1
        IR1   = 0
        IR2   =     IRN
        IR3   = 2 * IRN

	ICALL = 0
        ILINE = 0
        DO I=KEY1,KEY2
           CALL UNSDAL(DOCF1,ICALL,NDOCT,I,PLIST,NREG,
     &                BTUF(P1),MAXKEYT,IREGSRED,NKEYMAX,LERR)
           IF (LERR .NE. 0) GOTO 11
           ILINE                = ILINE + 1
           BTUF(P2+ILINE-1)     = I
           BTUF(P3+IR1+ILINE-1) = PLIST(ICOL0)
           BTUF(P3+IR2+ILINE-1) = PLIST(ICOL1)
           BTUF(P3+IR3+ILINE-1) = PLIST(ICOL2)
        ENDDO

11      CLOSE(NDOCT)

C       DO THE SHUFFLING, BY RANDOM INTERCHANGE

        DO  I=0,ILINE-1
C         CREATE RANDOM IVAL IN RANGE 0...ILINE-1
	  CALL  RANDOM_NUMBER(OUT)
	  IVAL              = OUT*FLOAT(ILINE-1)

          B1                = BTUF(P2+    IVAL) 
          B2                = BTUF(P3+IR1+IVAL)
          B3                = BTUF(P3+IR2+IVAL)
          B4                = BTUF(P3+IR3+IVAL)

          BTUF(P2    +IVAL) = BTUF(P2    +I) 
          BTUF(P3+IR1+IVAL) = BTUF(P3+IR1+I)
          BTUF(P3+IR2+IVAL) = BTUF(P3+IR2+I)
          BTUF(P3+IR3+IVAL) = BTUF(P3+IR3+I)

          BTUF(P2    +I)    = B1 
          BTUF(P3+IR1+I)    = B2
          BTUF(P3+IR2+I)    = B3
          BTUF(P3+IR3+I)    = B4
        ENDDO

        DO I=0,ILINE-1
           PLIST(1) = I + 1
           PLIST(2) = BTUF(P2    +I)
           PLIST(3) = BTUF(P3+IR1+I)
           PLIST(4) = BTUF(P3+IR2+I)
           PLIST(5) = BTUF(P3+IR3+I)
           CALL SAVD(NDOCT,PLIST,5,IRTFLG)
        ENDDO

        CALL SAVDC
        CLOSE (NDOCT)

	END
