C++*********************************************************************
C
C FILGEN.F          ALTERED NOV 87 FOR NEW FILE FORMAT al
C                   LONG FILENAMES ADDED DEC 88    ARDEAN LEITH
C                   USED OPFILE NOV 00             ARDEAN LEITH
C                   SGI LEAK ON INTERNAL FMT     AUB 02 ArDean Leith
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
C    FILGEN(FILNAM,LUN1)
C
C    PARAMETERS:     FILNAM    CHAR. VARIABLE FOR FILE NAME
C                    NLET      NUMBER OF CHARS. IN FILE NAME
C                    LUN1      LOGICAL UNIT NUMBER FOR OPENING FILES
C
C    PURPOSE:    LISTS FILE PARAMETERS ON A SERIES OF IMAGES
C
C--*********************************************************************

	SUBROUTINE FILGEN(FILNAM,NLET,LUN1)

        INCLUDE 'CMBLOCK.INC'

        LOGICAL        NOFIND
        CHARACTER *(*) FILNAM

        CALL GETFILENUM(FILNAM(1:NLET),IFIRST,IDIG,.FALSE.,IRTFLG)

11      IF (IMGNUM .LT. 0) THEN
C         SPECIAL CASE FOR NO DIGITS AT END OF FILE-NAME
          MAXIM = 0
          CALL OPFILEC(0,.FALSE.,FILNAM,LUN1,'Z',IFORM,NSAM1,NROW1,NDUM,
     &                   MAXIM,' ',.TRUE.,IRTFLG)

          IF (IRTFLG .EQ. 0) THEN
C            FILE FOUND. PRINT OUT FILE INFORMATION USING FILDAT
             CALL FILDAT(LUN1,NSAM1)
             CLOSE(LUN1)
          ENDIF
          RETURN
        ENDIF

C       IDIG IS NUMBER OF CONSECUTIVE DIGITS AT END OF THE FIRST FILE NAME
        LASTFI = 10**IDIG - 1
        IGO    = NLET - IDIG + 1

20      IFOUND = 0
	NOFIND = .TRUE.
        NUMNOT = 0
C       NUMNOT COUNTS HOW MANY SUCCESSIVE FILES HAVE NOT BEEN FOUND

	DO  IFILE = IFIRST,LASTFI
C         CREATE NEXT FILE NAME
          CALL INTTOCHAR(IFILE,FILNAM(IGO:NLET),NNN,IDIG)

          MAXIM = 0
          CALL OPFILEC(0,.FALSE.,FILNAM,LUN1,'Z',IFORM,NSAM1,NROW1,NDUM,
     &                   MAXIM,' ',.TRUE.,IRTFLG)

          IF (IRTFLG .EQ. 0) THEN
C            FILE FOUND. PRINT OUT FILE INFORMATION USING FILDAT
C            THIS MULTIPLE LISTING CAN NOT SET PARAMETER VALUES

             CALL FILDAT(LUN1,NSAM1)
	     NOFIND = .FALSE.
             NUMNOT = 0
             CLOSE(LUN1)

           ELSE
C            FILE NOT FOUND
             NUMNOT = NUMNOT +1
             IF (NUMNOT .GT. 10) THEN
C               STOP AFTER 10 NON-EXISTING FILES
                IFOUND = 1
                EXIT
             ENDIF
           ENDIF
        ENDDO

        IF (NOFIND) WRITE(NOUT,*) '*** NO SUCH FILES'
        WRITE(NOUT,*) ' '
        RETURN
           
999     WRITE(NOUT,*) '*** ERROR IN FILENAME OR PGM.'
        RETURN

        END
