
C++*********************************************************************
C
C SAVDN1.F  
C                             LONG FILENAMES      JAN   89 al
C                             OPENDOC PARAMETERS  DEC 2000 ARDEAN LEITH
C                             INCORE OPENDOC      JUL 2003 ARDEAN LEITH
C                             IRTFLG = 9          JUL 2003 ARDEAN LEITH
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
C    SAVDN1(NDOC,DOCNAM,DLIST,NLISTT,NRUN,IAP)
C
C    PURPOSE:    SUBROUTINE TO SAVE PARAMETERS IN DOCUMENT FILE WHICH
C                MAY BE ALREADY OPENED, CALLED INSIDE A PROGRAM
C
C    PARAMETERS:   NDOCT   LUN NUMBER OF FILE TO SAVE REGISTERS   (SENT)
C                  DOCNAM  NAME FOR DOC FILE 
C                              (SPIDER EXTENSION NOT NECESSARY)   (SENT)
C                  DLIST   ARRAY CONTAINING FLOATING PT. NUMBERS  (SENT)
C                                 TO BE SAVED.
C                              (FIRST NUMBER IS KEY)
C	           NLISTT  NUMBER OF ELEMENTS IN ARRAY            (SENT)
C	                     (<0 IS FLAG TO NOT ECHO OPEN INFO)
C                  NRUN    0 IF FIRST CALL (OPENS FILE), ELSE 1   (SENT)
C                  IAP     1 IF OPEN/APPEND, 0 IF OPEN/REWIND     (SENT)
C
C--*******************************************************************

	SUBROUTINE SAVDN1(NDOCT,DOCNAM,DLIST,NLISTT,NRUN,IAP)

        INCLUDE 'CMBLOCK.INC'

	REAL,DIMENSION(*) :: DLIST
	CHARACTER(LEN=*)  :: DOCNAM
	LOGICAL           :: ADDEXT,APPEND,NEWFILE                

        SAVE NDOC

	IF (NRUN .LE. 0) THEN
C         FIRST CALL OF THIS ROUTINE, OPEN DOC FILE NOW

          ADDEXT = .TRUE.
          IDOT   = INDEX(DOCNAM,'.',.TRUE.) 
          IF ((IDOT .GT. 0)  .AND.  
     &        (INDEX(DOCNAM(IDOT:),'/') .LE. 0) .AND.
     &        (INDEX(DOCNAM(IDOT:),CHAR(92)) .LE. 0)) ADDEXT = .FALSE.
   
          APPEND = (IAP .NE. 0) 
C         IF UNWANTED, DO NOT ECHO FILE OPENING INFO
          IF (NLISTT .LT. 0)  IRTFLG = -9
          CALL OPENDOC(DOCNAM,ADDEXT,NLET,NDOCT,NDOC,.FALSE.,' ',
     &                 .FALSE.,APPEND,.TRUE.,NEWFILE,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN
        ENDIF

C       PUT REGISTERS IN DOC FILE
        NLIST  = ABS(NLISTT) - 1
        IF (NLIST .LE. 0) RETURN

        IKEY = DLIST(1)

        CALL LUNDOCWRTDAT(NDOC,IKEY,DLIST(2),NLIST,IRTFLG)

        RETURN
	END

