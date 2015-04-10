
C++*********************************************************************
C
C SAVD.F                      ALTERED          JAN 18 1988 ARDEAN LEITH
C                             LONG FILE NAMES     JAN   88 ARDEAN LEITH
C                             OPENDOC PARAMETERS  DEC 2000 ARDEAN LEITH
C                             NDOCIN              JUL 2003 ARDEAN LEITH
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
C  SAVD(NDOCIN,DLIST,NLIST,IRTFLG)
C
C  PURPOSE: SAVE PARAMETERS IN DOCUMENT FILE, CALLED INSIDE A PROGRAM, 
C           SOLICITS FILENAME & OPEN DOC FILE ON FIRST CALL.  IF FILE
C           EXISTS IT WILL OPEN FOR APPEND, NOT REPLACE.
C
C  PARAMETERS:
C       NDOCIN   LUN NUMBER OF FILE TO SAVE REGISTERS             (SENT)
C	DLIST    ARRAY CONTAINING FLOATING PT. NUMBERS            (SENT)
C                      TO BE SAVED.
C	NLIST    NUMBER OF ELEMENTS IN ARRAY                      (SENT)
c       IRTFLG   ERROR FLAG (0 IS NORMAL)                     (RETURNED)
C
C  NOTE:         HAS SAVDC ENTRY POINT!!!!
C
C--*********************************************************************

	SUBROUTINE SAVD(NDOCIN,DLIST,NLIST,IRTFLG)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        REAL,DIMENSION(*)      :: DLIST
        CHARACTER(LEN=MAXNAM)  :: DOCNAM
	LOGICAL                :: NEWFILE,APPEND,OPENED

        SAVE          OPENED,APPEND,NDOC

	DATA          OPENED/.FALSE./,APPEND/.TRUE./

	IF (.NOT. OPENED) THEN
C          GET NAME AND OPEN THE DOC FILE, SETS NDOC
           CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOCIN,NDOC,.TRUE.,
     &                  'DOCUMENT',
     &                 .FALSE.,APPEND,.TRUE.,NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           OPENED = .TRUE.
        ENDIF

        IKEY = DLIST(1)
	CALL LUNDOCWRTDAT(NDOC,IKEY,DLIST(2),NLIST-1,IRTFLG)

	RETURN

	ENTRY SAVDC
	APPEND = .TRUE.
	OPENED = .FALSE.
        RETURN

	ENTRY SAVD_OPENED
	OPENED = .TRUE.
	RETURN

	ENTRY SAVD_NOAPPEND
	APPEND = .FALSE.
	RETURN

	END

