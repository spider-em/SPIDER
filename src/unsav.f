
C++*********************************************************************
C
C UNSAV.F              DOCNAM LENGTHENED DEC 9 88 al
C                      REPLACES UNSAVD ArDean Leith
C                      USED LUNDOC JUNE 99 ArDean Leith
C                      OPENDOC PARAMETERS CHANGED DEC 2000 ARDEAN LEITH
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
C PURPOSE: RETRIEVE PARAMETERS FROM A SINGLE SPECIFED KEY FOUND IN
C          A DOC. FILE.  OPENS FILE IF NECESSARY.   
C                
C UNSAV(DOCNAM,NOPENT,NDOC,IKEY,PLIST,NLIST,IRTFLG,NSS)
C
C PARAMETERS:
C	DOCNAM		NAME OF DOC FILE - CHAR. ARRAY           (SENT)
C	NOPENT          SWITCH =<0,  SUPRESSES ERROR MESSAGES    (SENT)
C	                SWITCH = 0,  FIRST TIME - OPEN DOC FILE
C			SWITCH = 1,  DOC FILE ALREADY OPENED
C	NDOC		LUN NUMBER FOR DOC FILE                  (SENT)
C	IKEY		KEY WANTED IN DOC FILE                   (SENT)
C	PLIST  		ARRAY OF RETRIEVED NUMBERS           (RETURNED)
C	NLIST  		NUMBER OF ELEMENTS IN ARRAY              (SENT)
C	IRTFLG	        0 = NO ERRORS IN DOC FILE
C			1 = ERROR IN DOC FILE 
C			2 = DOC FILE NOT FOUND               (RETURNED)
C	NSS		0 = ALL KEYS ARE SEARCHED, THEN REWOUND  (SENT)
C			1 = SEQUENTIAL MODE; FIRST MATCHING
C			    KEY IS USED, DOC IS NOT REWOUND
C			2 = SEQUENTIAL MODE; FIRST MATCHING
C			    KEY IS USED, DOC REWOUND & SEARCHED 
C                           AGAIN IF NOT FOUND ON FIRST PASS
C
C--*********************************************************************
		
        SUBROUTINE UNSAV(DOCNAM,NOPENT,NDOC,IKEY,PLIST,NLIST,IRTFLG,NSS)

	INCLUDE 'CMLIMIT.INC' 
	INCLUDE 'CMBLOCK.INC' 
	
	REAL, DIMENSION(*) ::     PLIST
        CHARACTER(LEN=*) ::       DOCNAM
	LOGICAL    ::             NEWFILE,FLAGER,TILLEND,GOBACK

        SAVE NDOCNIC

C	NSS FLAG = 0, ALL KEYS ARE SEARCHED, AT END DOC IS REWOUND
        TILLEND = (NSS .EQ. 0)

C       SET ERROR RETURN FLAG
	IRTFLG   = 1

        NOPEN  = NOPENT
        FLAGER = .TRUE.
        IF (NOPEN .LT. 0) THEN
C          NEGATIVE NOPEN SUPPRESSES ERROR MESSAGES ON RETRIEVAL
           FLAGER = .FALSE.
           NOPEN  = ABS(NOPEN)
        ENDIF

        IF (NLIST .LE. 0) THEN
	   CALL ERRT(101,'*** REGISTER LIST EMPTY IN UNSAV',IDUM)
	   RETURN

        ELSEIF (NOPEN. EQ. 0) THEN
C          MUST OPEN FILE FIRST
           IF (.NOT. FLAGER) IRTFLG = -9
           CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOC,NDOCNIC,.FALSE.,' ',
     &                  .TRUE.,.TRUE.,.TRUE.,NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           CALL  LUNDOCSAYHDR(NDOCNIC,NOUT,IRTFLG)
        ENDIF

C       RETRIEVE IKEY'S VALUES
        GOBACK = (NOPEN .NE. 0) .AND. (NSS .EQ. 2)
        CALL  LUNDOCREDDAT(NDOCNIC,IKEY,PLIST,NLIST,IGOT,
     &                     TILLEND,GOBACK,IRTFLG)

        IF (IRTFLG .EQ. 0) THEN
C          DESIRED KEY WAS FOUND OK
           IF (IGOT .LT. NLIST .AND. FLAGER) THEN
C             DID NOT FIND AS MANY REGISTERS AS WANTED FOR THIS KEY
              WRITE(NOUT,8887) NLIST,IKEY,IGOT
8887	      FORMAT
     &        (' *** WANTED:',I3,' REG. FOR KEY:',I7,' ONLY GOT:',I3)
              IRTFLG = 1
           ENDIF

        ELSE
C          DID NOT EVEN FIND THIS KEY
           IF (FLAGER) THEN
              WRITE(NOUT,8889) IKEY
8889          FORMAT(' *** KEY:',I7,' NOT FOUND')
           ENDIF
           IRTFLG = 1
        ENDIF

	IF (TILLEND .AND. NDOCNIC .GT. 0) REWIND NDOC

	RETURN
	END
