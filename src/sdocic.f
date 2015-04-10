
C++*********************************************************************
C
C SDOCIC.F     USED LUNDOC                          JUN  99 ARDEAN LEITH
C               REMOVED ALLOCIT                     MAY  00 ARDEAN LEITH
C               LUNDOCREDALL PARAMETERS CHANGED     DEC  00 ARDEAN LEITH
C               ADDED 'SD IC' CAPABILITY            JUN  03 ARDEAN LEITH
C               IPQ(ILOC) = NLIST                   JUN  04 ARDEAN LEITH
C               IPQ(ILOC) OVERFLOW TRAP             MAY  09 ARDEAN LEITH
C               NEXTKEY SET                         MAY  09 ARDEAN LEITH
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
C   SDOCIC(SDNEW,SDCOPY)
C
C   PURPOSE:  ROUTINE DUPLICATES FUNCTION OF 'SD'
C             BUT USES IN-CORE STORAGE, WHICH MAKES SUBSEQUENT READ 
C             ACCESSES VERY QUICK.  IT SOLICITS DOC FILE NAME.  THE
C             REGISTERS ARE SPECIFED ON THE COMMAND LINE
C   
C   PARAMETERS: SDNEW     FLAG TO CREATE NEW FILE               (SENT)
C               SDCOPY    FLAG TO COPY TO DISK FILE             (SENT)
C
C   USAGE:    'SD IC NEW'      -- FIRST ACCESS: ALLOCCATE DOC. ARRAY 
C             'SD IC 11,[reg]' -- WRITE DIRECTLY INTO DOC. ARRAY
C             'SD IC COPY'     -- COPY ARRAY TO DISK BASED DOC. FILE
C
C--*********************************************************************

	SUBROUTINE SDOCIC(SDNEW,SDCOPY)

        USE DOCIC_INFO

        INCLUDE 'CMBLOCK.INC' 
C       'CMLIMIT.INC' IS AVAILABLE FROM: DOCIC_INFO

        REAL, DIMENSION(:), POINTER :: IPQ
        CHARACTER(LEN=MAXNAM)       :: DOCNAM
        LOGICAL                     :: COMOUT,SDNEW,SDCOPY
        LOGICAL                     :: NEWFILE
        CHARACTER(LEN=1)            :: NULL

        NULL = CHAR(0)

        DATA LUNDOCT/70/

        CALL FILERD(DOCNAM,NLET,NULL,'FILE OR ARRAY NAME~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       SEE IF THIS FILE IS ALREADY OPEN IN-CORE
        CALL ISDOCINCORE(DOCNAM(1:NLET),NIC,NEWNIC,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (SDNEW) THEN
C          CREATE NEW EMPTY IN-CORE DOC. FILE ----------------------

           IF (NIC .GT. 0) THEN
C             THIS INCORE FILE ALREADY EXISTS
              WRITE(NOUT,90) DOCNAM(1:NLET)
90            FORMAT(' *** ',A,' : ALREADY IN-CORE',/)
              RETURN
    
           ELSEIF (NEWNIC .EQ. 0) THEN
C             NO SPACE AVAILABLE IN DOC NAME LIST
              WRITE(NOUT,91) 
91            FORMAT(' *** NO SPACE AVAILABLE IN DOC. NAME LIST, ',
     &               ' CLOSE ANOTHER INCORE FILE FIRST')
              CALL ERRT(100,'SDOCIC',IDUM)
              RETURN
           ENDIF

           CALL RDPRIS(MAXREG,MAXY,NOT_USED,
     &                 'NUMBER OF REGISTERS & KEYS ALLOWED',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          ALLOCATE AN RUNTIME SRRAY FOR DOC. FILE CONTENTS
           MEMWANT = (MAXREG + 1) * MAXY
           ALLOCATE(IPQ(MEMWANT),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
               CALL ERRT(102,'INCORE DOC. ALLOCATION FAILED',MEMWANT)
               RETURN
           ENDIF

C          ZERO THE WHOLE ARRAY
           IPQ = 0

C          KEEP ARRAY POINTER IN LUNDOC
           LOCDOC(NEWNIC)%IPT => IPQ
           OLDNAM(NEWNIC)     = DOCNAM(1:NLET)
           NLETOLDNAM(NEWNIC) = NLET
           NUMKEYS(NEWNIC)    = MAXY
           NUMCOLS(NEWNIC)    = MAXREG + 1
           NEXTKEY(NEWNIC)    = 1
           NICLAST            = NEWNIC
 	   RETURN
 
        ELSEIF (SDCOPY) THEN
C          COPY INCORE DOC FILE TO PHYSICAL OR IN-CORE --------------

C          OPEN OUTPUT DOC FILE
           CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUNDOCT,LUNDOC,.TRUE.,
     &             'OUTPUT DOCUMENT',.FALSE.,.FALSE.,.TRUE.,
     &             NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0)RETURN
           IF (LUNDOC .LT. 0) THEN
              CALL ERRT(101,'CAN NOT COPY TO INCORE ARRAY',NE)
              RETURN
           ENDIF

C          INCORE POINTER TO THIS DOC. FILE CONTENTS 
           IPQ => LOCDOC(NIC)%IPT

C          GET ARRAY SIZE SAVED WHEN FILLED
           MAXX   = NUMCOLS(NIC)  
           NUMREG = MAXX - 1
           MAXY   = NUMKEYS(NIC)
           ILOC   = 1
           DO  IKEY=1,MAXY
              IF (IPQ(ILOC) .GT. 0) THEN
C               KEY IS IN USE
                CALL LUNDOCWRTDAT(LUNDOC,IKEY,IPQ(ILOC+1),NUMREG,IRTFLG)
              ENDIF
              ILOC = ILOC + MAXX
           ENDDO
           CLOSE(LUNDOCT)
           RETURN     ! ---------------------------------------------
        ENDIF

C       PARSE REGISTER LINE, CHECK FOR ',' OR 'X'

        IX = 6
        CALL REG_DOC_PARSE(FCHAR(IX:),COMOUT,IKEY,NLIST,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       REGISTER LIST OK, CHECK INCORE ARRAY SIZE WHEN CREATED
        MAXX   = NUMCOLS(NIC)  
        MAXY   = NUMKEYS(NIC)
        IF (IKEY .GT. MAXY) THEN
           CALL ERRT(102,'OVERFLOW, MAX. KEY ALLOWED',MAXY)
           RETURN
        ELSE IF (NLIST .GT. (MAXX - 1)) THEN
           CALL ERRT(102,'MAX. NUMBER OF REGISTERS',MAXREG)
           RETURN
        ELSEIF (COMOUT) THEN
           CALL ERRT(100,'IN-CORE COMMENT KEYS NOT ALLOWED',NE)
           RETURN
        ENDIF

C       RETRIEVE THE INCORE POINTER  
        IPQ => LOCDOC(NIC)%IPT

C       WANT TO SAVE REGISTERS IN INCORE ARRAY 
        ILOC      = (IKEY-1) * MAXX + 1
        IPQ(ILOC) = NLIST

C       RETRIEVE DATA FROM REGISTER(S) LISTED IN NSEL INTO: IPQ(ILOC+1)
        CALL REG_GET_NSELA(NLIST,IPQ(ILOC + 1),IRTFLG)

        IF (.NOT. SILENT) WRITE(NOUT,*) ' '

	RETURN
	END


