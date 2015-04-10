
C++*********************************************************************
C
C UNSDAL.F          ADAPTED FOR EXTENDED FILE NAMES NOV 88 ArDean Leith
C                   ALTERED FOR RUNTIME USE SEPT 96 ArDean Leith
C                   USED LUNDOC JUNE 99 ArDean Leith
C                   OPENDOC PARAMETERS CHANGED DEC 2000 ARDEAN LEITH
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
C  UNSDAL(DOCNAM,ICALL,NDOC,IKEY,PLIST,NLIST,DBUF,MAXKEY,MAXREG,
C          NKEY,IRTFLG)
C
C  PURPOSE:  RECOVERS ARRAY OF REGISTERS FROM DOCUMENT FILES.  IF IKEY
C            IS <0 RECOVERS THE LINE OF REGISTERS FROM THE COMMENT
C            LINE HAVING THAT COMMENT KEY ONLY. 
C            CALLED WITHIN ROUTINES THAT USE A DEDICATED DOC. FILE
C
C  PARAMETERS:
C         DOCNAM        DOC. FILE NAME INCLUDING EXTENSION        (SENT)
C         ICALL         FLAG SHOWING DOC FILE IS IN DBUF NOW (SENT/RET.)                 
C                       ICALL = 1 (GET LINE FROM OPEN ARRAY)
C                       ICALL = 0 (MUST OPEN & READ FILE) 
C         NDOC          LOGICAL UNIT FOR DOC FILE                 (SENT)
C         IKEY          KEY FOR LINE OF REG. RETURNED IN PLIST    
C                       IF ZERO DOES NOT FIND ANY KEY             (SENT)
C         PLIST         ARRAY TO HOLD OUTPUT FROM KEY IKEY    (RETURNED)
C         NLIST         NUMBER OF REGISTERS RETRIEVED ON A        (SENT)
C                       LINE (1ST VALUE RETRIEVED ON EACH
C                       LINE IS KEY NO. NOT A REGISTER)
C                       MAXIMUM IS 6 (JUNE 99)          
C         DBUF          BUFFER TO RETRIEVE ARRAY FROM DOC.
C                       FILE (CLEARED TO ZEROS AT START)      (RETURNED)
C         MAXKEY        NUMBER OF HIGHEST LINE THAT CAN BE        (SENT) 
C                       RETRIEVED IN ARRAY DBUF 
C         MAXREG        ONE PLUS MAX. NUMBER OF REGISTERS PER        
C                       LINE THAT CAN BE RETRIEVED IN DBUF.
C                       (FIRST POSITION ON LINE CONTAINS KEY)     (SENT)
C         NKEY          NUMBER OF HIGHEST KEY FOUND IN FILE   (RETURNED)
C         IRTFLG        ERROR FLAG (O IS NORMAL RETURN)      (SENT/RET.)
C
C     
C  TYPICAL DOC FILE LINES:
C        KEY #REGS/LINE    VALUES ........
C  COL: 123456789 123456789 123456789 123456789 123456789 1234565789
C          1 4   20.070000   17.379999   17.379999   17.379999
C       999994   20.070000   17.379999   17.379999   17.379999
C         -1 5   21.070000   12.379999   12.379999   16.379999
C        ; COMMENT LINE (PREVIOUS LINE IS A CONTINUATION FOR KEY 99999)
C        : 1 4   20.070000   17.379999   17.379999   17.379999
C        ; COMMENT LINE (PREVIOUS LINE IS A COMMENT KEY LINE)
C
C  WHEN RETRIEVED THE LINE FOR KEY 1 IS PLACED IN LINE 1  OF DBUF AND
C  THE FIRST VALUE ON LINE ONE IS THE KEY NUMBER: 1 AND THE SECOND VALUE
C  OF LINE ONE IN DBUF IS: 20.07
C
C  THE LINE OF DBUF CORRESPONDING TO IKEY IS RETURNED IN PLIST
C  IF IKEY IS NOT FOUND AN ERROR MESSAGE IS PRINTED BUT DBUF
C  IS STILL RETURNED OK.
C
C  NOTE: AS OF 6/22/96, WE ARE ALLOWING KEYS UP TO 99999. I.E USING THE
C	  FIRST COLUMN.     ML
C  NOTE: AS OF 6/17/99, WE ARE ALLOWING KEYS UP TO 999999. AL
C
C--*******************************************************************

        SUBROUTINE UNSDAL(DOCNAM,ICALL,NDOC,IKEY,PLIST,NLIST,
     &                    DBUF,MAXKEY,MAXREG,NKEY,IRTFLG)

	INCLUDE 'CMBLOCK.INC' 

	DIMENSION       PLIST(*)
        DIMENSION       DBUF(MAXREG,MAXKEY)

	CHARACTER *(*)  DOCNAM
        CHARACTER *80   RECLIN
	LOGICAL ::      NEWFILE,WARNIT
     
C       NEGATIVE VALUE OF IRTFLG SUPPRESSES TERMINAL OUTPUT OF TITLE
	NPR    = IRTFLG
        WARNIT = .TRUE.

C       SET ERROR RETURN FLAG
	IRTFLG = 1

	IF (NLIST .LE. 0) THEN
	   WRITE(NOUT,*) '*** NUMBER OF REGISTERS NOT SPECIFIED'
	   RETURN

        ELSEIF (NLIST .GT. MAXREG) THEN
           WRITE(NOUT,101) (MAXREG-1)
101        FORMAT(' *** UNABLE TO UNSAVE DOCUMENT FILE'/
     &            ' *** CURRENTLY ONLY ',I3,' REGISTERS ALLOWED')
           RETURN
        ENDIF

        IF (IKEY .LT. 0) THEN
C          DESIRE A COMMENTED KEY

C          OPEN THE DOC FILE USING EXTENSION FROM DOCNAM
           CALL OPENDOC(DOCNAM,.FALSE.,NLET,NDOC,NDOCNIC,.FALSE.,' ',
     &               .TRUE.,.TRUE.,.TRUE.,NEWFILE,NPR)
           IF (NPR .NE. 0) RETURN
           NGOT = NLIST
           CALL LUNDOCGETCOM(NDOCNIC,IKEY,PLIST,NGOT,.FALSE.,IRTFLG)
           CLOSE(NDOC)
           RETURN
        ENDIF

C       ICALL=0 MEANS THAT CURRENTLY, NO CORE IMAGE OF DOCUMENT EXISTS.
C       OPEN FILE, DUMP DOCUMENT INTO CORE, PICK REGISTERS SELECTED,
C       SWITCH ICALL TO 1, CLOSE DOCUMENT FILE, AND RETURN.
C       ICALL=1 MEANS THAT CORE IMAGE IS AVAILABLE, AND ACCESS IS
C       AUTHORIZED BY CALLING ROUTINE WHICH DID THE COMPARISON BETWEEN
C       OLD FILE AND NEW FILE NAME.

        IF (ICALL .GT. 0) GOTO 591

C       -----------------------------------------------------------

C       OPEN THE DOC FILE USING EXTENSION FROM DOCNAM
        CALL OPENDOC(DOCNAM,.FALSE.,NLET,NDOC,NDOCNIC,.FALSE.,' ',
     &               .TRUE.,.TRUE.,.TRUE.,NEWFILE,NPR)
        IF (NPR .NE. 0) RETURN

C       ECHO FIRST HEADER FROM FILE
        CALL LUNDOCSAYHDR(NDOCNIC,NOUT,IRTFLG)
 
C       CLEAR DBUF RETURNED ANSWER BUFFER
        DO  I=1, MAXREG
           DO J = 1,MAXKEY
              DBUF(I,J) = 0.0
           ENDDO
        ENDDO

        NKEY   = 0
        NMAX   = MAXREG - 1
C -----------------------------------------------------------------

510     CALL LUNDOCREDNXT(NDOCNIC,IKEYT,PLIST,NMAX,IGO,ICOUNT,IRTFLG)
        IF (IRTFLG .EQ. 1) THEN
C          ERROR READING THIS LINE OF FILE, IGNORE THE ERROR
C          READ NEXT LINE OF DOC FILE
	   GOTO 510
        
        ELSEIF (IKEYT .GT. MAXKEY) THEN
C          KEY THAT WILL NOT FIT IN DBUF SENDS ERROR MSG.
           IF (WARNIT) THEN 
              WRITE(NOUT,9901) MAXKEY
9901          FORMAT(' KEYS GREATER THAN: ',I7,' NOT RETRIEVED')
              WARNIT = .FALSE.
           ENDIF

        ELSEIF (IRTFLG .NE. 2) THEN
C          KEY FITS IN DBUF OK, PUT KEY IN FIRST COL OF DBUF
           DBUF(1,IKEYT) = IKEYT
C          PUT VALUES IN REMAINING COLS. OF THIS LINE OF DBUF
           DO J = 2,ICOUNT+1
              DBUF(J,IKEYT) = PLIST(J-1)
           ENDDO
           IF (NKEY .LT. IKEYT) NKEY = IKEYT

C          READ NEXT LINE OF DOC FILE
	   GOTO 510
        ENDIF

C ---------------------------------------------------------------

C       END OF DOCUMENT FILE FOUND. SWITCH ICALL TO 1.
	ICALL = 1
        CLOSE(NDOC)

C ---------------------------------------------------------------

591     IF (IKEY .GT. 0) THEN
           IF (IKEY .GT. MAXKEY .OR. INT(DBUF(1,IKEY) + 0.5) .EQ. 0)THEN
	      WRITE(NOUT,8889) IKEY
8889	      FORMAT(' *** KEY:',I7,'  NOT FOUND')
              RETURN
	   ENDIF

	   DO  K = 1,NLIST
              PLIST(K) = DBUF(K+1,IKEY)
           ENDDO
        ENDIF

        IRTFLG = 0

	RETURN
	END


