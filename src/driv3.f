
C++*********************************************************************
C
C    DRIV3                  CHANGED NDOC        JULY  2000 ArDean Leith
C                           ADDED 'SY DOC'      APRIL 2005 J. LEBARRON       
C                           ADDED 'UD NEXT'     FEB   2007 ArDean Leith       
C                           ADDED 'UD FIND'     JUN   2008 ArDean Leith       
C                           ADDED 'UD MAX'      MAY   2009 ArDean Leith       
C                           ADDED 'UD MAX'      MAY   2009 ArDean Leith       
C                           ADDED 'SD H'        FEB   2010 ArDean Leith       
C                           REMOVED 'SD S'      OCT   2010 ArDean Leith       
C                           ADDED 'SY' GROEL    NOV   2011 ArDean Leith       
C                           REMOVED 'LD *'      MAY   2013 ArDean Leith       
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C   DRIV3(MAXDIM)
C
C   PURPOSE:   MAIN DRIVER FOR ROUTINES REMOVED FROM DRIVER IN MAR 93
C              CONTAINS ROUTINES ACCESSING DOC FILES 
C
C  PARAMETERS: MAXDIM     MAX LENGTH OF COMMON BUFFER
C
C  HANDLES:   'SD', 'UD', 'LD', 'SY'
c
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE DRIV3(MAXDIM)
          
	INCLUDE 'CMLIMIT.INC' 
	INCLUDE 'CMBLOCK.INC' 

        CHARACTER(LEN=MAXNAM)         :: DOCNAM
        CHARACTER(LEN=2)              :: CSYM
        CHARACTER(LEN=1)              :: NULL = CHAR(0)
        LOGICAL                       :: SDNEW,SDCOPY,ENDIT
        LOGICAL                       :: NEWFILE,HEADER

C       SINCE "UD" LEAVES DOC FILE OPEN DO NOT USE NDOC FOR OTHER USES!
        INTEGER, PARAMETER            :: NDOC    = 4
        INTEGER, PARAMETER            :: NDOC2   = 70
        INTEGER, PARAMETER            :: NDOCOUT = 72

C       CARRY OUT THE OPERATION
        SELECT CASE(FCHAR(1:2))

        CASE('SD')
C          SAVE DOC --------------------------------------------- SD **

100  	   IF (FCHAR(4:5) == 'SH') THEN
C             RANDOMLY SHUFFLE LINES IN DOCFILE ------------ SD SHUFFLE
              CALL SHUFFLEDOC(MAXDIM)

           ELSEIF (FCHAR(4:4) == 'S') THEN
C             SORT DOC FILE -------------------------------------- SD S
              !!CALL SORTDOC(MAXDIM)
              CALL ERRT(101,
     &          'OPERATION OBSOLETE, USE <DOC SORT> INSTEAD',IDUM)

           ELSEIF (FCHAR(4:4) == 'C') THEN
C             SAVE DOC CLUSTER FILE ------------------------------ SD C
C             TRANSFER COORDINATES FROM CLUSTER FILE TO DOCUMENT
              CALL RDCLUS

           ELSEIF (FCHAR(4:4) == 'N') THEN
C             SAVE DOC NONLINEAR MAPPING ----------------------- SD NLM
C             TRANSFER COORDINATES FROM NLM 2D FILE TO DOCUMENT
              CALL ERRT(101,'OPERATION NO LONGER SUPPORTED',IDUM)

           ELSEIF (FCHAR(4:5) == 'IC') THEN
C             SAVE TO INCORE DOC FILE -------------------------- SD IC *
              SDNEW  = (FCHAR(7:7) == 'N') 
              SDCOPY = (FCHAR(7:7) == 'C')  
              ENDIT  = (FCHAR(6:6) == 'E' .OR. 
     &                  FCHAR(7:7) == 'E') 
	      IF (ENDIT) THEN
                  CALL UDOCIC(NDUM,ENDIT)
              ELSE
                 CALL SDOCIC(SDNEW,SDCOPY) 
              ENDIF

           ELSE
C             SAVE REGISTERS TO DOC FILE -------------------------- SD
              CALL FILERD(DOCNAM,NLETD,DATEXC,'DOCUMENT',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 5000
              HEADER = (FCHAR(4:4) .NE. 'H')
              CALL SAVDOCQ(DOCNAM(1:NLETD),NLETT,HEADER,IRTFLG) 
           ENDIF


        CASE('UD')
C          UNSAVE DOCUMENT ------------------------------------- UD ???

200	   IF (FCHAR(4:4) == 'I') THEN
C             UNSAVE DOCUMENT - IN CORE ------------------------ UD IC?
              ENDIT  = (FCHAR(6:6) == 'E' .OR. FCHAR(7:7) == 'E') 
	      CALL UDOCIC(NDOC2,ENDIT)

  	   ELSEIF (FCHAR(4:6) == 'F E'   .OR. 
     &             FCHAR(4:5) == 'FE'    .OR.
     &             FCHAR(4:8) == 'FINDE'  .OR.
     &             FCHAR(4:9) == 'FIND E') THEN
C             UNSAVE DOCUMENT FIND - END USE ---------------- UD FIND E
	      CALL UDOCIC(NDOC2,.TRUE.)

  	   ELSEIF (FCHAR(4:7) == 'NXTE'   .OR. 
     &             FCHAR(4:8) == 'NEXTE'  .OR.
     &             FCHAR(4:8) == 'NXT E'  .OR.
     &             FCHAR(4:9) == 'NEXT E') THEN
C             UNSAVE DOCUMENT SERIAL - END USE -------------- UD NEXT E
	      CALL UDOCIC(NDOC2,.TRUE.)

  	   ELSEIF (FCHAR(4:4) == 'F') THEN
C             UNSAVE DOCUMENT FIND - IN CORE ------------------ UD FIND
	      CALL UDOCIC(NDOC2,.FALSE.)

  	   ELSEIF (FCHAR(4:7) == 'NEXT' .OR.
     &             FCHAR(4:6) == 'NXT') THEN
C             UNSAVE DOCUMENT SELECT - IN CORE ---------------- UD NEXT
	      CALL UDOCIC(NDOC2,.FALSE.)

           ELSEIF (FCHAR(4:4) == 'N') THEN
C             UNSAVE DOCUMENT KEYCOUNT -------------------------- UD N 
	      CALL FILERD(DOCNAM,NLETD,DATEXC,'DOCUMENT',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 5000
	      CALL RDDOCN(DOCNAM(1:NLETD),NDOC2,.FALSE.)

           ELSEIF (FCHAR(4:4) == 'M') THEN
C             UNSAVE MAX/MIN IN REG COL -------------------------- UD M 
	      CALL FILERD(DOCNAM,NLETD,DATEXC,'DOCUMENT',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 5000
	      CALL RDDOCN(DOCNAM(1:NLETD),NDOC2,.TRUE.)

           ELSE
C             UNSAVE DOCUMENT - REGULAR --------------------------- UD
              CALL UDOC(FCHAR,NDOC)
           ENDIF


        CASE('LD')
C          LIST DOCUMENT FILE ------------------------------------- LD
           CALL ERRT(101,'OBSOLETE BUGGY OPERATION',NDUM)
           GOTO 5000

C          CALL FILERD(DOCNAM,NLETD,DATEXC,'DOCUMENT',IRTFLG)
C	   IF (IRTFLG .NE. 0) GOTO 5000
C	   LATER TO BE SWITCHED BY COPT
C	   ISEQ = 0  
C	   CALL LISTDC(FCHAR,DOCNAM(1:NLETD),NFIL,NDOC2,ISEQ)


        CASE('SY') 
C          SYMMETRY DOCUMENT FILE --------------------------------- SY
           CALL FILERD(DOCNAM,NLET,NULL,'SYMMETRY DOCUMENT',
     &                 IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 5000

           CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOCOUT,NICDOCOUT,.FALSE.,
     &                  ' ',.FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 5000

           CALL RDPRMC(CSYM,NLET,.TRUE.,
     &               'SYMMETRY TYPE? (C/CI/T/O/I/G)',NULL,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CLOSE(NDOCOUT)
              GOTO 5000
           ENDIF

           IF (CSYM(1:1) == 'C') THEN
              CALL RDPRI1S(IFOLD,NOT_USED,
     &                     'ROTATIONAL SYMMETRY FOLDNESS',IRTFLG)
              IF (IRTFLG .NE. 0)  GOTO 5000
           ENDIF

           CALL SYMANG(CSYM,IFOLD,NICDOCOUT,IRTFLG)
           CLOSE(NDOCOUT)

        END SELECT

5000    RETURN
	END

