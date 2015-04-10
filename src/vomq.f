C++*********************************************************************
C
C VOMQ.F          REWRITTEN                         JUL 03 ArDean Leith
C                 'AP SH' SUPPORT                   MAY 04 ArDean Leith
C                 MIRROR FLAG ADDED TO OUTPUT       DEC 04 ArDean Leith
C                 MIRROR FLAG BUG                   APR 06 ArDean Leith
C                 LABELS IN DOC FILE                FEB 13 ArDean Leith
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
C VOMQ
C
C PURPOSE: CREATES SELECTION DOC FILE FROM 'AP SH' OR 'AP REF' OUTPUT
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE VOMQ

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        REAL, POINTER          :: APMQBUF(:,:)
        CHARACTER(LEN=MAXNAM)  :: FILNAM,FILPAT,DOCNAM,GRPDOCFIL
        CHARACTER(LEN=80    )  :: COMMEN,COMMEN2
        CHARACTER(LEN=1)       :: NULL = CHAR(0)
        REAL                   :: DLIST(4)
        INTEGER                :: ITOTAL
        LOGICAL                :: NEWFILE

        INTEGER, PARAMETER     :: LUND1  = 80
        INTEGER, PARAMETER     :: LUND2T = 81
        INTEGER, PARAMETER     :: LUND3T = 82


        CALL RDPRM(CCTHRESH,NOT_USED,'CCC THRESHOLD')

C       'AP SH' DOC FILE COLUMNS
C       PSI,THE,PHI, MIR-REF#, EXP#, INPLANE<, SX,SY, NPROJ, DIFF, CCROT,
C             INPLANE<,SXNEW,SYNEW,MIR-NEW
C       'AP MQ' DOC FILE COLUMNS
C       MIR-REF#, CCROT, INPLANE<, SX, XY, EXP#

C       RETRIEVE  APMQ/APSH  DOC FILE CONTENTS (CLOSES LUND1)
        MAXYT = 0
        MAXXT = 0
        CALL GETDOCDAT('ALIGNMENT DOC FILE',.TRUE.,DOCNAM,LUND1,
     &                .TRUE.,MAXXT,MAXYT,APMQBUF,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        MAXREF = 0
        CALL RDPRI1S(MAXREF,NOT_USED,
     &               'NUMBER OF REFERENCES USED IN ALIGNMENT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       FIND NAME FOR OUTPUT DOC FILE 
        NMAX = 0
        CALL FILELIST(.TRUE.,NDUM,FILPAT,NLET,IDUM,NMAX,IDUM,
     &       'TEMPLATE FOR REFERENCE SELECTION DOC FILES',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       OPEN GROUP OUTPUT DOC FILE 
        CALL OPENDOC(GRPDOCFIL,.TRUE.,NLETG,LUND2T,LUND2,.TRUE.,
     &              'REFERENCE SUMMARY OUTPUT',.FALSE.,.TRUE.,
     &              .TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999


C                 123456789 123456789 123456789 123456789 
        COMMEN = 'Ref.  # OF IMAGES'               
        CALL LUNDOCPUTCOM(LUND2,COMMEN(1:20),IRTFLG)

C                 123456789 123456789 123456789 123456789 123456789
        COMMEN = 'Key   EXP IMAGE         CC     <0 IS MIRRORED'               

C                 123456789 123456789 123456789 123456789 123456789
        COMMEN2= '  MATCHES FOR REFERENCE:xxxxxx'               

C	# OF IMAGES FOR THIS GROUP
        ITOTAL = 0.0

C       LOOP OVER REF  IMAGES SETS FROM APMQ/APSH DOC FILE
	DO IREF = 1,MAXREF

           CALL FILGET(FILPAT,FILNAM,NLET,IREF,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999

C          OPEN CURRENT GROUP REF. SELECTOR OUTPUT FILE, NO OPEN ECHO
           IRTFLG = -9
           CALL OPENDOC(FILNAM,.TRUE.,NLETD,LUND3T,LUND3,.FALSE.,
     &              ' ',.FALSE.,.TRUE.,.TRUE.,NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999
        
           WRITE(COMMEN2(25:),'(I6)',ERR=999) IREF
           CALL LUNDOCPUTCOM(LUND3,COMMEN2(1:30),IRTFLG)
           CALL LUNDOCPUTCOM(LUND3,COMMEN(1:45),IRTFLG)
       
C          LOOP OVER ALL IMAGES IN APMQ/APSH  DOC FILE
           IKEY = 0
	   DO IMG = 1,MAXYT

             IREFNOW = 0
             ICOUNT  = APMQBUF(1,IMG)
             IF (ICOUNT > 0 .AND. ICOUNT < 11) THEN
C               KEY USED IN APMQBUF, APMQ STYLE INPUT

                IREFNOW = ABS(APMQBUF(2,IMG))
                CCCNOW  = APMQBUF(3,IMG) 
                IMGNOW  = APMQBUF(7,IMG)
                FMIRNOW = APMQBUF(2,IMG)
            
             ELSEIF (ICOUNT > 0) THEN
C               KEY USED IN APMQBUF, APSH STYLE INPUT
                IREFNOW = ABS(APMQBUF(5,IMG))
                CCCNOW  = APMQBUF(12,IMG) 
                IMGNOW  = APMQBUF(6,IMG)
                IF (ICOUNT >= 15) THEN
                    FMIRNOW = APMQBUF(16,IMG)
                ELSE
                    FMIRNOW = 0.0
                ENDIF
             ENDIF

             IF (IREFNOW == IREF .AND. CCCNOW > CCTHRESH) THEN
C               THIS KEY IS PRESENT IN APMQBUF
                IKEY     = IKEY + 1

C               ADD LINE TO CURRENT GROUP REF. SELECTOR OUTPUT FILE
                DLIST(1) = IMGNOW
                DLIST(2) = CCCNOW
                DLIST(3) = 1.0
                IF (FMIRNOW < 0.0) DLIST(3) = -1.0

                CALL LUNDOCWRTDAT(LUND3,IKEY,DLIST,3,IRTFLG)
             ENDIF
          ENDDO

          IF (IREF <= 1 .AND. VERBOSE) WRITE(NOUT,90)
90        FORMAT('  Reference number    No. of images')

          IF (VERBOSE) WRITE(NOUT,91) IREF,IKEY
91        FORMAT(1X,I6,8X,I8)

C         CLOSE CURRENT GROUP REF. SELECTOR OUTPUT FILE
          CLOSE(LUND3T)

          ITOTAL = ITOTAL + IKEY

C         ADD LINE TO CURRENT GROUP TOP-LEVEL OUTPUT FILE
          DLIST(1) = IKEY
          CALL LUNDOCWRTDAT(LUND2,IREF,DLIST,1,IRTFLG)

        ENDDO
         
        WRITE(NOUT,93) ITOTAL
93      FORMAT('  Total number of images above CCC threshold: ',I9)

999     IF (ASSOCIATED(APMQBUF)) DEALLOCATE(APMQBUF)
        NULLIFY(APMQBUF)

	CLOSE(LUND2T)
        CLOSE(LUND3T)

        END
