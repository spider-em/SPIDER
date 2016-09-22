
C ++********************************************************************
C
C  NEW                                           MAR 2014 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  DOCSEPARATE
C                                                                 
C  PURPOSE: USE VALUE FOUND IN A SPECIFED DOC FILE REG. AS AN INDEX TO
C           SEPARATE A DOCUMENT FILE INTO MULTIPLE DOCUMENT FILES.
C   
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE DOCSEPARATE()

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC' 

      CHARACTER(LEN=80)        :: PROMPT
      CHARACTER(LEN=160)       :: MSG
      CHARACTER(LEN=MAXNAM)    :: DOCNAM,DOCOUT,FILPAT
      CHARACTER(LEN=1)         :: NULL = CHAR(0)

      INTEGER, PARAMETER       :: NDOC     = 80
      INTEGER, PARAMETER       :: NDOCOUT  = 81
      INTEGER                  :: NLET,IRTFLG,NOT_USED,ICOUNT,IGOT
      INTEGER                  :: ICOMM, MYPID, MPIERR,  NDOCT,  NLETP 
      INTEGER                  :: IGRPCOL, NUMFIL, NLETM, LASTGRP, IKEY 
      INTEGER                  :: NLIST, IGRP, NDOCOUTT,IHIGRP,NGRP
      LOGICAL                  :: GETTEMPLATE 

      INTEGER                  :: ILIST(NIMAXPLUS)

      INTEGER, PARAMETER       :: NMAXDL = 80
      REAL                     :: DLIST(NMAXDL)
      REAL                     :: UNUSED,FIGOT

      LOGICAL                  :: ADDEXT,GETNAME,PUTMSG
      LOGICAL                  :: ISOLD,APPEND,MESSAGE,NEWFILE

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C     OPEN INPUT DOC FILE
      ADDEXT  = .TRUE.
      GETNAME = .TRUE.
      ISOLD   = .TRUE.
      APPEND  = .FALSE.
      MESSAGE = .TRUE.

      CALL OPENDOC(DOCNAM,ADDEXT,NLET,NDOC,NDOCT,GETNAME,
     &            'INPUT DOCUMENT',ISOLD,APPEND,MESSAGE,
     &             NEWFILE,IRTFLG)
      IF (IRTFLG == -1) RETURN


C               123456789 123456789 123456789 123456789 123456789 123456789
      PROMPT = 'REGISTER COL. USED FOR SEPARATION (0 IS KEY)' 
      NLETP  = 44
      CALL RDPRI1S(IGRPCOL,NOT_USED,PROMPT(1:NLETP),IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999
      IF (IGRPCOL < 0) THEN
         CALL ERRT(102,'BAD REGISTER',IGRPCOL)
         RETURN
      ENDIF

C     GET TEMPLATE FOR OUTPUT DOC FILE
C               123456789 123456789 123456789 123456789 1234567890123456789
      PROMPT = 'OUTPUT DOCUMENT FILE TEMPLATE' 

      GETTEMPLATE = .TRUE.
      CALL FILELIST(GETTEMPLATE,NDOCOUT,FILPAT,NLETP,
     &              ILIST,0,NUMFIL,PROMPT(1:29),IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      CALL RDPRMC(MSG,NLETM,.TRUE.,'OPTIONAL DOC FILE HEADER',
     &            NULL,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999
      PUTMSG = (NLETM > 0)

      IHIGRP  = -1           ! HIGHEST GROUP NUMBER
      LASTGRP = -1           ! LAST GROUP IN USE
      NGRP    = 0
      INUMBR  = 0            ! DIMENSION NIMAXPLUS FROM CMLIMIT

      ADDEXT  = .TRUE.
      GETNAME = .FALSE.
      ISOLD   = .FALSE.
      APPEND  = .TRUE.
      MESSAGE = .TRUE.
      IRTFLG  = -8         ! NO IC USE

      DO
C        READ NEXT LINE FROM DOC FILE
         CALL LUNDOCREDNXT(NDOC,IKEY,DLIST,NMAXDL,UNUSED,NLIST,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9998    ! 2=EOF
         IF (NLIST <= 0) CYCLE

         IGRP = DLIST(IGRPCOL)     ! TRUNCATES
         IF (IGRP < 0) THEN
            CALL ERRT(102,'ILLEGAL INDEX:',IGRP)
            GOTO 9999
         ENDIF
     
         IF (IGRP .NE. LASTGRP) THEN
C           NEW GROUP, OPEN NEW (POSSIBLY EXISTING) DOC FILE

            CLOSE(NDOCOUT)

            CALL FILGET(FILPAT,DOCOUT,NLETP,IGRP,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
        
            CALL OPENDOC(DOCOUT,ADDEXT,NLET,NDOCOUT,NDOCOUTT,GETNAME,
     &            '  ',ISOLD,APPEND,MESSAGE,
     &             NEWFILE,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

            IF (NEWFILE) THEN
               INUMBR(IGRP) = 0    ! SETS KEY FOR DOC FILE USAGE
               IF (PUTMSG) THEN
                  CALL LUNDOCPUTCOM(NDOCOUTT,MSG(1:NLETM),IRTFLG)
               ENDIF
            ENDIF
            LASTGRP = IGRP
            IHIGRP  = MAX(IHIGRP,IGRP)   ! HIGHEST GROUP ENCOUNTERD
            NGRP    = NGRP + 1
            
         ENDIF

C        COUNT NUMBER OF LINES IN EACH DOC FILE
         INUMBR(IGRP) = INUMBR(IGRP) + 1

C        PUSH LINE INTO THIS DOC. FILE
         CALL LUNDOCWRTDAT(NDOCOUTT,IKEY,DLIST,NLIST,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
      ENDDO


9998  CALL REG_SET_NSEL(1,1,FLOAT(NGRP),FLOAT(IHIGRP),0,0,0,IRTFLG)
      CLOSE (NDOC)
      CLOSE(NDOCOUT)

C     OPEN SUMMARY DOC FILE
      ADDEXT  = .TRUE.
      GETNAME = .TRUE.
      ISOLD   = .FALSE.
      APPEND  = .TRUE.
      MESSAGE = .TRUE.

      CALL OPENDOC(DOCNAM,ADDEXT,NLET,NDOC,NDOCT,GETNAME,
     &            'SUMMARY DOCUMENT',ISOLD,APPEND,MESSAGE,
     &             NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      IF (DOCNAM(1:1) == '*' .AND. NLET == 1) RETURN

      CALL LUNDOCPUTCOM(NDOCT,'         VALUE,  OCCURANCES',IRTFLG)
      IKEY = 0

      DO IGRP = 1,NGRP
         IGOT = INUMBR(IGRP) 
         IF (IGOT > 0) THEN

C           PUSH LINE INTO SUMMARY DOC. FILE 
            IKEY     = IKEY + 1
            DLIST(1) = IGRP          
            DLIST(2) = IGOT          
            CALL LUNDOCWRTDAT(NDOCT,IKEY,DLIST,2,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
         ENDIF
      ENDDO

9999  CLOSE (NDOC)
      CLOSE(NDOCOUT)

      END



