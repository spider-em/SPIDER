C++*********************************************************************
C
C VORA.F                ROTATE DOCUMENT FILE     07/12/93                         07/12/93  PP
C                       ADDED 'VO RAS'           MAY 01    ARDEAN LEITH
C                       -K DOES NOT PRINT        JUL 01    ARDEAN LEITH
C                       MISSING KEY BUG          OCT 2012  ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  VORA()
C
C  PURPOSE:  ROTATE DOCUMENT FILE  
C            MULTIPLIES INVERTED INPUT ANGLES BY ANGLES FROM THE 
C            ANG DOC FILE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE VORA()

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INCLUDE 'F90ALLOC.INC'
        REAL,  POINTER         :: DOCBUF(:,:)

        REAL                   :: FI1(3),FI2(3),FIO(4)
        CHARACTER (LEN=MAXNAM) :: DOCNAM,DOCOUT
        CHARACTER (LEN=1)      :: NULL = CHAR(0)
        CHARACTER(LEN=40)      :: COMMEN
   
        INTEGER                :: NDOCINT,MAXX,MAXY,IRTFLG,NOT_USED
        INTEGER                :: NICDOCOUT,ICOMM,MYPID,MPIERR,IANG
        INTEGER                :: NLET,NDOC,NLIST,IKEY,NEWKEY,ICOUNT
        INTEGER                :: NDOUT
        REAL                   :: FANG,FVAL,DLIST(3)
        LOGICAL                :: ASKNAM,GETSIZ,NEWFILE,ADDEXT,ISOLD
        LOGICAL                :: WRTCOM,APPEND

        INTEGER,PARAMETER      :: NDOCIN  = 81
        INTEGER,PARAMETER      :: NDOCOUT = 82


        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       OPEN EXISTING DOC FILE
C       MAXX IS 1 + NUM OF REGISTERS SINCE DOCBUF CONTAINS KEY ALSO
        MAXX    = 0
        MAXY    = 0
        NDOCINT = NDOCIN
        ASKNAM  = .TRUE.
        GETSIZ  = .TRUE.

        CALL GETDOCDAT('ANGULAR DOCUMENT',ASKNAM,DOCNAM,
     &                  NDOCINT,GETSIZ,MAXX, MAXY,DOCBUF,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN


        IF (FCHAR(4:6) == 'RAS') THEN

            CALL  RDPRM3S(FI1(1),FI1(2),FI1(3),NOT_USED,
     &                   'ROTATION ANGLES - PSI, THETA & PHI',IRTFLG)
            IF (IRTFLG .NE. 0)  RETURN

11          CALL  RDPRM2S(FANG,FVAL,NOT_USED,
     &                   'ANGLE NUMBER(e.g. PSI IS 1) & VALUE',IRTFLG)
            IF (IRTFLG .NE. 0)  RETURN

            IANG = FANG
            IF (IANG < 0 .OR. IANG > 3) THEN
               CALL ERRT(102,'ANGLE NUMBER MUST BE 0...3',IANG)
               GOTO 11
            ENDIF
        ELSE

            FI1(1) = 0.0
            FI1(2) = 0.0
            FI1(3) = -999999

            CALL RDPRM3S(FI1(1),FI1(2),FI1(3),NOT_USED,
     &             'ROTATION ANGLES - PSI, THETA & PHI',IRTFLG)
            IF (FI1(3) == -999999) THEN

               FI1(3) = 0
               CALL RDPRM1S(FI1(3),NOT_USED,
     &             'ROTATION ANGLE - PHI',IRTFLG)
            ENDIF
            IANG = 0
        ENDIF

        ADDEXT  = .TRUE.
        ASKNAM  = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .FALSE.
        WRTCOM  = .TRUE.


        CALL OPENDOC(DOCOUT,ADDEXT,NLET,
     &               NDOCOUT,NICDOCOUT,ASKNAM,'OUTPUT DOCUMENT',
     &               ISOLD,APPEND,WRTCOM,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C                 123456789 123456789 123456789 1234567890
        COMMEN = '           PSI,       THETA,         PHI'
        CALL LUNDOCPUTCOM(NDOCOUT,COMMEN(1:40),IRTFLG)

C       DOCBUF HAS ICOUNT IN FIRST COL
        NLIST  = MAXX

        NEWKEY = 0
        DO IKEY = 1,MAXY
           ICOUNT = DOCBUF(1,IKEY)        ! NO. OF REG. FOR THIS KEY 

           IF (ICOUNT  >= 3) THEN
C             GOT VALID DOC FILE DATA LINE WITH 3 ANGLES

C             CONVERT ANGLES
              CALL CALD(FI1, DOCBUF(2,IKEY), DLIST)

C             CAN SET ONE ANGLE TO SPECIFIED VALUE
              IF (IANG > 0) DLIST(IANG) = FVAL

C             PUSH DLIST INTO OUTPUT DOC. FILE
              CALL LUNDOCWRTDAT(NICDOCOUT,IKEY,DLIST,3,IRTFLG)

           ENDIF
        ENDDO

C       DEALLOCATE DOC. FILE MEMORY
9995    IF (ASSOCIATED(DOCBUF)) DEALLOCATE(DOCBUF)
        IF (MYPID <= 0) CLOSE(NDOUT)
        IF (MYPID <= 0) CLOSE(NDOC)

        END
