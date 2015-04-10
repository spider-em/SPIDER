
C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
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
C                                                                      *
C  MRGETINFO                                                                    *
C                                                                      *
C  PURPOSE: GETS INFORMATION NEEDED TO ALIGN IMAGES.                                                           *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C INPUT:
C     NOUT= OUTPUT DEVICE # (SCREEN)
C     NDAT= FIXED DEVICE # (DRIVE)
C OUTPUT:
C     XYPTS(LV,LS,3)= ARRAY HOLDING COORDS OF MARKERS IN EACH PROJECTION
C     IREF= INDEX OF REFERENCE IMAGE IN NUMBER, XYPTS, ETC.
C     NUMBER(LV)= ARRAY WITH INDEXES MATCHING THOSE OF XYPTS, ANGLE,
C           NUMPTS, PTACTIVE, AND SCALE HOLDING FILE NUMBERS
C     ANGLE(3,LV)= EULER ANGLES FOR A VIEW (PSI,THETA,PHI)
C     SCALE(LV)= SCALE FOR EACH IMAGE. SET TO 1.0 IF NEW.
C     TSHIFT(2,LV)= SHIFT FOR EACH IMAGE. SET TO 0.0 IF NEW.
C     FIRST= LOGICAL; .T.=START FROM SCRATCH, .F.=USE OLD OUTPUT TO
C           BEGIN REFINEMENT
C     NTTVL= TOTAL NUMBER OF POINTS FOUND IN THE ALL IMAGES COMBINED
C COMMON VARIABLES PASSED OUT:
C     PTACTIVE(LS,LV)= BOOLEAN ARRAY TELLING IF A POINT IS BEING USED IN
C           ALIGNMENT OR NOT FOR A VIEW
C     NUMPTS(LV)= HIGHEST NUMBERED MARKER TO BE FOUND IN AN IMAGE
C     NTVW= TOTAL NUMBER OF VIEWS. ALSO HIGHEST VALID INDEX OF XYPTS,
C           ANGLE, ETC.
C     NTPT= HIGHEST NUMBER OF MARKERS USED
C     CIR(2)= COORDS OF CENTER OF VIEW. NSAM/2+1, NROW/2+1
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRGETINFO(XYPTS, IREF, NUMBER, ANGLE, SCALE, TSHIFT,
     &              FIRST,NTTVL,PARAMQ,PTACTIVE,NUMPTS,NTVW,NTPT,CIR)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      PARAMETER (LV=300)
      PARAMETER (LS=256)
      PARAMETER (MAXREG=7)
      PARAMETER (MAXKEY=256)

      LOGICAL           PARAMQ(4), PTACTIVE(LS,LV)
      INTEGER           NUMPTS(LV)
      REAL              CIR(2)

      CHARACTER(LEN=MAXNAM) :: XYIN,XYIN1,XYIN2,CORFIL,CORIN,DOCTIT
      CHARACTER(LEN=MAXNAM) :: SERNAME

      DIMENSION         XYPTS(2,LS,LV),ANGLE(3,LV),SCALE(LV)
      DIMENSION         TSHIFT(2,LV)
      DIMENSION         DBUF(MAXREG,MAXKEY),PLIST(MAXREG)
      INTEGER           NUMBER(LV),MXY(2)
      CHARACTER * 1     NULL,ANSW
      LOGICAL           FAKE,FIRST

      DATA DTOR/0.017453292/
      DATA LUNDOC,LUN1/10,79/

      NULL = CHAR(0)

C      ASK FOR DATA FILE
       NMAX = LV

       CALL FILELIST(.TRUE.,LUNDOC,XYIN1,NA,NUMBER,NMAX,NVW,
     &     'PREFIX OF MARKER DOC FILES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C      NVW - TOTAL NUMBER OF IMAGES

       NLETP   = NA
       SERNAME = XYIN1

       CALL RDPRMI(IRFR,IDUM,NOT_USED,'REFERENCE FILE NUMBER')

C      GET CENTER OF IMAGE
       CALL RDPRIS(MXY(1),MXY(2),NOT_USED,
     &      'X,Y IMAGE DIMENSIONS',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CIR(1) = MXY(1) / 2 + 1
       CIR(2) = MXY(2) / 2 + 1

      CALL RDPRMC(ANSW,NC,.TRUE.,'REFINE SCALE (VS REFERENCE)? (Y/N)',
     &      NULL,IRTFLG)
      PARAMQ(3) = (ANSW(:1) .EQ. 'Y') 

      CALL RDPRMC(ANSW,NC,.TRUE.,
     &      'REFINE TILT ANGLE (VS REFERENCE)? (Y/N)', NULL,IRTFLG)
      PARAMQ(2) = (ANSW(:1) .EQ. 'Y') 

      CALL RDPRMC(ANSW,NC,.TRUE.,'REFINE IN-PLANE ROTATION? (Y/N)',
     &      NULL,IRTFLG)
      PARAMQ(1) = (ANSW(:1) .EQ. 'Y') 

      CALL RDPRMC(ANSW,NC,.TRUE.,
     &      'REFINE SHIFT (VS REFERENCE)? (Y/N)',NULL,IRTFLG)
      PARAMQ(4) = (ANSW(:1) .EQ. 'Y') 

      CALL RDPRMC(ANSW,NC,.TRUE.,'USE PREVIOUS CORRECTIONS? (Y/N)',
     &           NULL,IRTFLG)
      FIRST = (ANSW(:1) .EQ. 'Y') 

      IF (FIRST) THEN
         CALL FILERD(CORFIL,NLTC,NULL,'CORRECTIONS INPUT DOC',
     &        IRTFLG)
         IF (IRTFLG.NE.0) RETURN
      ENDIF

      NTPT  = 0
      IVIEW = 0
      NTVW  = 0
      NTTVL = 0

      IRT2  = 0    ! SAY OPENDOC INFO

      NLETFIL = 0
      DO IVIEWT=1,NVW
C       GO THROUGH AND PULL ALL X,Y COORDS FROM DOC FILES

        CALL FILGET(XYIN1,XYIN2,NLETFIL,NUMBER(IVIEWT),IRTFLG1)
        IF (IRTFLG1 .NE. 0) CYCLE   ! BAD FILE NAME GENERATION

        IF (NUMBER(IVIEWT) .EQ. IRFR) IREF = IVIEWT
        XYIN2(NLETP+1:NLETP+1) = CHAR(0)

        IOPEN = 0
        IKEY  = 1
        CALL FILNAMANDEXT(XYIN2,DATEXC,XYIN,NLET,.TRUE.,IER)

        IF (IVIEWT > 1) IRT2  = -9    ! DO NOT SAY ANY MORE OPENDOC INFO
        CALL UNSDAL(XYIN,IOPEN,LUN1,IKEY,PLIST,3,DBUF,
     &              MAXKEY,MAXREG,NKEY2,IRT2)

C       UNSDAL CLOSES THE FILE
C       IRT2 WILL BE 1 IF KEY 1 IS MISSING, SO OK IF IRT2=0 OR 1
        IF (IRT2 .NE. 0 .AND. IRT2 .NE. 1) CYCLE

        IVIEW = IVIEW + 1
        NTVW  = NTVW  + 1
        IF (NTVW .GT. LV) THEN
           CALL ERRT(102,'NUMBER OF VIEWS CAN NOT EXCEED',LV)
           RETURN
        ENDIF

C       GET X,Y COORDS OF EACH MARKER FROM ALL VIEWS. A VIEW
C       CAN LACK SOME MARKERS
        NUMPTS(IVIEWT) = 0
        NUSED          = 0

        DO I=1,NKEY2
C         IF COORD OUTSIDE RANGE, THE POINT IS FAKE AND NOT USED
          FAKE = DBUF(2,I).LT.0 .OR. DBUF(2,I) .GT. (2*CIR(1)) .OR.
     &           DBUF(3,I).LT.0 .OR. DBUF(3,I) .GT. (2*CIR(2))
C         OR WHEN THE KEY IS MISSING ...
     &    .OR.  DBUF(1,I).EQ.0.0

          IF (FAKE) THEN
	    WRITE(NOUT,*)  '  FAKE POINT !!  ',IVIEWT,I
            PTACTIVE(I,IVIEWT) = .FALSE.
            XYPTS(1,I,IVIEWT)  = 0.0
            XYPTS(2,I,IVIEWT)  = 0.0
          ELSE
            PTACTIVE(I,IVIEWT) = .TRUE.
            NTTVL              = NTTVL + 1
            XYPTS(1,I,IVIEWT)  = DBUF(2,I )- CIR(1)
            XYPTS(2,I,IVIEWT)  = DBUF(3,I) - CIR(2)
            IF (I .GT. NUMPTS(IVIEW)) NUMPTS(IVIEW) = I
          ENDIF
        ENDDO

        IF (NUMPTS(IVIEW) .GT. NTPT) NTPT = NUMPTS(IVIEW)

C       GET INITIAL ANGLES FROM DOC FILE
        IOPEN = 0
        IKEY  = -1
        IRT2  = -9
        CALL UNSDAL(XYIN,IOPEN,LUN1,IKEY,PLIST,3,DBUF,
     &                 MAXKEY,MAXREG,NKEY2,IRT2)
c       UNSDAL CLOSES THE FILE

        ANGLE(1,IVIEWT)  = PLIST(1) * DTOR
        ANGLE(2,IVIEWT)  = PLIST(2) * DTOR
        ANGLE(3,IVIEWT)  = PLIST(3) * DTOR
        SCALE(IVIEWT)    = 1.000000
        TSHIFT(1,IVIEWT) = 0.0
        TSHIFT(2,IVIEWT) = 0.0
      ENDDO

      IF (.NOT. FIRST) GOTO 999

      IOPEN = 0
      IKEY  = 1
      IC    = 1
      IRT2  = 0
      CALL FILNAMANDEXT(CORFIL,DATEXC,CORIN,NLET,.TRUE.,IER)

      CALL UNSDAL(CORIN,IOPEN,LUN1,IKEY,PLIST,3,DBUF,
     &              MAXKEY,MAXREG,NKEY2,IRT2)

C     IRT2 WILL BE 1 IF KEY 1 IS MISSING, SO OK IF IRT2=0 OR 1
      IF (IRT2 .NE. 0 .AND. IRT2 .NE. 1) THEN
         WRITE(NOUT)' NO CORRECTIONS FILE. CONTINUING AS FIRST RUN'
         FIRST = .FALSE.
         GOTO 999
      ENDIF

C     COMMENTS TEMP MODIFIED TO READ IN MARKER ANGLES FILE
      DO  J=1,NKEY2
        IVIEW           = DBUF(1,J)
        SCALE(IVIEW)    = DBUF(2,J)
        ANGLE(2,IVIEW)  = DBUF(3,J)*DTOR
        ANGLE(1,IVIEW)  = DBUF(4,J)*DTOR
        TSHIFT(1,IVIEW) = DBUF(5,J)
        TSHIFT(2,IVIEW) = DBUF(6,J)
      ENDDO

 999  RETURN
      END
