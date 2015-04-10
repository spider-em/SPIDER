C++*********************************************************************
C
C RADAV.F    ADDED 'RO SD'                      MAY 2012  ARDEAN LEITH
C            ADDED RITHALF                      NOV 2013  ARDEAN LEITH
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
C  RADAV(LUN1,LUN2)
C
C  PURPOSE: ROTATIONAL AVERAGING INTO A SINGLE LINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

      SUBROUTINE RADAV()

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=MAXNAM)   :: FILNAM,FILDOC
      CHARACTER(LEN=80)       :: COMMEN
      INTEGER                 :: MAXIM1,ITYPE,NX1,NY1,NZ1,NX2,NY2,NZ2
      INTEGER                 :: MAXIM2,IRMAX,IRMAX2,IRMAX3
      INTEGER                 :: NOUTANG,NLET,IRTFLG,IERR,KEY,IHALF
      LOGICAL                 :: NEWFILE
      REAL                    :: DLIST(3)
      REAL,ALLOCATABLE        :: BUF(:)
      CHARACTER(LEN=1)        :: HALF

      INTEGER, PARAMETER      :: LUN1   = 21
      INTEGER, PARAMETER      :: LUN2   = 22
      INTEGER, PARAMETER      :: LUNDOC = 80

C     OPEN INPUT FILE
      MAXIM1 = 0
      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,NX1,NY1,NZ1,
     &               MAXIM1,'INPUT',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IRMAX  = NX1/2 + MOD(NX1,2)
      IRMAX2 = NY1/2 + MOD(NY1,2)
      IRMAX3 = NZ1/2 + MOD(NZ1,2)

      IF (NZ1 == 1)THEN
         NX2 = MIN(IRMAX,IRMAX2)
      ELSE
         NX2 = MIN(IRMAX,IRMAX2,IRMAX3)
      ENDIF

      NY2 = 1
      NZ2 = 1
           
C     OPEN AN OUTPUT FILE
      MAXIM2 = 0 
      ITYPE  = 1
      CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NX2,NY2,NZ2,
     &  	     MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      IRMAX  =  NX2
  
      ALLOCATE(BUF (IRMAX), STAT=IERR)
      IF (IERR .NE. 0) THEN
      	 CALL ERRT(46,'RADAV; BUF',IRMAX)
         GOTO 9999
      ENDIF

      IF (NZ1 == 1) THEN
C        IMAGE

C        CHECK FOR HALF IMAGE AVERAGE (L or R)
         HALF = FCHAR(7:7)

         CALL CRCSE1(LUN1,LUN2, NX1,NY1, IRMAX,BUF,HALF)
      ELSE
C        VOLUME
         CALL CRCSE3(LUN1,LUN2, NX1,NY1,NZ1, IRMAX,BUF)
      ENDIF

 

      IF (FCHAR(4:4) .NE. 'S') GOTO 9999

C     OPEN OUTPUT DOC FILE (FOR APPENDING)
      NOUTANG = LUNDOC
      CALL OPENDOC(FILDOC,.TRUE.,NLET,LUNDOC,NOUTANG,.TRUE.,
     &             'OUTPUT  DOCUMENT',.FALSE.,.TRUE.,.TRUE.,
     &             NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      IF (NEWFILE) THEN

C        LABELS FOR COLUMNS IN OUTPUT DOC FILE
         COMMEN = '        CONTENTS:  ROTATIONAL AVERAGE'
         CALL LUNDOCPUTCOM(LUNDOC,COMMEN,IRTFLG)
         COMMEN = '        AVERAGE       RADIUS      SPATIAL-FREQ'
         CALL LUNDOCPUTCOM(LUNDOC,COMMEN,IRTFLG)
      ENDIF

      DO KEY = 1,IRMAX

C        RADIAL AVERAGE
         DLIST(1) = BUF(KEY)

C        RADIUS
         DLIST(2) = KEY

C        SPATIAL FREQ.
         DLIST(3) = FLOAT(KEY) / (2.0 * IRMAX) 

         CALL LUNDOCWRTDAT(LUNDOC,KEY,DLIST,3,IRTFLG)
      ENDDO

9999  IF (ALLOCATED(BUF))   DEALLOCATE(BUF)
      CLOSE(LUNDOC)
      CLOSE(LUN1)
      CLOSE(LUN2)

      END
