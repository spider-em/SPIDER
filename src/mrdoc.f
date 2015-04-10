
C ++********************************************************************
C                                                                      *
C MRDOC                                                                   *
C          PROMPTS & DOC FILE HEADERS IMPROVED   FEB 2014 ARDEAN LEITH
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
C  MRDOC(SCALE,SHIFT,ANGLE,NUMBER,P3D,NTVW,NTPT)                       *
C                                                                      *
C  PURPOSE:  RECORDS THE ROTATIONS, SCALE, AND SHIFTS OF THE IMAGES    *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRDOC(SCALE,SHIFT,ANGLE,NUMBER,P3D,NTVW,NTPT)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      PARAMETER (LV=300)
      PARAMETER (LS=256)

      CHARACTER(LEN=MAXNAM) :: SERNAME,CORDOC,ZDOC
      CHARACTER(LEN=1)      :: NULL = CHAR(0)

      LOGICAL               :: ADDEXT,ASKNAM,ISOLD,APPEND
      LOGICAL               :: MESSAGE,NEWFILE

      REAL                  :: SCALE(LV),SHIFT(2,LV),ANGLE(3,LV)
      REAL                  :: P3D(3,LS)
      REAL                  :: BUFLOC(7)
      INTEGER               :: NUMBER(LV)

      DATA PI/3.141592654/,NCOR/12/,NXYZ/13/
      
      ASKNAM   = .TRUE.
      ADDEXT   = .TRUE.
      ISOLD    = .FALSE.
      APPEND   = .FALSE. 
      MESSAGE  = .TRUE. 
      CALL OPENDOC(CORDOC,ADDEXT,NLET,NCOR,NCORT,ASKNAM,
     &            'CORRECTIONS OUTPUT DOC FILE',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN   
        
      CORDOC = 'FOR IMAGE SERIES: ' // CORDOC(1:NLET)
      CALL LUNDOCPUTCOM(NCORT,CORDOC,IRTFLG)

      CORDOC = 
     & 'VIEW    SCALE         THETA          PSI         ' //
     & ' XSHIFT        YSHIFT'
      CALL LUNDOCPUTCOM(NCORT,CORDOC,IRTFLG)

      DO  IVIEW=1,NTVW
         BUFLOC(1) = SCALE(IVIEW)
         BUFLOC(2) = ANGLE(2,IVIEW)*180/PI
         BUFLOC(3) = ANGLE(1,IVIEW)*180/PI
         BUFLOC(4) = SHIFT(1,IVIEW)
         BUFLOC(5) = SHIFT(2,IVIEW)
         !CALL SAVDN1(NCOR,CORDOC,BUFLOC,NLIST,NRUN,0)
         CALL LUNDOCWRTDAT(NCORT,NUMBER(IVIEW),BUFLOC,5,IRTFLG)
      ENDDO

      CLOSE(NCOR)

      

      ASKNAM   = .TRUE.
      ADDEXT   = .TRUE.
      ISOLD    = .FALSE.
      APPEND   = .FALSE. 
      MESSAGE  = .TRUE. 
      CALL OPENDOC(ZDOC,ADDEXT,NLET,NCOR,NCORT,ASKNAM,
     &            '3-D COORDS OUTPUT DOC FILE',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN   
        
      CORDOC = 'FOR IMAGE SERIES: ' // ZDOC(1:NLET)
      CALL LUNDOCPUTCOM(NCORT,CORDOC,IRTFLG)

      CALL LUNDOCPUTCOM(NCORT,
     &     'POINT    XCOOR        YCOOR        ZCOOR',IRTFLG)

      DO IPOINT=1,NTPT
           BUFLOC(1) = P3D(1,IPOINT)
           BUFLOC(2) = P3D(2,IPOINT)
           BUFLOC(3) = P3D(3,IPOINT)
           !CALL SAVDN1(NXYZ,ZDOC,BUFLOC,NLIST,NRUN,0)
           CALL LUNDOCWRTDAT(NCORT,IPOINT,BUFLOC,3,IRTFLG)
      ENDDO

      CLOSE(NCOR)

      END

