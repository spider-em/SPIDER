
C ++********************************************************************
C                                                                      
C  MRDOCCOR    
C       PROMPTS AND DOC FILE HEADERS IMPROVED   FEB 2014 ARDEAN LEITH 
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
C MRDOCCOR(SCALE,SHIFT,PSI,THETA,PHI) 
C   
C PURPOSE: MARKER BASED ALIGNMENT - DOUBLE TILTED IMAGES
C
C PARAMETERS: 
C 
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRDOCCOR(SCALE,SHIFT,PSI,THETA,PHI)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=MAXNAM) :: CORDOC
      CHARACTER(LEN=1)      :: NULL = CHAR(0)
      INTEGER               :: IRTFLG,NCOR,NCORT,NLET
      REAL                  :: BUFLOC(7),SHIFT(3)

      LOGICAL               :: ADDEXT,ASKNAM,ISOLD,APPEND
      LOGICAL               :: MESSAGE,NEWFILE

      DATA    PI/3.141592654/, NCOR/12/

      DPSI   = PSI*180/PI
      DTHETA = THETA*180/PI
      DPHI   = PHI*180/PI

      WRITE(NOUT,123)
 123  FORMAT('  TO CORRECT, ROTATE USING  <VO RA> ',/,
     &       '     PSI,      THETA,      PHI')

      WRITE(NOUT,121) DPSI,DTHETA,DPHI
 121  FORMAT(2X,F7.2,4X,F7.2,4X,F7.2)


      WRITE(NOUT,124)
 124  FORMAT(/,'  SCALE,   XSHIFT,  YSHIFT,  ZSHIFT')

      WRITE(NOUT,122) SCALE,SHIFT(1),SHIFT(2),SHIFT(3)
 122  FORMAT(2X,F5.3,3X,F6.1,3X,F6.1,3X,F6.1,/)

      ASKNAM   = .TRUE.
      ADDEXT   = .TRUE.
      ISOLD    = .FALSE.
      APPEND   = .FALSE. 
      MESSAGE  = .TRUE. 
      CALL OPENDOC(CORDOC,ADDEXT,NLET,NCOR,NCORT,ASKNAM,
     &            'CORRECTIONS OUTPUT DOC FILE',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      CALL LUNDOCPUTCOM(NCORT,
     &   ' TO CORRECT, ROTATE USING  <VO RA> by : PSI,   THETA,   PHI', 
     &   IRTFLG)

      CALL LUNDOCPUTCOM(NCORT, 
     &   '         PSI,         THETA,         PHI',IRTFLG)

      BUFLOC(1) = PSI   * 180/PI
      BUFLOC(2) = THETA * 180/PI
      BUFLOC(3) = PHI   * 180/PI

      !CALL SAVDN1(NCOR,CORDOC,BUFLOC,NLIST,NRUN,0)
      CALL LUNDOCWRTDAT(NCORT,1,BUFLOC,3,IRTFLG)

      CALL LUNDOCPUTCOM(NCORT,
     &     '        SCALE,       XSHIFT,       YSHIFT,        ZSHIFT',
     &     IRTFLG) 

      BUFLOC(1) = SCALE
      BUFLOC(2) = 0.0
      BUFLOC(3) = 0.0
      BUFLOC(4) = 0.0
      !CALL SAVDN1(NCOR,CORDOC,BUFLOC,NLIST,NRUN,0)
      CALL LUNDOCWRTDAT(NCORT,2,BUFLOC,4,IRTFLG)

9999  CLOSE(NCOR)

      END
