C ++********************************************************************
C                                                AUTHOR: PAWEL PENCZEK  
C  MRK3                                                                
C          USED PARAMETERS INSTEAD OF COMMON     DEC 2000 ARDEAN LEITH 
C          USED OPENDOC                          DEC 2000 ARDEAN LEITH 
C          UNUSED SHIFT DEFAULTED                JAN 2000 ARDEAN LEITH 
C          PROMPTS & DOC FILE HEADERS IMPROVED   FEB 2014 ARDEAN LEITH  
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
C MRK3(MAXDIM)
C
C PURPOSE:  ALIGNS A 3D MODEL TO A REFERENCE 3D MODEL, BOTH CREATED
C           FROM OUTPUT FILES CREATED BY MKRALIGN. FINDS ALL
C           PARAMETERS TO TRANSFORM ONE MODEL TO ANOTHER (SCALE,
C           PHI, THETA, PSI) AND WILL CALCULATE THE NEW
C           PARAMETERS FOR THE DATA SET ASSOCIATED WITH THE MODEL.
C
C *********************************************************************

      SUBROUTINE MRK3(MAXDIM)

      INCLUDE     'CMLIMIT.INC'
      INCLUDE     'CMBLOCK.INC'

      PARAMETER (LS=2000,NLIST=5)

      DIMENSION                P3D1(3,LS), P3DREF(3,LS), ERRORP(LS)
      DIMENSION                CGV(3), CGR(3), BUF(NLIST)
      DIMENSION                SHIFT(3),ROT(3,3)
      CHARACTER(LEN=MAXNAM) :: ERRFILE,DOCFILE
      CHARACTER(LEN=1)      :: NULL = CHAR(0)

      LOGICAL               :: ADDEXT,ASKNAM,ISOLD,APPEND
      LOGICAL               :: MESSAGE,NEWFILE

      DATA  NDOUT/56/
      DATA  NDINT/55/

C     GET REFERENCE 3D MODEL
C     POINTS IN FILES ARE ALREADY ABOUT THE CENTER OF GRAVITY

       CALL OPENDOC(DOCFILE,.TRUE.,NLET,NDINT,NDIN,.TRUE.,
     &             'FIRST SERIES MARKER INPUT DOC',
     &            .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL LUNDOCREDSLC(NDIN,.FALSE.,IDUM,P3DREF,3,LS,
     &         .TRUE.,.FALSE.,1,3, 1,LS, NTOLD,MAXY1,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      CLOSE(NDINT)

      CALL OPENDOC(DOCFILE,.TRUE.,NLET,NDINT,NDIN,.TRUE.,
     &             'SECOND SERIES MARKER INPUT DOC',
     &             .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL LUNDOCREDSLC(NDIN,.FALSE.,IDUM,P3D1,3,NTOLD,
     &         .TRUE.,.TRUE.,1,3, 1,NTOLD, NTPT,MAXY2,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      CLOSE(NDINT)

C     CENTER POINTS ABOUT CENTER OF GRAVITY, GET COORDS
C     OF CENTER OF GRAVITY FOR EACH IMAGE

      CALL MRCG3(P3D1,  CGV,NTPT,LS)

      CALL MRCG3(P3DREF,CGR,NTPT,LS)

      WRITE(NOUT,1001)  CGV
1001  FORMAT(/,'  CENTER OF GRAVITY OF SET TO BE ALIGNED:',/,
     &		3(3X,1PE10.3))

      WRITE(NOUT,1002)  CGR
1002  FORMAT(/,'  CENTER OF GRAVITY OF REFERENCE SET:',/,
     &		3(3X,1PE10.3))

      WRITE(NOUT,1003)
1003  FORMAT(/,
     &  '  WARNING!  BOTH CENTERS OF GRAVITY SHOULD BE NEAR ZERO.',/,
     &  '            THEY ARE NOT CORRECTED IN THIS PROGRAM.')

      CALL MRQUATER(P3DREF,P3D1,ROT,NTPT,LS)

      CALL MRROTATE(ROT,PHI,THETA,PSI)

      CALL MRSCALE3(P3DREF,P3D1,SCALE,NTPT,LS)

      CALL MRERROR(P3DREF,P3D1,ERRORP,NTPT,
     &             PHI,THETA,PSI,SHIFT,SCALE,ERG,NOUT)

C     SHIFT IS NOT SET
      SHIFT(1) = 0.0
      SHIFT(2) = 0.0
      SHIFT(3) = 0.0

      CALL MRDOCCOR(SCALE,SHIFT,PSI,THETA,PHI)

      CALL MRNEWANGLE(PSI,THETA,PHI)

C     OPEN ERROR DOC. FILE

      ASKNAM   = .TRUE.
      ADDEXT   = .TRUE.
      ISOLD    = .FALSE.
      APPEND   = .FALSE. 
      MESSAGE  = .TRUE. 
      CALL OPENDOC(ERRFILE,ADDEXT,NLET,NDOUT,NDOUTT,ASKNAM,
     &            'ALIGNMENT / ERRORS OUTPUT DOC',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL LUNDOCPUTCOM(NDOUTT,
     &   '          X             Y            Z           ERR/PT', 
     &   IRTFLG)


C     STORE ALIGNED COORDINATES, ERROR PER POINT, AND THE TOTAL ERROR
      DO I=1,NTPT
	  DO M=1,3
	     BUF(M) = P3D1(M,I)
	  ENDDO
	  BUF(4) = ERRORP(I)

          CALL LUNDOCWRTDAT(NDOUTT,I,BUF,4,IRTFLG)

      ENDDO

      CALL LUNDOCPUTCOM(NDOUTT,'          TOTAL_ERRORR', IRTFLG)

      BUF(1) = ERG
      !CALL SAVDN1(NDOUT,ERRFILE,BUF,2,1,0)
      CALL LUNDOCWRTDAT(NDOUTT,-1,BUF,1,IRTFLG)

9999  CLOSE(NDOUT)

      END
