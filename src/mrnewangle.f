
C ++********************************************************************
C                                                                      
C MRNEWANGL   NLIST FIXED                        AUG 2000 ARDEAN LEITH             
C             CAN HANDLE MISSING KEYS            JUL 2001 ARDEAN LEITH 
C             INCORE LUNDOC                      JUL 2003 ARDEAN LEITH 
C             PROMPTS & DOC FILE HEADERS         FEB 2014 ARDEAN LEITH
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
C                                                                      *
C  MRNEWANGLE(PSI,THETA,PHI)                                           *
C                                                                      *
C  PURPOSE:  USING ANGLES TO CHANGE VOLUME, AND THE TILT ANGLES
C            OF THE SECOND IMAGE SERIES, CALCULATES THE EULERIAN
C            ANGLES OF THIS SECOND SERIES SO IT MAY BE USED IN
C            CONJUNCTION WITH THE ORIGINAL SERIES TO CREATE A
C            3D IMAGE
C                                                                      *
C  PARAMETERS:   PSI,THETA,PHI    ANGLES                               *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRNEWANGLE(PSI,THETA,PHI)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      PARAMETER (MAXKEY=9992)
 
      REAL                  :: A1(3),A2(3),A3(3)
      CHARACTER(LEN=1)      :: NULL = CHAR(0)
      CHARACTER(LEN=MAXNAM) :: ANGFIL,CORFIL
      LOGICAL               :: ADDEXT,ASKNAM,ISOLD,APPEND
      LOGICAL               :: MESSAGE,NEWFILE

      COMMON /DOC_BUF/ DBUF(4,MAXKEY)

      DATA  DTOR/0.017453292/

      DATA  NDOCT,NDOUT/55,56/

C     PREPARE ARRAY FOR ANGLES
      A1(1) = PSI/DTOR
      A1(2) = THETA/DTOR
      A1(3) = PHI/DTOR

C     OPEN INPUT DOC FILE,  OPENDOC CAN CHANGE NDOC
      CALL OPENDOC(ANGFIL,.TRUE.,NLET,NDOCT,NDOC,.TRUE.,
     &             'SECOND SERIES ANGULAR INPUT DOC',
     &             .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     OPEN OUTPUT DOC FILE BUT DO NOT WRITE ANYTHING
      ASKNAM   = .TRUE.
      ADDEXT   = .TRUE.
      ISOLD    = .FALSE.
      APPEND   = .FALSE. 
      MESSAGE  = .TRUE. 
      CALL OPENDOC(CORFIL,ADDEXT,NLET,NDOUT,NDOUTT,ASKNAM,
     &            'CORRECTED ANGULAR OUTPUT DOC',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      CALL LUNDOCPUTCOM(NDOUTT,
     &   '         PHI          THETA         PSI ', IRTFLG)


C     THE FOLLOWING DOES THE SAME AS 'VO RA'

C     RETRIEVE DBUF FROM DOC FILE
      CALL LUNDOCREDALL(NDOC,DBUF,4,MAXKEY,.TRUE.,NGOT,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      A2(1) = 0.0
      A2(3) = 0

      DO K = 1,NGOT
         IF (DBUF(1,K) > 0) THEN
C           KEY EXISTS IN DOC FILE
            A2(2) = DBUF(3,K)

C           ROTATE THE ANGLES
            CALL CALD(A1,A2,A3)

C           OUTPUT THE ANGLES
            !CALL SAVDN1(NDOUT,CORFIL,A3,4,1,0)
            CALL LUNDOCWRTDAT(NDOUTT,K,A3,3,IRTFLG)
          ENDIF
      ENDDO

9999  WRITE(NOUT,*) ' '
      CLOSE(NDOUT)
      CLOSE(NDOCT)

      END
