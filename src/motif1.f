
C++*********************************************************************
C
C    MOTIF1.F                                      JUL  08 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C  PURPOSE:   HANDLES 'LO' OPERATIONS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE MOTIF1()

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 

      CHARACTER(LEN=MAXNAM) :: FILNAM

      INTEGER               :: LUN1   = 21
      INTEGER               :: LUN2   = 22
      INTEGER               :: LUN3   = 23
      INTEGER               :: LUN4   = 24
      INTEGER               :: LUN5   = 25
      INTEGER               :: LUNDOC = 70

      LOGICAL               :: NEWFILE
      CHARACTER (LEN=1)     :: ANS
      CHARACTER (LEN=80)    :: COMMEN
      INTEGER               :: LUNDOCT
      INTEGER               :: NPEAKS = 10
      REAL                  :: PHI   = 0.0
      REAL                  :: THETA = 0.0
      REAL                  :: PSI   = 0.0
      LOGICAL, SAVE         :: MESSAGE = .TRUE. 

      SELECT CASE(FCHAR(1:4))

      CASE('LO I')

C          --------- INITIALIZE FOR LOCATE A MOTIF -------------- 'LO I'

           CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &        'CREATE SEARCH FFTS (S), MASK FFT (M), OR BOTH (B)',
     &        CDUM,IRTFLG)

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',
     &                  ITYPE,NSAM,NROW,NSLICE,MAXIM,
     &                  'SEARCH VOLUME',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           NNROW   = MOD(NROW,2)
           NNSLICE = MOD(NSLICE,2)
           IF (NNROW .NE. 0 .OR. NNSLICE .NE. 0) THEN
              WRITE(NOUT,*) ' FOR EFFICIENCY PLEASE PAD IMAGE TO EVEN ',
     &                      'NUMBER OF ROWS AND SLICES AND START AGAIN'
              CALL ERRT(101,'INVALID ROW OR SLICE DIMENSION',NE)
              GOTO 9999
           ELSEIF (ITYPE .NE. 3) THEN
              CALL ERRT(101,'ONLY WORKS ON VOLUMES',NE)
              GOTO 9999
           ENDIF

C          EXTRA COLUMN SPACE FOR FOURIER TRANSFORM
           LSE   = NSAM + 2 - MOD(NSAM,2)

           IF (ANS .EQ. 'S' .OR. ANS .EQ. 'B') THEN
C             OPEN LARGE FFT OUTPUT FILE
              MAXIM = 0
              ITYPE = -22    ! FOURIER
              CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'U',
     &                     ITYPE,LSE,NROW,NSLICE,MAXIM,         
     &                     'SEARCH VOLUME FFT OUTPUT',
     &                     .FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

C             OPEN LARGE SQUARE FFT OUTPUT FILE
              MAXIM = 0
              ITYPE = -22    ! FOURIER
              CALL OPFILEC(0,.TRUE.,FILNAM,LUN3,'U',
     &                     ITYPE,LSE,NROW,NSLICE,MAXIM,    
     &                    'SEARCH VOLUME SQ. FFT OUTPUT',
     &                    .FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
           ENDIF

           IF (ANS .EQ. 'M' .OR. ANS .EQ. 'B') THEN
C             OPEN  ROTATED MASK INPUT FILE
              MAXIM = 0
              CALL OPFILEC(0,.TRUE.,FILNAM,LUN4,'O',
     &                     ITYPE,NSAMM,NROWM,NSLICEM,MAXIM, 
     &                     'ROTATED MOTIF MASK INPUT',
     &                     .FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

C             OPEN PADDED ROTATED MASK FFT OUTPUT FILE
              MAXIM = 0
              ITYPE = -22    ! FOURIER
              CALL OPFILEC(0,.TRUE.,FILNAM,LUN5,'U',
     &                     ITYPE,LSE,NROW,NSLICE,MAXIM,
     &                     'PADDED MOTIF MASK FFT OUTPUT',
     &                     .TRUE.,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ENDIF

           CALL MOTIF_INIT(ANS,LUN1,LUN2,LUN3,LUN4,LUN5,
     &                     NSAM,NROW,NSLICE, NSAMM,NROWM,NSLICEM )

      CASE('LO L')

C          ------- CREATE LSD FOR LOCATING A MOTIF ------------ 'LO LSD'

C          OPEN ROTATED MASK INPUT FILE
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',
     &                  ITYPEM,NSAMM,NROWM,NSLICEM,MAXIM,
     &                  'ROTATED MOTIF MASK INPUT',
     &                  .FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          OPEN PADDED ROTATED MASK FFT INPUT FILE
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',
     &                  ITYPE,LSE,NROW,NSLICE,MAXIM,
     &                  'PADDED MOTIF MASK FFT INPUT',
     &                  .TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          OPEN SEARCH VOLUME FFT INPUT FILE
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN3,'O',
     &                  ITYPET,LSET,NROWT,NSLICET,MAXIM,
     &                  'SEARCH VOLUME FFT INPUT',
     &                  .TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          ENSURE SIZES ARE SAME SIZE
           CALL SIZCHK(UNUSED,LSE, NROW, NSLICE, ITYPE ,
     &                        LSET,NROWT,NSLICET,ITYPET,IRTFLG)
	   IF (IRTFLG .NE. 0)  GOTO 9999

C          OPEN SEARCH VOLUME SQUARE FFT INPUT FILE
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN4,'O',
     &                 ITYPE,LSET,NROWT,NSLICET,MAXIM,
     &                 'SEARCH VOLUME SQ. FFT INPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          ENSURE SIZES ARE SAME SIZE
           CALL SIZCHK(UNUSED,LSE, NROW, NSLICE, ITYPE, 
     &                        LSET,NROWT,NSLICET,ITYPET,IRTFLG)
	   IF (IRTFLG .NE. 0)  GOTO 9999

C          OPEN NEW LSD OUTPUT FILE
           NSAM    = LSE - 2
           ITYPE   = 3
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN5,'U',
     &                  ITYPE,NSAM,NROW,NSLICE,MAXIM,
     &                  'LSD VOLUME OUTPUT',.FALSE.,IRTFLG)
	   IF (IRTFLG .NE. 0)  GOTO 9999

           CALL MOTIF_LSD(LUN1,LUN2,LUN3,LUN4,LUN5, 
     &                    NSAM,NROW,NSLICE, NSAMM,NROWM,NSLICEM,
     &                    IRTFLG)

      CASE('LO  ')
C          --------- LOCATE A MOTIF ------------------------------ 'LO'

           CALL RDPRM3S(PHI,THETA,PSI,NOT_USED,
     &                  'PHI, THETA, & PSI',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',
     &                  ITYPEM,NSAMM,NROWM,NSLICEM,
     &                  MAXIM,'MOTIF INPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (NSLICEM < 2) THEN
              CALL ERRT(101,'ONLY WORKS ON VOLUMES',NE)
              GOTO 9999
           ENDIF

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',
     &                  ITYPEK,NSAMK,NROWK,NSLICEK,MAXIM,
     &                  'MASK FOR MOTIF INPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          ENSURE MASK AND MOTIF ARE SAME SIZE
           CALL SIZCHK(UNUSED,NSAMM,NROWM,NSLICEM,ITYPEM,
     &                        NSAMK,NROWK,NSLICEK,ITYPEK,IRTFLG)
	   IF (IRTFLG .NE. 0)  GOTO 9999

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN3,'O',
     &                  ITYPE,LSE,NROW,NSLICE,MAXIM,
     &                  'SEARCH VOLUME FFT INPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (NSLICE .LE. 1) THEN
              CALL ERRT(101,'OPERATION ONLY WORKS ON VOLUMES',NE)
              GOTO 9999
           ENDIF

           MAXIM  = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN4,'O',
     &                  ITYPEL,NSAML,NROWL,NSLICEL,MAXIM,
     &                  'LOCAL STANDARD DEV. INPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          ENSURE SEARCH VOLUME AND LSD VOLUME ARE SAME SIZE
           NSAM   = LSE - 2
           CALL SIZCHK(UNUSED,NSAML,NROWL,NSLICEL,0, 
     &                        NSAM, NROW, NSLICE, 0,IRTFLG)
	   IF (IRTFLG .NE. 0)  GOTO 9999

           CALL RDPRI1S(NPEAKS,NOT_USED,'NUMBER OF PEAKS WANTED',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          OPEN PEAK OUTPUT FILE
           MAXIM = 0
           ITYPE = 3
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN5,'U',
     &                  ITYPE,NSAM,NROW,NSLICE,MAXIM,
     &                  'PEAK OUTPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          OPEN OUTPUT DOC FILE (FOR APPENDING)
           LUNDOCT = LUNDOC
           CALL OPENDOC(FILNAM,.TRUE.,NLET,LUNDOC,LUNDOCT,.TRUE.,
     &           'PEAK OUTPUT DOCUMENT',.FALSE.,.TRUE.,MESSAGE,
     &            NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           MESSAGE = .FALSE.     ! NO COMMENT ON NEXT DOC FILE USAGE
           
C          COMMEN = '      X,            Y,             Z,          CC'
c          CALL LUNDOCPUTCOM(LUNDOCT,COMMEN,IRTFLG)

           CALL MOTIF(LUN1,LUN2,LUN3,LUN4,LUN5,LUNDOCT,
     &              NSAMM,NROWM,NSLICEM,  NSAM,NROW,NSLICE, NPEAKS,
     &              PHI,THETA,PSI)

C       ------------------- DEFAULT ------------------------------ '??'
      CASE DEFAULT
        CALL  ERRT(101,'UNKNOWN OPERATION',NE)
      END SELECT

9999  CLOSE(LUN5)
      CLOSE(LUN4)
      CLOSE(LUN3)
      CLOSE(LUN2)
      CLOSE(LUN1)
      CLOSE(LUNDOCT)

      END
