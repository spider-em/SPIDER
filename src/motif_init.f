C++*********************************************************************
C
C MOTIF_INIT        NEW                               JUL 08 ARDEAN LEITH
C
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
C                  
C  MOTIF_INIT(ANS,LUNVOL,LUNVOLFFT,LUNVOLSQFFT, LUNMSK,LUNMSKFFT, 
C             NSAM,NROW,NSLICE, NSAMM,NROWM,NSLICEM)

C  PURPOSE: CREATE FFT SEARCH VOLUME, FFT OF SQUARED SEARCH VOLUME
C           AND FFT OF ROTATED PADDED MASK, AS REQUESTED. 
C           THESE FILES ARE NEEDED FOR: 'LO LSD' AND 'LO'

C PARAMETERS:
C       ANS                          OP TYPE FLAG              (INPUT)
C       LUNVOL...                    I/O UNITS                 (INPUT)
C       NSAM,NROW,NSLICE             VOLUME SIZE               (INPUT)
C       NSAMM,NROWM,NSLICEM          MASK SIZE                (INPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE MOTIF_INIT(ANS,LUNVOL,LUNVOLFFT,LUNVOLSQFFT,
     &                        LUNMSK,LUNMSKFFT, 
     &                        NSAM,NROW,NSLICE, NSAMM,NROWM,NSLICEM)

#ifdef SP_LIBFFTW3
        USE TYPE_KINDS
        INTEGER(KIND=I_8)         :: IPLAN = 0  !STRUCTURE POINTER 
#endif

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER, INTENT(IN)       :: LUNVOL,LUNVOLFFT,LUNVOLSQFFT
	INTEGER, INTENT(IN)       :: LUNMSK,LUNMSKFFT
	INTEGER, INTENT(IN)       :: NSAM,NROW,NSLICE 
	INTEGER, INTENT(IN)       :: NSAMM,NROWM,NSLICEM 
	CHARACTER, INTENT(IN)     :: ANS 
 
        LOGICAL                   :: SPIDER_SIGN  = .TRUE.
        LOGICAL                   :: SPIDER_SCALE = .TRUE.

C       ALLOCATABLE ARRAYS
	REAL, ALLOCATABLE         :: BUFVOL(:,:,:)


C       EXTRA COLUMN SPACE FOR FOURIER TRANSFORM
        LSE = NSAM + 2 - MOD(NSAM,2)

        MWANT = LSE*NROW*NSLICE 
 	ALLOCATE (BUFVOL(LSE,NROW,NSLICE), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'MOTIF_INIT; BUFVOL...',MWANT)
           GOTO 9999
        ENDIF

        IF (ANS .EQ. 'M' .OR. ANS .EQ. 'B') THEN

C          PAD ROTATED MASK INTO SAME SIZE AS FFT OF LARGE VOLUME ------
           CALL REDNPADVOL(LUNMSK, 0.0, NSAMM,NROWM,NSLICEM,
     &                                  LSE,NROW,NSLICE,
     &                                  BUFVOL, IRTFLG)
	   IF (IRTFLG .NE. 0)  GOTO 9999
           
C          FIND NO. OF NON-ZERO PIXELS INSIDE MASK
           PIX_IN_MASK = SUM(BUFVOL(1:NSAMM, 1:NROWM, 1:NSLICEM))
           X51 = 1.0 / PIX_IN_MASK        ! FOR SPEED

           IF (VERBOSE) THEN
              WRITE(NOUT,*) ' Non-zero mask pixels:', PIX_IN_MASK
              !write(6,*)    ' Non-zero mask pixels:',PIX_IN_MASK
           ENDIF

           CALL REG_GET_USED(NSEL_USED)
           IF (NSEL_USED .GT. 0) THEN
C             OUTPUT X51 TO SPIDER'S REGISTERS
              CALL REG_SET_NSEL(1,1, X51, 0.0, 0.0, 0.0, 0.0,IRTFLG)
           ENDIF

C          FFT ON THIS PADDED ROTATED MASK 
           INV = +1
           CALL FMRS(BUFVOL, NSAM,NROW,NSLICE, IPLAN,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,'MOTIF_INIT; FFT ERROR ON PADDED MASK',NE)
              GOTO 9999
           ENDIF

C          OUTPUT FFT OF THIS PADDED ROTATED MASK VOLUME      
           CALL WRTVOL(LUNMSKFFT,LSE,NROW,1,NSLICE,BUFVOL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

        IF (ANS .NE. 'M') THEN
C           WISH TO DO FFT'S ON SEARCH VOLUME -------------------------

C          PAD SEARCH VOLUME INTO FFT VOLUME ---------------------- 
           CALL READV(LUNVOL,BUFVOL, LSE,NROW, NSAM,NROW,NSLICE)

           !fmin = minval(bufvol(1:nsam, 1:nrow, 1:nslice))
           !fmax = maxval(bufvol(1:nsam, 1:nrow, 1:nslice))
           !write(6,*) ' Search volume:              ',fmin,'...',fmax

C          FFT ON SEARCH VOLUME
           INV = +1
           CALL FMRS(BUFVOL, NSAM,NROW,NSLICE, IPLAN,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,'MOTIF_INIT; FFT ERROR ON VOLUME',NE)
              GOTO 9999
           ENDIF

C          OUTPUT FFT OF SEARCH VOLUME      
           CALL WRTVOL(LUNVOLFFT,LSE,NROW,1,NSLICE,BUFVOL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          PAD LARGE VOLUME INTO FFT VOLUME ----------------------

C          IT SAVES SPACE TO RE-READ?
           CALL READV(LUNVOL,BUFVOL, LSE,NROW, NSAM,NROW,NSLICE)

           !fmin = minval(bufvol(1:nsam, 1:nrow, 1:nslice))
           !fmax = maxval(bufvol(1:nsam, 1:nrow, 1:nslice))
           !write(6,*) ' Search volume:               ',fmin,'...',fmax

c$omp      parallel do private(i,j,k)
           DO K = 1,NSLICE
              DO J = 1,NROW
                 DO I = 1,NSAM
                    BUFVOL(I,J,K) = BUFVOL(I,J,K) ** 2  
                 ENDDO
              ENDDO
           ENDDO
c$omp      end parallel do

           !fmin = minval(bufvol(1:nsam, 1:nrow, 1:nslice))
           !fmax = maxval(bufvol(1:nsam, 1:nrow, 1:nslice))
           !write(6,*) ' Squared search vol:         ',fmin,'...',fmax

C          FFT ON SQUARED SEARCH VOLUME -------------------------
           INV = +1
           CALL FMRS(BUFVOL, NSAM,NROW,NSLICE, IPLAN,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,'MOTIF_INIT FFT ERROR ON SQUARED VOLUME',NE)
              GOTO 9999
           ENDIF

C          OUTPUT FFT OF SQUARED SEARCH VOLUME      
           CALL WRTVOL(LUNVOLSQFFT,LSE,NROW,1,NSLICE,BUFVOL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

9999    IF (ALLOCATED(BUFVOL)) DEALLOCATE(BUFVOL)

        END

