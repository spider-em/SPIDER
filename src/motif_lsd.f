cpgi$g opt=O3

C++*********************************************************************
C
C MOTIF_LSD           NEW                          JUN 08 ARDEAN LEITH
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
C  MOTIF_LSD(LUNMSK,LUNMSKFFT,LUNVOL,LUNVOLSQ,LUNLSD, 
C                   NSAM,NROW,NSLICE, IRTFLG)
C
C  PURPOSE:  PREPARE LSD VOLUME FOR MOTIF SEARCH.  USES FOUR VOLUMES
C            OF SPACE. COULD CUT THE MEMORY BY WRITING OUT TO DISK
C            AND OPERATING LINE BY LINE?
C  
C PARAMETERS:
C       LUNMSK....                   I/O UNITS                 (INPUT)
C       NSAM,NROW,NSLICE,            VOLUME SIZE               (INPUT)
C       NSAMM,NROWM,NSLICEM,         MASK SIZE                 (INPUT)
C       IRTFLG                       ERROR FLAG                (OUTPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE MOTIF_LSD(LUNMSK,LUNMSKFFT, LUNVOL,LUNVOLSQ,LUNLSD, 
     &                   NSAM,NROW,NSLICE, NSAMM,NROWM,NSLICEM, IRTFLG)

#ifdef SP_LIBFFTW3
        USE TYPE_KINDS
        INTEGER(KIND=I_8)     :: IPLAN = 0  !STRUCTURE POINTER 
#endif

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER, INTENT(IN)  :: LUNMSK,LUNMSKFFT 
	INTEGER, INTENT(IN)  :: LUNVOL,LUNVOLSQ,LUNLSD 
	INTEGER, INTENT(IN)  :: NSAM,NROW,NSLICE
 
c       LOGICAL              :: SPIDER_SIGN  = .FALSE.
c       LOGICAL              :: SPIDER_SCALE = .FALSE.
        LOGICAL              :: SPIDER_SIGN  = .TRUE.
        LOGICAL              :: SPIDER_SCALE = .TRUE.

        DOUBLE PRECISION     :: X51,TEMP

C       ALLOCATABLE ARRAYS
	REAL, ALLOCATABLE    :: BUFMP(:,:,:),   BUFVOL(:,:,:)
	REAL, ALLOCATABLE    :: BUFVOLSQ(:,:,:),BUFMUL(:,:,:)

C       EXTRA COLUMN SPACE FOR FOURIER TRANSFORM
        LSE = NSAM + 2 - MOD(NSAM,2)

C       MAKE VOLUMES FOR PADDED IMAGES WITH FOURIER ROW LENGTH
	ALLOCATE(BUFVOL  (LSE,NROW,NSLICE), 
     &           BUFMP   (LSE,NROW,NSLICE),
     &           BUFMUL  (LSE,NROW,NSLICE),
     &           BUFVOLSQ(LSE,NROW,NSLICE),
     &           STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = 4 * LSE*NROW*NSLICE    
           CALL ERRT(46,'MOTIF_LSD; BUFVOL...',MWANT)
           GOTO 9999
        ENDIF

C       READ MOTIF MASK VOLUME (BUFMP MUCH LARGER, HAS PADDING SPACE) 
        CALL READV(LUNMSK,BUFMP,LSE,NROW, NSAMM,NROWM,NSLICEM)

C       SUM OF MASK PIXELS (ONLY USE FOR: MOTIF MASK VOL. HERE)
        PIX_IN_MASK = SUM(BUFMP(1:NSAMM, 1:NROWM, 1:NSLICEM))
        X51 = 1.0 / PIX_IN_MASK        ! FOR SPEED

C       READ FFT OF SEARCH VOLUME 
        CALL REDVOL(LUNVOL, LSE,NROW, 1,NSLICE, BUFVOL,IRTFLG)
	IF (IRTFLG .NE. 0)  GOTO 9999

C       INPUT FFT OF PADDED ROTATED MASK VOLUME 
        CALL REDVOL(LUNMSKFFT, LSE,NROW, 1,NSLICE, BUFMP,IRTFLG)
	IF (IRTFLG .NE. 0)  GOTO 9999
        
#ifdef NEVER
        write(6,*)  ' Inverse FFT on padded rotated mask volume'
        inv =  -1
        call fmrs(bufmp, nsam,nrow,nslice, iplan,
     &            spider_sign,spider_scale, inv,irtflg)
        if (irtflg .ne. 0) stop
        fmin = minval(bufmp(1:nsam, 1:nrow, 1:nslice))
        fmax = maxval(bufmp(1:nsam, 1:nrow, 1:nslice))
        write(6,*) ' Padded rot, mask vol:',fmin,'...',fmax
        call writev(lunlsd,bufmp, lse,nrow, nsam,nrow,nslice)
#endif
                
        !write(6,*)  ' Mult big vol fft * padded rot mask fft'
C       MULTIPLY FFT OF VOLUME WITH COMPLEX CONJUGATE OF PADDED MASK FFT
        NCOMP = LSE * NROW * NSLICE / 2
        CALL CMULCON(BUFVOL,BUFMP,BUFMUL,NCOMP)

        !write(6,*)  ' Inverse FFT on multiplied volume'
C       INVERSE FFT ON MUTIPLIED VOLUME
        INV =  -1
        CALL FMRS(BUFMUL, NSAM,NROW,NSLICE, IPLAN,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(101,'MOTIF_LSD; FFT ERROR ON BUFMUL',NE)
           RETURN
        ENDIF
#ifdef NEVER
        fmin = minval(bufmul(1:nsam, 1:nrow, 1:nslice))
        fmax = maxval(bufmul(1:nsam, 1:nrow, 1:nslice))
        write(6,*)' Fftd search*mask :  ',fmin,'...',fmax
        call writev(lunlsd,bufmul, lse,nrow, nsam,nrow,nslice)
#endif

        !write(6,*)  ' Normalize the multiplied vol:'
C       NORMALIZE  P2 = (P1 / (No. OF NON-ZERO PIXELS INSIDE MASK))**2 
        NSIZE = NSAM * NROW * NSLICE
c$omp   parallel do private(i,j,k)
        DO K =       1,NSLICE
           DO J =    1,NROW
              DO I = 1,NSAM
                 BUFMUL(I,J,K) = (BUFMUL(I,J,K) * X51) **2
              ENDDO
           ENDDO
        ENDDO
c$omp   end parallel do

#ifdef NEVER
        fmin = minval(bufmul(1:nsam, 1:nrow, 1:nslice))
        fmax = maxval(bufmul(1:nsam, 1:nrow, 1:nslice))
        write(6,*)' Normalized bufmul :  ',fmin,'...',fmax
        call writev(lunlsd,bufmul, lse,nrow, nsam,nrow,nslice)
#endif

C       LOAD FFT OF SEARCH VOLUME SQUARED
        CALL REDVOL(LUNVOLSQ, LSE,NROW,1,NSLICE,BUFVOLSQ,IRTFLG)
	IF (IRTFLG .NE. 0)  GOTO 9999

        !write(6,*)  ' Multiplying search vol sq'
C       MULTIPLY FFT OF SQUARE OF SEARCH VOLUME WITH COMPLEX CONJUGATE 
C       OF FFT OF ROT PADDED MASK VOLUME
        CALL CMULCON(BUFVOLSQ,BUFMP,BUFVOL,NCOMP)

        !write(6,*)  ' Inverse FFT of mult search vol sq'
C       INVERSE FFT ON SQUARE OF SEARCH VOLUME TIMES COMPLEX CONJUGATE
        INV = -1
        CALL FMRS(BUFVOL, NSAM,NROW,NSLICE, IPLAN,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(101,'MOTIF_LSD; REVERSE FFT ERROR ON BUFVOL',NE)
           GOTO 9999
        ENDIF

#ifdef NEVER
        fmin = minval(bufvol(1:nsam, 1:nrow, 1:nslice))
        fmax = maxval(bufvol(1:nsam, 1:nrow, 1:nslice))
        write(6,*)' Bufvolsq:  ',fmin,'...',fmax
        call writev(lunlsd,bufvol, lse,nrow, nsam,nrow,nslice)
#endif

C       From Bimal's comments in SPIDER script:
C       Get rid of SQRT of -ve # and division by zero (while dividing
C       the CCF by local standard deviation). If you find CCC > 1.0
C       and the locations of the motif outside the outline of the
C       searched large volume then look at the file _25 and you will
C       find pixel values quite small ~ < 1e-10 but not equal to
C       zero. In this case, change the mask threshold in the
C       following operation from zero (0) to something like 1e-10.
C       this may solve the problem. However, you may need to change
C       this parameter to find the correct one to be used. This
C       problem of getting -ve variance or getting incorrect variance
C       when pixels under the mask have same/very_close values is
C       due to the way variance is calculated using Fourier transform. 

        !write(6,*) ' x51:',x51
      
C       NORMALIZE, THRESHOLD, TAKE SQRT, AND INVERT 
        DO K =       1,NSLICE
           DO J =    1,NROW
              DO I = 1,LSE
C               NORMALIZE
                BUFVOL(I,J,K) = BUFVOL(I,J,K) * DBLE(X51) - 
     &                          DBLE(BUFMUL(I,J,K)) 

C               PREVENT DIVISION BY ZERO
                IF (BUFVOL(I,J,K) < 1E-10) BUFVOL(I,J,K) = 9E+20 

C               INVERT TO SAVE TIME IF LSD IS REUSED IN MOTIF.F
                BUFVOL(I,J,K) = DBLE(X51) / 
     &                          DBLE(SQRT(DBLE(BUFVOL(I,J,K))))
              ENDDO
           ENDDO
        ENDDO

        IF (VERBOSE) THEN
           FMIN = MINVAL(BUFVOL(1:NSAM, 1:NROW, 1:NSLICE))
           FMAX = MAXVAL(BUFVOL(1:NSAM, 1:NROW, 1:NSLICE))
           WRITE(6,*)    ' LSD vol range: ',FMIN,'...',FMAX
           WRITE(nout,*) ' LSD vol range: ',FMIN,'...',FMAX
        ENDIF
     
C       OUTPUT LOCAL STANDARD DEVIATION VOLUME (FROM PADDED BUFFER)     
        CALL WRITEV(LUNLSD,BUFVOL, LSE,NROW, NSAM,NROW,NSLICE)
 
9999    IF (ALLOCATED(BUFVOL))   DEALLOCATE(BUFVOL)
        IF (ALLOCATED(BUFMP))    DEALLOCATE(BUFMP)
        IF (ALLOCATED(BUFMUL))   DEALLOCATE(BUFMUL)
        IF (ALLOCATED(BUFVOLSQ)) DEALLOCATE(BUFVOLSQ)

        END



C++*********************************************************************
C
C  CMULCON.F                                WRITTEN JUL 08 ARDEAN LEITH
C
C **********************************************************************
C
C  PURPOSE: COMPLEX MULTIPLICATION OF 2 ARRAYS OF COMPLEX NUMBERS
C           USUALLY USED AS KLUDGE TO AVOID EQUIVALENCING 
C
C PARAMETERS:
C       BUFA         FIRST ARRAY                (INPUT)
C       BUFB         SECOND ARRAY               (INPUT)
C       BUFOUT       RETURN ARRAY               (OUTPUT)
C       NCOMP        ARRAY LENGTH               (INPUT)
C
C DANGER:   BUFOUT CAN NOT BE SAME AS: BUFA OR BUFB 
C
C **********************************************************************

       SUBROUTINE CMULCON(BUFA,BUFB,BUFOUT,NCOMP)

       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: NCOMP 

       COMPLEX, INTENT(IN)  :: BUFA(*), BUFB(*)
       COMPLEX, INTENT(OUT) :: BUFOUT(*)

       INTEGER              :: I 

C      COMPLEX 1 FOURIER MULTIPLICATION WITH 2 CONJUGATE
c$omp  parallel do private(i)
       DO I = 1,NCOMP
          BUFOUT(I) = BUFA(I) * CONJG(BUFB(I)) 
       ENDDO
c$omp  end parallel do

       END

