C++*********************************************************************
C
C MOTIF            NEW FROM PROCEDURE              JUN 08 ARDEAN LEITH
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
C MOTIF(LUNMOT,LUNMSK,LUNVOLFFT,LUNLSD,LUNPKS,LUNDOC, 
C       NSAMM,NROWM,NSLICEM, NSAM,NROW,NSLICE, NPEAKS, PHI,THETA,PSI)
C
C PURPOSE:  MAIN REPETITIVE UNIT IN MOTIF SEARCH
C           MASKS THE MOTIF AND PADS TO SEARCH VOLUME SIZE.  FFT'S
C           THE PADDED MOTIF AND COMPLEX MULTIPLIES BY THE FFT OF THE
C           SEARCH VOLUME.  REVERSE FFT'S THE PRODUCT AND MULTIPLIES BY
C           THE LOCAL STANDARD DEVIATION VOLUME.  SEARCHES FOR PEAKS
C           IN THE PRODUCT VOLUME AND LISTS PEAKS IN DOC FILE
C  
C PARAMETERS:
C       LUNMSK...                    I/O UNITS                 (INPUT)
C       NSAMM,NROWM,NSLICEM          MASK SIZE                 (INPUT)
C       NSAM,NROW,NSLICE,            VOLUME SIZE               (INPUT)
C       NPEAKS                       NUMBER OF PEAKS WANTED    (INPUT)
C       PHI,THETA,PSI                MOTIF ROTATION            (INPUT)
C
C FLOW:
C                                             FFT
C   motif (m) --------\                       v
C           x & pad > |-> padded-motif (mp) ----> motif-FFT (mp) -`
C   mask  (k) --------/                                           |
C                                                                 | 
C                                         '-----------------------'                |
C                                        v  
C                     FFT                 `------->\  
C                      v                  CMULCON  |-> vol-mult (mul)                  
C   search-vol (vol) -----> vol-FFT (vol) -------->/            |
C                                                               |
C                                                               | < FFT
C                                /- vol_lsd (mp)                |
C  peaks-doc <-- peaks (mul) <--| < PEAK_SEARCH_SCALE_LSD       |                                   
C                                \- vol-peak-norm (mul) <-------'                                   
C    
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE MOTIF(LUNMOT,LUNMSK,LUNVOLFFT,LUNLSD,LUNPKS,LUNDOC, 
     &                   NSAMM,NROWM,NSLICEM, NSAM,NROW,NSLICE, NPEAKS,
     &                   PHI,THETA,PSI)

#ifdef SP_LIBFFTW3
        USE TYPE_KINDS
        INTEGER(KIND=I_8)            :: IPLAN = 0    !STRUCTURE POINTER 
        INTEGER(KIND=I_8)            :: IPLANR = 0   !STRUCTURE POINTER 
#endif

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER, INTENT(IN)          :: LUNMSK,LUNMOT,LUNVOLFFT 
	INTEGER, INTENT(IN)          :: LUNLSD,LUNPKS,LUNDOC 
	INTEGER, INTENT(IN)          :: NSAMM,NROWM,NSLICEM 
	INTEGER, INTENT(IN)          :: NSAM,NROW,NSLICE
	INTEGER, INTENT(IN)          :: NPEAKS
	REAL, INTENT(IN)             :: PHI,THETA,PSI
 

        REAL                         :: SMALLEST,VAL,PEAKHTT
        INTEGER                      :: ISLICE,IROW,ISAM,LOCSMALL
        INTEGER                      :: IXT,IYT,IZT, MID,I

c       LOGICAL                      :: SPIDER_SIGN  = .FALSE.
c       LOGICAL                      :: SPIDER_SCALE = .FALSE.
        LOGICAL                      :: SPIDER_SIGN  = .TRUE.
        LOGICAL                      :: SPIDER_SCALE = .TRUE.

c       character (len=maxnam)       :: filnam

C       AUTOMATIC ARRAYS
	REAL                         :: PEAKHT(NPEAKS)
	INTEGER                      :: IPEAKX(NPEAKS)
	INTEGER                      :: IPEAKY(NPEAKS)
	INTEGER                      :: IPEAKZ(NPEAKS)
        INTEGER                      :: ILOCS(3), ILOCS1(1)
        REAL                         :: DLIST(7)

C       ALLOCATABLE ARRAYS
	REAL, ALLOCATABLE            :: BUFMP(:,:,:)
	REAL, ALLOCATABLE            :: BUFVOL(:,:,:)
	REAL, ALLOCATABLE            :: BUFMUL(:,:,:)
	REAL, ALLOCATABLE            :: BUFK(:),BUFM(:)


C       EXTRA COLUMN SPACE FOR FOURIER TRANSFORM
        LSE  = NSAM + 2 - MOD(NSAM,2)

C       MAKE BUFMP & BUFVOL FOR PADDED IMAGES WITH FOURIER ROW LENGTH
	ALLOCATE(BUFVOL(LSE,NROW,NSLICE), 
     &           BUFMP( LSE,NROW,NSLICE),
     &           BUFMUL(LSE,NROW,NSLICE),
     &           BUFM(NSAMM*NROWM*NSLICEM),
     &           BUFK(NSAMM*NROWM*NSLICEM),
     &           STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = 3 * LSE*NROW*NSLICE +  2 * NSAMM*NROWM*NSLICEM       
           CALL ERRT(46,'MOTIF; BUFVOL...',MWANT)
           GOTO 9999
        ENDIF

        !write(6,*) ' Motif size:          ',nsamm,nrowm,nslicem
        !write(6,*) ' Search volume size:  ',nsam,nrow,nslice

C       PROCESS MOTIF VOLUME SO THAT THE AVERAGE INSIDE 
C       THE MASK = 0 AND STANDARD DEVIATION INSIDE THE MASK = 1

C       READ SMALL MOTIF VOLUME 
        CALL REDVOL(LUNMOT, NSAMM,NROWM, 1,NSLICEM, BUFM,IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 9999

        !fmin = minval(bufm)
        !fmax = maxval(bufm)
        !write(6,*)' Small motif vol: ',fmin,'...',fmax

C       READ SMALL MOTIF MASK VOLUME  
        CALL REDVOL(LUNMSK, NSAMM,NROWM, 1,NSLICEM, BUFK,IRTFLG)
	IF (IRTFLG .NE. 0)  GOTO 9999

        !fmin = minval(bufk)
        !fmax = maxval(bufk)
        !write(6,*)' Small motif mask vol: ',fmin,'...',fmax

C       MULTIPLY MOTIF * MOTIF MASK
        NPIXM    = NSAMM * NROWM * NSLICEM  ! NUMBER OF MASK/MOTIF PIXELS

        SUMX40   = 0.0                      ! SUM OF MOTIF PIXELS
        SQSUMX45 = 0.0                      ! SUM OF MOTIF PIXELS SQUARED
        SUMX50   = 0.0                      ! SUM OF PIXELS IN MASK

c$omp   parallel do private(i),reduction(+:sumx40,sqsumx45,sumx50)
        DO I = 1,NPIXM
           IF (BUFK(I) > 0.5) THEN
              SUMX40   = SUMX40   + BUFM(I)     ! SUM OF MOTIF PIXELS
              SQSUMX45 = SQSUMX45 + BUFM(I)**2  ! SUM OF MOTIF PIXELS SQUAR
              SUMX50   = SUMX50   + 1           ! # OF PIXELS IN MASK
           ENDIF
       ENDDO

       IF (SUMX50 .LE. 0) THEN
           CALL ERRT(101,'NO NON-ZERO PIXELS INSIDE MASK',NE)
           GOTO 9999
        ENDIF
        !WRITE(NOUT,*) ' Non-zero pixels inside mask:', SUMX50
        !write(6,*)    ' non-zero pixels inside mask:', sumx50

C       AVG INSIDE MASK         
        AVGI = SUMX40 / SUMX50       ! SUM OF PIXELS INSIDE / # PIXELS  

C       SD INSIDE MASK        
        SDEV = SQRT((SQSUMX45- (SUMX40 * SUMX40 / SUMX50)) /(SUMX50 -1))
        VMULTX48 = 1.0 / SDEV             ! FOR SPEED

        !write(6,901)avgi, vmultx48
901     format('  avgi, vmultx48:              ',3(1pg13.2))

C       PAD ROTATED MOTIF INTO SIZE OF LARGE VOLUME, MASK, & NORMALIZE
c$omp   parallel do private(islice,irow,isam,iloc)
        DO       ISLICE = 1,NSLICEM
           DO    IROW   = 1,NROWM
              DO ISAM   = 1,NSAMM
                  ILOC = (ISLICE-1)*NSLICEM*NROWM + 
     &                   (IROW-1)*NROWM + ISAM

                  IF (BUFK(ILOC) < 0.5) THEN
                     BUFMP(ISAM,IROW,ISLICE) = 0.0
                  ELSE
                     BUFMP(ISAM,IROW,ISLICE) = 
     &                    (BUFM(ILOC) - AVGI) * VMULTX48
                  ENDIF
              ENDDO
 
C             FILL REMAINING PADDING COLS 
              BUFMP(NSAMM+1:LSE, IROW, ISLICE) = 0.0
           ENDDO

C          FILL REMAINING PADDING ROWS 
           BUFMP(1:LSE, NROWM+1:NROW, ISLICE) = 0.0
        ENDDO
c$omp   end parallel do

C       FILL REMAINING PADDING SLICES 
        BUFMP(1:LSE, 1:NROW, NSLICEM+1:NSLICE) = 0.0

        !fmin = minval(bufmp(1:nsam, 1:nrow, 1:nslice))
        !fmax = maxval(bufmp(1:nsam, 1:nrow, 1:nslice))
        !write(6,*)' Masked,scaled,padded motif vol: ',fmin,'...',fmax

C       FFT ON MASKED, SCALED, PADDED MOTIF VOLUME
        INV = +1
        CALL FMRS(BUFMP, NSAM,NROW,NSLICE, IPLAN,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(101,'MOTIF FFT ERROR ON BUFMP',NE)
           RETURN
        ENDIF

C       SET F(0,0) ELEMENT = ZERO FOR NORMALIZATION AS IN REAL SPACE 
        BUFMP(1,1,1) = 0.0
        BUFMP(2,1,1) = 0.0

C       READ FFT OF SEARCH VOLUME  
        CALL REDVOL(LUNVOLFFT, LSE,NROW, 1,NSLICE, BUFVOL,IRTFLG)

C       MULTIPLY SEARCH VOLUME FFT WITH COMPLEX CONJUGATE OF PADDED MASK
        NCOMP = LSE * NROW * NSLICE / 2
        CALL CMULCON(BUFVOL,BUFMP,BUFMUL,NCOMP)

C       INVERSE FFT ON MULTIPLIED SEARCH VOLUME
        INV = -1
        CALL FMRS(BUFMUL, NSAM,NROW,NSLICE, IPLANR,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(101,'MOTIF REVERSE FFT ERROR ON BUFMUL',NE)
           RETURN
        ENDIF

        !fmin = minval(bufmul(1:nsam, 1:nrow, 1:nslice))
        !fmax = maxval(bufmul(1:nsam, 1:nrow, 1:nslice))
        !write(6,*) ' fftd search*padded mask:',fmin,'...',fmax

C       LSD SCALING & PEAK LOCATION ----------------------------

C       READ PRECOMPUTED LSD IMAGE INTO PADDED SPACE OF BUFMP
        CALL READV(LUNLSD,BUFMP, LSE,NROW, NSAM,NROW,NSLICE)
	IF (IRTFLG .NE. 0)  GOTO 9999

        !fmin = minval(bufmp(1:nsam, 1:nrow, 1:nslice))
        !fmax = maxval(bufmp(1:nsam, 1:nrow, 1:nslice))
        !write(6,*)' lsd vol: ',fmin,'...',fmax

C       MULTIPLY BY CORRESPONDING ELEMENT OF LOCAL STANDARD DEV. VOLUME 
C       LSD VOLUME WAS DIVIDED BY NUMBER OF NON-ZERO PIXELS IN MASK         

        SMALLEST = -HUGE(SMALLEST)    ! HEIGHT SMALLEST PEAK
        LOCSMALL = 1                  ! LOCATION OF SMALLEST PEAK
        IGOT     = 0                  ! CURRENT NUMBER OF PEAKS
        PEAKHT   = SMALLEST           ! INITIALIZE PEAK ARRAY

        DO       ISLICE=1,NSLICE
           DO    IROW  =1,NROW
              DO ISAM  =1,NSAM

                VAL = BUFMUL(ISAM,IROW,ISLICE) * BUFMP(ISAM,IROW,ISLICE)
 
                BUFMUL(ISAM,IROW,ISLICE) = VAL

                IF (VAL .GT. SMALLEST) THEN
C                  REPLACE CURRENT SMALLEST DATA
                   !write(6,*) 'val,smallst,igot:',val,smallest,igot

                   PEAKHT(LOCSMALL) = VAL
                   IPEAKX(LOCSMALL) = ISAM
                   IPEAKY(LOCSMALL) = IROW
                   IPEAKZ(LOCSMALL) = ISLICE
                   IF (IGOT .LT. NPEAKS) IGOT = IGOT + 1

C                  UPDATE CURRENT SMALLEST POINTER
                   ILOCS1    = MINLOC(PEAKHT)
                   LOCSMALL  = ILOCS1(1)
                   SMALLEST  = PEAKHT(LOCSMALL)

                ENDIF      ! END OF: IF (VAL .GT. SMALLEST)
             ENDDO         ! END OF: DO ISAM
           ENDDO           ! END OF: DO IROW
        ENDDO              ! END OF: DO ISLICE

C      OUTPUT SCALED PEAK VOLUME NSAM x NROW x NSLICE     
       CALL WRITEV(LUNPKS, BUFMUL, LSE,NROW, NSAM,NROW,NSLICE)

       IF (VERBOSE) THEN
          FMIN = MINVAL(BUFMUL(1:NSAM, 1:NROW, 1:NSLICE))
          FMAX = MAXVAL(BUFMUL(1:NSAM, 1:NROW, 1:NSLICE))
          WRITE(NOUT,*) ' Scaled peak vol range:',FMIN,'...',FMAX
          !WRITE(6,*)    ' Scaled peak vol range:',FMIN,'...',FMAX
       ENDIF

c      do i = 1,npeaks
c        write(6,91)i,IPEAKX(I),IPEAKY(I),IPEAKZ(I),peakht(i)
c      enddo

C      PEAK SORTING AND OUTPUT ----------------------------

C      SORT THE PEAKS BY HEIGHT
       CALL SORTR3I(PEAKHT,IPEAKX,IPEAKY,IPEAKZ,IGOT)

C      SAVE PEAKS IN DOC FILE
       DLIST(1) = PHI
       DLIST(2) = THETA
       DLIST(3) = PSI

       DO I = IGOT,1,-1
           IKEY     = IGOT - I + 1

           DLIST(4) = IPEAKX(I)
           DLIST(5) = IPEAKY(I)
           DLIST(6) = IPEAKZ(I)

C          DON'T WRITE CCC IN SCIENTIFIC FORMAT. UNIX 'sort' COMMAND 
C          WILL FAIL.    CCC < 0.05 DOES NOT MEAN MUCH
           IF (PEAKHT(I) < 0.05) PEAKHT(I) = -99.0 
           DLIST(7) = PEAKHT(I)

           CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,7,IRTFLG)

#ifdef DEBUG
           write(6,91)ikey,ipeakx(i),ipeaky(i),ipeakz(i),peakht(i)
91         format(' ',i3,'(',i5,',',i5,',',i5,'):',1pg13.4)
#endif

       ENDDO

 
9999   IF (ALLOCATED(BUFVOL)) DEALLOCATE(BUFVOL)
       IF (ALLOCATED(BUFMUL)) DEALLOCATE(BUFMUL)
       IF (ALLOCATED(BUFMP))  DEALLOCATE(BUFMP)
       IF (ALLOCATED(BUFM))   DEALLOCATE(BUFM)
       IF (ALLOCATED(BUFK))   DEALLOCATE(BUFK)

       CALL FLUSHRESULTS()

       END



#ifdef NEVER
C++*********************************************************************
C
C  PEAKLOC.F                                WRITTEN JUL 08 ARDEAN LEITH
C
C **********************************************************************
C
C  PURPOSE: SEARCH FOR NPEAKS PEAKS IN BUFPK SCALED BY LOCAL STD. DEV. 
C
C PARAMETERS:
C       BUFPK        PEAK ARRAY                         (INPUT)
C       BUFLSD       LSD  ARRAY                         (INPUT)
C       NPEAKS       PEAK ARRAY LENGTH                  (INPUT)
C       IGOT         PEAKS FOUND                        (OUTPUT)
C       IX,IY,IZ     PEAK LOCATION ARRAYS               (OUTPUT)
C       PEAKHT       PEAK ARRAY                         (OUTPUT)
C       IRTFLG       ERROR FLAG                         (OUTPUT)
C
C
C **********************************************************************

       SUBROUTINE PEAKLOC(BUFPK,BUFLSD,
     &                    NSAMP,NROWP,NSLICEP, NSAM,NROW,NSLICE,
     &                    NPEAKS, IGOT,PEAKHT,IX,IY,IZ, IRTFLG)

       IMPLICIT NONE

       REAL, INTENT(IN)     :: BUFPK(NSAMP,NROWP,NSLICEP)
       REAL, INTENT(IN)     :: BUFLSD(NSAMP,NROWP,NSLICEP)
       INTEGER, INTENT(IN)  :: NSAMP,NROWP,NSLICEP 
       INTEGER, INTENT(IN)  :: NSAM,NROW,NSLICE 
       INTEGER, INTENT(IN)  :: NPEAKS 
       INTEGER, INTENT(OUT) :: IGOT 

       REAL, INTENT(OUT)    :: PEAKHT(NPEAKS)
       INTEGER, INTENT(OUT) :: IX(NPEAKS),IY(NPEAKS),IZ(NPEAKS)
       INTEGER, INTENT(OUT) :: IRTFLG 

       REAL                 :: SMALLEST,VAL,PEAKHTT
       INTEGER              :: ISLICE,IROW,ISAM,LOCSMALL
       INTEGER              :: ILOC(1)
       INTEGER              :: IXT,IYT,IZT, MID,I
    
       SMALLEST = -HUGE(SMALLEST)    ! HEIGHT SMALLEST PEAK
       LOCSMALL = 1                  ! LOCATION OF SMALLEST PEAK
       IGOT     = 0                  ! CURRENT NUMBER OF PEAKS
       PEAKHT   = SMALLEST           ! INITIALIZE PEAK ARRAY

       DO ISLICE = 1,NSLICE
          DO IROW = 1,NROW
             DO ISAM = 1,NSAM
C               FIND PEAK VALUE
ccc                VAL = BUFPK(ISAM,IROW,ISLICE) * BUFLSD(ISAM,IROW,ISLICE)
                VAL = BUFPK(ISAM,IROW,ISLICE)

                IF (VAL .GT. SMALLEST) THEN
C                  REPLACE CURRENT SMALLEST DATA
c                  write(6,*) 'val,smallst,igot:',val,smallest,igot

                   PEAKHT(LOCSMALL) = VAL
                   IX(LOCSMALL)     = ISAM
                   IY(LOCSMALL)     = IROW
                   IZ(LOCSMALL)     = ISLICE
                   IF (IGOT .LT. NPEAKS) IGOT = IGOT + 1

C                  UPDATE CURRENT SMALLEST POINTER
                   ILOC      = MINLOC(PEAKHT)
                   LOCSMALL  = ILOC(1)
                   SMALLEST  = PEAKHT(LOCSMALL)
                ENDIF
            ENDDO
         ENDDO
       ENDDO

C      SORT THE PEAKS
       CALL SORTR3I(PEAKHT,IX,IY,IZ,IGOT)

C      REVERSE THE SORT
       MID = IGOT / 2
       DO I=1,MID
          IXT          = IX(I)
          IX(I)        = IX(IGOT-I+1) 
          IX(IGOT-I+1) = IXT
          IYT          = IY(I)
          IY(I)        = IY(IGOT-I+1) 
          IY(IGOT-I+1) = IYT
          IZT          = IZ(I)
          IZ(I)        = IZ(IGOT-I+1) 
          IZ(IGOT-I+1) = IZT
          PEAKHTT      = PEAKHT(I)
          PEAKHT(I)    = PEAKHT(IGOT-I+1) 
          PEAKHT(IGOT-I+1) = PEAKHTT
       ENDDO
        
       IRTFLG = 0
       IF (IGOT < NPEAKS) IRTFLG = 1

       END        


      !write(6,*) ' Pixel(3,1,1),(4,1,1):',bufvol(3,1,1),bufvol(4,1,1)
      !write(6,*) ' Pixel(100,1,1),(101,1,1):',
      !&             bufvol(100,1,1),bufvol(101,1,1)

      maxim = 0
      itype = 3
      call opfilec(0,.false.,'jnkpadedmotif',51,'U',itype,
     &            lse,nrow,nslice, maxim,' ',.false.,irtflg)
      call wrtvol(51,lse,nrow,1,nslice,bufmp,irtflg)
      close(51)

        CALL PEAKLOC(BUFMUL,BUFMP,
     &             LSE,NROW,NSLICE, NSAM,NROW,NSLICE,
     &             NPEAKS,PEAKHT,IPEAKX,IPEAKY,IPEAKZ, IRTFLG)

      do i = 1,npeaks
        write(6,91)IPEAKX(I),IPEAKY(I),IPEAKZ(I),peakht(i)
91      format(' peak(',i5,',',i5,',',i5,'):',1pg13.4)
      enddo


      npeaks = 10
      do i = 1,npeaks
        ilocs = maxloc(bufmul)
        pdum  = bufmul(ilocs(1),ilocs(2),ilocs(3))

        write(6,91)pdum, ilocs(1),ilocs(2),ilocs(3)
        bufmul(ilocs(1),ilocs(2),ilocs(3)) = fmin
      enddo

c        ixt   = isx - ilocs(1) + 2
c        iyt   = isy - ilocs(2) + 2
c        izt   = isz - ilocs(3) + 2
c        if (ixt .le. 0) ixt = ixt + nsam   
c        if (iyt .le. 0) iyt = iyt + nrow  
c        if (izt .le. 0) izt = izt + nslice  
c        write(6,91)pdum, ixt,iyt,izt
c        write(6,*) ' ---------------- '


#endif

