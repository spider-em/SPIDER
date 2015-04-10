
C++*********************************************************************
C
C JOHIST.F                 NEW                    MAR  2003 ArDean Leith
C                          ENTROPY2 REG BUG       APR  2003 ArDean Leith
C                          VMIN, VMAX             FEB  2004 ArDean Leith
C                          DBLE(NPIX)             DEC  2005 ArDean Leith
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
C    JOHIST(LUN1,LUN2,NSAM,NROW,NSLICE,NBINS, 
C           FMIN1,FMAX1,FMIN2,FMAX2,IRTFLG)
C
C    PURPOSE:    MUTUAL SHARED INFO FOR 2 INPUT IMAGES.
C
C    PARAMETERS:  LUN1       IO UNIT NUMBER OF IMAGE FILE
C                 LUN2       IO UNIT NUMBER OF IMAGE FILE
C                 NSAM,NROW  DIMENSIONS OF IMAGE              
C                 NSLICE     DIMENSIONS OF IMAGE
C                 FMIN,FMAX  IMAGES MIN & MAX                    (SENT)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE JOHIST(LUN1,LUN2,NSAM,NROW,NSLICE,NBINS,
     &                 FMIN1,FMAX1,FMIN2,FMAX2,VMIN,VMAX,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

#ifndef SP_32
      INTEGER *8 NPIX
      INTEGER *8 FREQ1(NBINS),FREQ2(NBINS)
      INTEGER* 8 FREQ(NBINS,NBINS) 
#else
      INTEGER *4 NPIX
      INTEGER *4 FREQ1(NBINS),FREQ2(NBINS)
      INTEGER* 4 FREQ(NBINS,NBINS) 
#endif

      REAL, DIMENSION(NSAM) :: REDBUF1,REDBUF2
      DOUBLE PRECISION      :: FNPIX

C     ZERO THE HISTOGRAM FREQUENCIES
      FREQ  = 0 
      FREQ1 = 0 
      FREQ2 = 0 

      HDIFF = VMAX - VMIN
      FF    = (NBINS - 1.0) / HDIFF

      NPIX   = NSAM * NROW * NSLICE
      FNPIX  = DBLE(1.0) / DBLE(NPIX)
      NREC   = NSLICE * NROW

C     GET HISTOGRAMS FROM IMAGE VALUES

      NPIX1  = 0
      NPIX2  = 0
      NPIX12 = 0

      DO  IREC=1,NREC
         CALL REDLIN(LUN1,REDBUF1,NSAM,IREC)
         CALL REDLIN(LUN2,REDBUF2,NSAM,IREC)

         DO  ISAM = 1,NSAM
C           HISTOGRAM THIS PIXEL 

C           FIND BIN1 NUMBER
            BVAL  = REDBUF1(ISAM)
            JBIN1 = INT((BVAL - VMIN) * FF) + 1.5

C           FIND BIN2 NUMBER
            BVAL  = REDBUF2(ISAM)
            JBIN2 = INT((BVAL - VMIN) * FF) + 1.5

            IF (JBIN1.GE.1 .AND. JBIN1.LE.NBINS) NPIX1 = NPIX1 + 1.0

            IF (JBIN2.GE.1 .AND. JBIN2.LE.NBINS) NPIX2 = NPIX2 + 1.0

            IF (JBIN1.GE.1 .AND. JBIN1.LE.NBINS  .AND.
     &          JBIN2.GE.1 .AND. JBIN2.LE.NBINS) THEN
C              WITHIN JOINT HISTOGRAM RANGE
               NPIX12 = NPIX12 + 1
               FREQ1(JBIN1)      = FREQ1(JBIN1) + 1
               FREQ2(JBIN2)      = FREQ2(JBIN2) + 1
               FREQ(JBIN1,JBIN2) = FREQ(JBIN1,JBIN2) + 1
            ENDIF
        ENDDO
      ENDDO

C     FIND ENTROPY
      ENTROPY1 = 0.0
      DO IBIN = 1,NBINS
         FT      = FREQ1(IBIN) * FNPIX
         IF (FT .GT. 0.0) ENTROPY1 = ENTROPY1 - FT * LOG(FT) 
      ENDDO

      ENTROPY2 = 0.0
      DO IBIN = 1,NBINS
         FT      = FREQ2(IBIN) * FNPIX
         IF (FT .GT. 0.0) ENTROPY2 = ENTROPY2 - FT * LOG(FT) 
      ENDDO

      ENTROPY = 0.0
      DO IBIN = 1,NBINS
         DO JBIN=1,NBINS
            FT = FREQ(IBIN,JBIN) * FNPIX
            IF (FT .GT. 0.0) ENTROPY = ENTROPY - FT * LOG(FT)
         ENDDO 
      ENDDO

      FMSI = ENTROPY1 + ENTROPY2 - ENTROPY

      WRITE(NOUT,*) ' '
      WRITE(NOUT,90) FMIN1,FMAX1,FMIN2,FMAX2,VMIN,VMAX,NPIX,NBINS,
     &               NPIX1,NPIX2,NPIX12ENTROPY1,ENTROPY2,ENTROPY,FMSI

90    FORMAT(
     &    ' FIRST FILE RANGE:  ',1PG11.4,'   .........      ',1PG11.4,/,
     &    ' 2ND   FILE RANGE:  ',1PG11.4,'   .........      ',1PG11.4,/,
     &    ' HISTOGRAM RANGE:   ',1PG11.4,'   .........      ',1PG11.4,/,
     &    ' PIXELS:            ',I11,    '  NO. OF BINS:    ',I11,/,
     &    ' FIRST FILE PIXELS: ',I11,    '  2ND FILE PIXELS:',I11,/,
     &    ' BOTH  FILE PIXELS: ',I11,/,
     &    ' FIRST FILE ENTROPY:  ',1PG11.4,/,
     &    ' 2ND   FILE ENTROPY:  ',1PG11.4,/,
     &    ' JOINT      ENTROPY:  ',1PG11.4,/,
     &    ' MUTUAL SHARED INFO:  ',1PG11.4/)

      CALL REG_GET_USED(NSEL_USED)

      IF (NSEL_USED .GT. 0) THEN
        CALL REG_SET_NSEL(1,4,FMSI,ENTROPY1,ENTROPY2,ENTROPY,0.0,IRTFLG)
      ENDIF

      RETURN
      END
  
