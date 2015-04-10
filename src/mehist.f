
C++*********************************************************************
C
C    MEHIST.F            NEW JUNE 05 ArDean Leith
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
C    MEHIST(LUN1,LUN2,NSAM,NROW,HMIN,HMAX)
C
C    PURPOSE: MAX. ENTROPY HISTOGRAM THRESHOLDING BASED ON THE 
C             ENTOPY OF HISTOGRAM.
C
C    REFERENCE: P.K. SAHOO, S. SOLTANI, K.C. WONG &, Y.C. CHEN 
C               "A SURVEY OF THRESHOLDING TECHNIQUES", 
C               COMPUTER VISION, GRAPHICS, AND IMAGE PROCESSING, 
C               VOL. 41, PP.233-260, 1988.
C
C    PARAMETERS:
C        LUN1      LOGICAL UNIT NUMBER OF INPUT FILE
C        LUN2      LOGICAL UNIT NUMBER OF OUTPUT FILE
C        NSAM      NUMBER OF SAMPLES
C        NROW      NUMBER OF ROWS
C        NSLICE    NUMBER OF SLICES
C
C--*********************************************************************

      SUBROUTINE MEHIST(LUN1,LUN2,NSAM,NROW,NSLICE,HMIN,HMAX)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      REAL, DIMENSION(NSAM*NROW)      :: REDBUF
      INTEGER, PARAMETER              :: NBINS = 256
      REAL, DIMENSION(NBINS)          :: FREQ,FNORMHIST,PT,HB,HW

C     CALCULATE MAXIMUM ENTROPY SPLIT OF A HISTOGRAM.

C     ZERO THE HISTOGRAM FREQUENCIES
      FREQ  = 0.0 

      HDIFF = HMAX - HMIN
      FF    = (NBINS - 1.0) / HDIFF

C     FIND HISTOGRAM OF IMAGE, PLACE HISTOGRAM IN FREQ 
      DO  ISLICE=1,NSLICE
         CALL REDVOL(LUN1,NSAM,NROW,ISLICE,ISLICE,REDBUF,IRTFLG)
         DO  ISAM = 1,NSAM*NROW
C           HISTOGRAM THIS PIXEL 

C           FIND BIN NUMBER
            JBIN = INT((REDBUF(ISAM) - HMIN) * FF) + 1.5

            IF (JBIN .GE. 1 .AND. JBIN .LE. NBINS) THEN
C              WITHIN HISTOGRAM RANGE
               FREQ(JBIN) = FREQ(JBIN) + 1.0
            ENDIF
         ENDDO
      ENDDO
  
C     NORMALIZE HISTOGRAM, MAKE THE SUM OF ALL BINS EQUAL TO 1.
      FNUMPIX   = 1.0 / (NSAM * NROW * NSLICE)
      FNORMHIST = FREQ * FNUMPIX

      PT(1) = FNORMHIST(1)
      DO  IBIN=2,NBINS
         PT(IBIN) = PT(IBIN - 1) + FNORMHIST(IBIN) 
      ENDDO

C     ENTROPY FOR BLACK AND WHITE PARTS OF THE HISTOGRAM
      EPSILON = TINY(EPSILON)

      DO IBIN=1,NBINS

        IF (PT(IBIN) > EPSILON) THEN 
C         BLACK ENTROPY
          HHB = 0.0
          DO I = 1,IBIN
             IF (FNORMHIST(I) > EPSILON) THEN 
                HHB = HHB -
     &          FNORMHIST(I)/ PT(IBIN) * LOG10(FNORMHIST(I)/PT(IBIN))
              ENDIF
           ENDDO
           HB(IBIN) = HHB
        ELSE
           HB(IBIN) = 0.0
        ENDIF

C       WHITE  ENTROPY
        PTW = 1 - PT(IBIN)
        IF (PTW > EPSILON) THEN
            HHW = 0.0
            DO I = IBIN+1,NBINS
               IF (FNORMHIST(I) > EPSILON) THEN
                  HHW = HHW -
     &               FNORMHIST(I) / PTW * LOG10(FNORMHIST(I) / PTW)
               ENDIF
            ENDDO
            HW(IBIN) = HHW
         ELSE  
            HW(IBIN) = 0.0
         ENDIF
      ENDDO       

C     FIND HISTOGRAM INDEX WITH MAXIMUM ENTROPY
      VALMAX  = HB(1) + HW(1)
      IBINMAX = 1

      DO IBIN = 2,NBINS
        VAL = HB(IBIN) + HW(IBIN)
        IF (VAL > VALMAX) THEN
          VALMAX  = VAL
          IBINMAX = IBIN
        ENDIF
      ENDDO

      THRESH = (IBINMAX-1) / FF + HMIN

      WRITE(NOUT,91) IBINMAX,FF
91    FORMAT('  IBINMAX:',I10,'  FF:',G10.3)

      WRITE(NOUT,90) HMIN,HMAX,THRESH
90    FORMAT('  RANGE:',G10.3,'...',G10.3,'  THRESHOLDED AT:',G10.3)

C     THRESHOLD IMAGE, PLACE IT REDBUF 
      DO  ISLICE=1,NSLICE
         CALL REDVOL(LUN1,NSAM,NROW,ISLICE,ISLICE,REDBUF,IRTFLG)
         DO  ISAM = 1,NSAM*NROW
C           THRESHOLD THIS PIXEL 

            IF (REDBUF(ISAM) .LT. THRESH) THEN
               REDBUF(ISAM) = 0.0
            ELSE
               REDBUF(ISAM) = 1.0
            ENDIF
         ENDDO
         CALL WRTVOL(LUN2,NSAM,NROW,ISLICE,ISLICE,REDBUF,IRTFLG)
      ENDDO
     
      END
