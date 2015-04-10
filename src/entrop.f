
C++*********************************************************************
C
C  ENTROP.F              NEW                     FEB  2003  ArDean Leith
C                        VMIN, VMAX              FEB  2004  ArDean Leith
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
C    ENTROP(LUN,IFORM,NSAM,NROW,NSLICE,ENTROPY,IRTFLG)
C
C    PURPOSE:    COMPUTE ENTROPY OF IMAGE FROM HISTOGRAM 
C
C    PARAMETERS:  LUN        IO UNIT NUMBER OF IMAGE FILE
C                 NSAM,NROW  DIMENSIONS OF IMAGE              
C                 NSLICE     DIMENSIONS OF IMAGE
C                 ENTROPY    ENTROPY                              (RET.)
C                 IRTFLG                                          (RET.)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE ENTROP(LUN,NSAM,NROW,NSLICE,ENTROPY,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      COMMON     FREQ
      REAL, DIMENSION(5000) :: FREQ

      DOUBLE PRECISION      :: FNUMEL

      COMMON /IOBUF/ REDBUF(NBUFSIZ)

C     MAKE SURE STATISTICS ARE CURRENT
      IF (IMAMI .NE. 1) CALL NORM3(LUN,NSAM,NROW,NSLICE,FMAX,FMIN,AV)

C     NBINS SHOULD BE <= 5000
      NBINS  = 256
      CALL RDPRI1S(NBINS,NOT_USED,'NUMBER OF BINS',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      IF (NBINS .GT. 5000) THEN
         CALL ERRT(101,'MAX NUMBER OF BINS =5000',NE)
         RETURN
      ENDIF

      WRITE(NOUT,*) ' IMAGE RANGE: ',FMIN,'.....',FMAX
      VMIN = FMIN
      VMAX = FMAX
      CALL RDPRM2S(VMIN,VMAX,NOT_USED,'HISTOGRAM RANGE',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
 
C     ZERO THE HISTOGRAM FREQUENCIES
      DO  K = 1,NBINS
         FREQ(K) = 0.0 
      ENDDO

      FNPIX  = NSAM * NROW * NSLICE
      HDIFF  = VMAX - VMIN
      IF (HDIFF .EQ. 0.0) THEN
         CALL ERRT(101,'BLANK IMAGE',NE)
         IRTFLG = 1
         RETURN
      ENDIF

      FF     = (NBINS - 1.0) / HDIFF
      FINPIX  = 0

      DO  IREC=1, NSLICE * NROW
         CALL REDLIN(LUN,REDBUF,NSAM,IREC)

         DO  ISAM = 1,NSAM
C           HISTOGRAM THIS PIXEL 

C           FIND BIN NUMBER
            BVAL = REDBUF(ISAM)
            JBIN = INT((BVAL - VMIN) * FF) + 1.5

            IF (JBIN.GE.1 .AND. JBIN.LE.NBINS) THEN
C              WITHIN HISTOGRAM RANGE
               FREQ(JBIN) = FREQ(JBIN) + 1.0
               FINPIX      = FINPIX + 1
           ENDIF
        ENDDO
      ENDDO

C     FIND ENTROPY
      ENTROPY = 0.0

      DO  IBIN = 1,NBINS
         FT      = FREQ(IBIN) / FNPIX
         IF (FT .GT. 0.0) ENTROPY = ENTROPY - FT * LOG(FT) 
      ENDDO

      WRITE(NOUT,*) ' '

      IBIG = HUGE(IBIG)

      IF (FNPIX .LT. IBIG) THEN
         NPIX  = NSAM * NROW * NSLICE
         INPIX = FINPIX
         WRITE(NOUT,90) FMIN,FMAX,VMIN,VMAX,INPIX,NPIX,NBINS,ENTROPY

90       FORMAT(
     &    ' FILE RANGE:        ',1PG11.4,'   .........     ',1PG11.4,/,
     &    ' HIST. RANGE:       ',1PG11.4,'   .........     ',1PG11.4,/,
     &    ' IMAGE PIXELS:      ',I11,    '   HIST. PIXELS: ',I11,/, 
     &    ' NO. OF BINS:       ',I11,/,
     &    ' ENTROPY:           ',1PG11.4/)
      ELSE

         WRITE(NOUT,91) FMIN,FMAX,VMIN,VMAX,FINPIX,FNPIX,NBINS,ENTROPY
91       FORMAT(
     &    ' FILE RANGE:        ',1PG11.4,'   .........     ',1PG11.4,/,
     &    ' HIST. RANGE:       ',1PG11.4,'   .........     ',1PG11.4,/,
     &    ' IMAGE PIXELS:      ',1PG11.4,'   HIST. PIXELS: ',1PG11.4,/, 
     &    ' NO. OF BINS:       ',I11,/,
     &    ' ENTROPY:           ',1PG11.4/)
      ENDIF

      CALL REG_GET_USED(NSEL_USED)

      IF (NSEL_USED .GT. 0) THEN
         CALL REG_SET_NSEL(1,1,ENTROPY, 0.0, 0.0, 0.0, 0.0,IRTFLG)
      ENDIF

      RETURN
      END

