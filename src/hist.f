
C++*********************************************************************
C
C HIST.F                   FMINT FOR HIST RANGE  JUN  2000  ArDean Leith
C                          NaN OK                JAN  2001  ArDean Leith
C                          OLD MODE BUG          JAN  2002  ArDean Leith
C                          FNUMEL OVERFLOW BUG   JUL  2002  ArDean Leith
C                          INTEGER OUTPUT        OCT  2002  ArDean Leith
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
C    HIST(LUN,LUNMA,LUNDOC,NX,NY,NZ,IDUST,HMIN,HMAX)
C
C    PURPOSE:    COMPUTE 128 PLACE HISTOGRAM FROM IMAGE RECORDS, 
C                DISPLAY HISTOGRAM ON LINE PRINTER, TERMINAL, OR
C                IN DOC. FILE AND, OPTIONALLY, REMOVE DATA THAT ARE 
C                OUT OF A SPECIFIED STATISTICAL RANGE.
C
C    PARAMETERS:  LUN        IO UNIT NUMBER OF IMAGE FILE
C                 LUNDMA     IO UNIT NUMBER FOR MASK FILE
C                 LUNDOC     IO UNIT NUMBER FOR DOCUMENT FILE
C                 NX,NY  DIMENSIONS OF IMAGE              
C                 NZ     DIMENSIONS OF IMAGE
C                 HMIN,HMAX  HIST. MIN & MAX                    (RET.)
C                 HSIG,HMODE HIST. S.D. & MODE (USED BY DUST)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE HIST(LUN,LUNMA,LUNDOC,NX,NY,NZ,
     &                HMIN,HMAX,HSIG,HMODE)

      INCLUDE 'CMBLOCK.INC'

      DIMENSION        :: FREQ(128)
      DIMENSION        :: PLIST(4)
      DOUBLE PRECISION :: HAV, HAV2, DTOP, FNUMEL
      CHARACTER *1     :: ANS
      LOGICAL          :: MASK
      CHARACTER *1     :: NULL = CHAR(0)

      IF (IFORM .LE. 0) THEN
C        CAN NOT DO HISTOGRAM ON FOURIER FILES
         CALL ERRT(101,'CAN NOT HISTOGRAM FOURIER FILES',NE)
         RETURN
      ENDIF

C     MAKE SURE STATISTICS ARE CURRENT
      IF (IMAMI .NE. 1) CALL NORM3(LUN,NX,NY,NZ,FMAX,FMIN,AV)
      HMAX = FMAX
      HMIN = FMIN

C     NBINS MUST BE <= 1000!!!  why?? 
      NBINS  = 128
      NDEV   = NDAT
      IDEV   = 1
      IF (FCHAR(4:4) == 'R' .OR. FCHAR(5:5) == 'R') THEN
C        ONE OF THE LETTERS IS 'R' ASK USER FOR HISTOGRAM RANGE
         CALL RDPRM2(HMIN,HMAX,NOT_USED,'HISTOGRAM RANGE (MIN, MAX)')

      ELSEIF (FCHAR(4:4) == 'M') THEN
C        ASK USER FOR HISTOGRAM SINK
         CALL RDPRMC(ANS,NLET,.TRUE.,
     &      'OUTPUT TO RESULTS FILE, DOC. FILE, OR TERMINAL? (R/D/T)',
     &       NULL,IRTFLG)
      ENDIF

      IF  (FCHAR(4:4) == 'T' .OR. FCHAR(5:5) == 'T' .OR.
     &    (FCHAR(4:4) == 'M' .AND. ANS == 'T')) THEN
C         OUTPUT TO TERMINAL, NOT RESULTS FILE
          NBINS = 70
          NDEV  = NOUT
          IDEV  = 0
      ENDIF

C     ZERO THE HISTOGRAM FREQUENCIES
      DO  K = 1,NBINS
         FREQ(K) = 0.0 
      ENDDO

C     MPI 3-MAR-80 ADDED TO AVOID "BEAT" PHENOMENON (REMOVED AUG 98)   
C     APPEARED TO EXPAND RANGE OF HIST. MIN/MAX FOR < 0...2 

      MASK = .FALSE.
      NREC = NZ * NY
      IF (FCHAR(4:4) == 'M') THEN
C        HISTOGRAM UNDER MASK ONLY
         MASK = .TRUE.
C        DETERMINE RANGE UNDER MASK
         CALL HISTMINMAX(LUN,LUNMA,NX,NREC,HMIN,HMAX)
      ENDIF

      HDIFF  = HMAX - HMIN
      FF     = (NBINS - 1.0) / HDIFF
      BSIZ   = 1.0 / FF
      FNUMEL = 0.0
      HAV    = 0.0
      HAV2   = 0.0

      CALL HISTIMAGE(LUN,LUNMA,NX,NREC,FREQ,
     &               FNUMEL,HAV,HAV2,HMIN,FF,NBINS,MASK)

      IF (FCHAR(4:4) == 'D' .OR. 
     &    FCHAR(5:5) == 'D' .OR.
     &   (FCHAR(4:4) == 'M' .AND. ANS == 'D')) THEN
C        OUTPUT TO DOC FILE
         PLIST(2) = HMIN
         DO IBIN     =1,NBINS
            PLIST(1) = IBIN
            PLIST(2) = HMIN + (BSIZ *(IBIN - 1))
            PLIST(3) = FREQ(IBIN)
            CALL SAVD(LUNDOC,PLIST,3,IRTFLG)
         ENDDO
         CALL SAVDC
         CLOSE(LUNDOC)

      ELSEIF (FCHAR(1:2) .NE. 'DU') THEN
C        OUTPUT TO RESULTS FILE OR TERMINAL 

         WRITE(NDEV,*) ' '
C        GRAPHX IN GRAPHS USES FMIN & FMAX FOR HISTOGRAM LABELS
         FMINT = FMIN
         FMIN  = HMIN
         FMAXT = FMAX
         FMAX  = HMAX
         CALL GRAPHS(NDEV,FREQ,NBINS,1,IDEV,1.0,IRTFLG)
         FMIN = FMINT
         FMAX = FMAXT
         WRITE(NDEV,*) ' '
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(100,'HIST',NE)
            RETURN
         ENDIF
      ENDIF

C     FIND MAXIMUM FREQUENCY OCCURING IN HISTOGRAM (HISMAX) & LOCATION
      HISMAX = FREQ(1)
      MAXBIN = 1
      DO  IBIN = 2,NBINS
         IF (FREQ(IBIN) .GE. HISMAX) THEN
            HISMAX = FREQ(IBIN)
            MAXBIN = IBIN
         ENDIF
      ENDDO

C     CONVERT LOCATION MAXBIN TO CORRESPONDING IMAGE INTENSITY (MODE)
      IF (MAXBIN == 1) THEN
C        SUB-BIN ESTIMATE OF MODE
         BMODE  = 0.5

      ELSEIF (MAXBIN == NBINS) THEN
C        SUB-BIN ESTIMATE OF MODE 
C        (MAYBE THIS SHOULD AVERAGE ALL BINS AT END WITH SAME NUMBER?)
         BMODE  = FLOAT(NBINS) - 0.5
      ELSE
C        SUB-BIN ESTIMATE OF MODE
         YM1    = FREQ(MAXBIN-1)
         YP1    = FREQ(MAXBIN+1)
         BMODE  = FLOAT(MAXBIN-1) + (YM1-YP1)*0.5/ (YM1+YP1-2.*HISMAX)
      ENDIF
      HMODE = HMIN + BMODE * BSIZ

      DTOP  = HAV2 - HAV * HAV / FNUMEL
      IF (DTOP .LT. 0.0) THEN
C        SQRT OF NEGATIVE NUMBER
         WRITE(NOUT,*) '*** in HIST, SQRT(',DTOP,') IMPOSSIBLE' 
         CALL ERRT(100,'HIST',NE)
         RETURN

      ELSEIF (FNUMEL == 1) THEN
         WRITE(NOUT,*) '*** CAN NOT DETERMINE S.D. (ONLY ONE PIXEL) ' 
         CALL ERRT(100,'HIST',NE)
         RETURN
      ENDIF

      HAV    = HAV  / FNUMEL
      HSIG   = DSQRT(DTOP / (FNUMEL - 1))
      FNPIX  = NX * NY * NZ
      FNBINS = NBINS

      WRITE(NDEV,*) ' '

      IBIG = HUGE(IBIG)
      IF (FNPIX .LT. IBIG) THEN
         NPIX  = NX * NY * NZ
         NUMEL = FNUMEL
         WRITE(NDEV,90) FMIN,FMAX,HMIN,HMAX,NPIX,NUMEL,
     &               NBINS,BSIZ,HAV,HMODE,HSIG

90       FORMAT(
     &    '  FILE RANGE:        ',1PG11.4,'   .........     ',1PG11.4,/,
     &    '  HISTOGRAM RANGE:   ',1PG11.4,'   .........     ',1PG11.4,/,
     &    '  TOTAL PIXELS:      ',I11,    '   PIXELS IN HIST.: ',I11,/,
     &    '  NO. OF BINS:       ',I11,    '   BIN SIZE:     ',1PG11.4,/,
     &    '  HIST. MEAN:        ',1PG11.4,'   HIST. MODE:   ',1PG11.4,/,
     &    '  HIST. S.D.:        ',1PG11.4/)
      ELSE
        WRITE(NDEV,91) FMIN,FMAX,HMIN,HMAX,FNPIX,FNUMEL,
     &               FNBINS,BSIZ,HAV,HMODE,HSIG

91      FORMAT(
     &  '  FILE RANGE:        ',1PG11.4,'   .........     ',1PG11.4,/,
     &  '  HISTOGRAM RANGE:   ',1PG11.4,'   .........     ',1PG11.4,/,
     &  '  TOTAL PIXELS:      ',1PG11.4,'   PIXELS IN HIST.: ',1PG11.4,/,
     &  '  NO. OF BINS:       ',1PG11.4,'   BIN SIZE:     ',1PG11.4,/,
     &  '  HIST. MEAN:        ',1PG11.4,'   HIST. MODE:   ',1PG11.4,/,
     &  '  HIST. S.D.:        ',1PG11.4/)
      ENDIF

      WRITE(NDEV,*) ' '

      IF (NDEV .NE. NOUT) THEN
         WRITE(NOUT,*) ' '
         IF (FNPIX .LT. IBIG) THEN
            WRITE(NOUT,90) FMIN,FMAX,HMIN,HMAX,NPIX,
     &                     NUMEL,NBINS,BSIZ,HAV,HMODE,HSIG
         ELSE
            WRITE(NOUT,91) FMIN,FMAX,HMIN,HMAX,NPIX,
     &                     FNUMEL,FNBINS,BSIZ,HAV,HMODE,HSIG
         ENDIF
         WRITE(NOUT,*) ' '
      ENDIF

      RETURN
      END

C     -------------------- HISTIMAGE --------------------------------

      SUBROUTINE HISTIMAGE(LUN,LUNMA,NX,NREC,FREQ,
     &                     FNUMEL,HAV,HAV2,HMIN,FF,NBINS,MASK)

      INCLUDE 'CMBLOCK.INC'

      DIMENSION        FREQ(128), BUFMASK(NX), REDBUF(NX)
      DOUBLE PRECISION HAV, HAV2, FNUMEL
      LOGICAL          MASK

C     GET HISTOGRAM FROM IMAGE VALUES

      IF (MASK) THEN
C        MASKED, HANDLES NaN for NOT MASK OK
         DO  IREC=1,NREC
            CALL REDLIN(LUN,REDBUF,NX,IREC)
            CALL REDLIN(LUNMA,BUFMASK,NX,IREC)
            DO  ISAM = 1,NX
               IF (BUFMASK(ISAM) .GE. 0.5) THEN
C                 HISTOGRAM THIS PIXEL, MASK HAS POSITIVE VALUE)

C                 FIND BIN NUMBER
                  BVAL = REDBUF(ISAM)
                  JBIN = INT((BVAL - HMIN) * FF) + 1.5

                  IF (JBIN.GE.1 .AND. JBIN.LE.NBINS) THEN
C                    WITHIN HISTOGRAM RANGE
                     FREQ(JBIN) = FREQ(JBIN) + 1.0
                     FNUMEL     = FNUMEL + 1
                     HAV        = HAV  + BVAL 
                     HAV2       = HAV2 + DBLE(BVAL) * DBLE(BVAL)
                  ENDIF
              ENDIF
           ENDDO
         ENDDO
      ELSE
C        NO MASK
         DO  IREC=1,NREC
            CALL REDLIN(LUN,REDBUF,NX,IREC)

            DO  ISAM = 1,NX
C              HISTOGRAM THIS PIXEL 

C              FIND BIN NUMBER
               BVAL = REDBUF(ISAM)
               JBIN = INT((BVAL - HMIN) * FF) + 1.5

               IF (JBIN.GE.1 .AND. JBIN.LE.NBINS) THEN
C                 WITHIN HISTOGRAM RANGE
                  FREQ(JBIN) = FREQ(JBIN) + 1.0
                  FNUMEL     = FNUMEL + 1
                  HAV        = HAV  + BVAL 
                  HAV2       = HAV2 + DBLE(BVAL) * DBLE(BVAL)
               ENDIF
           ENDDO
         ENDDO
      ENDIF

      RETURN
      END

C     -------------------- HISTMINMAX --------------------------------

      SUBROUTINE HISTMINMAX(LUN,LUNMA,NX,NREC,
     &                      HMIN,HMAX)

      INCLUDE 'CMBLOCK.INC'


      DIMENSION        BUFMASK(NX), REDBUF(NX)

C     DETERMINE RANGE UNDER MASK
      DO IREC=1, NREC
         CALL REDLIN(LUN,  REDBUF,NX,IREC)
         CALL REDLIN(LUNMA,BUFMASK,NX,IREC)

         DO  ISAM = 1,NX
            IF (BUFMASK(ISAM) .GE. 0.5) THEN
C              PIXEL HAS POSITIVE MASK VALUE
               BVAL = REDBUF(ISAM)
               IF (BVAL .LT. HMIN) HMIN = BVAL
               IF (BVAL .GT. HMAX) HMAX = BVAL
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END

