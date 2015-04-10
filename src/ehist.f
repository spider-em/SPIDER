
C++*********************************************************************
C
C    EHIST.F             FILENAMES LENGTHENED JAN 89 al
C                        REMOVED RHIST AUG 96 al
C                        VOLUMES APRIL 2002 al
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
C    EHIST(FILNAM,LUN1,LUN2,NSAM,NROW)
C
C    PURPOSE:      HISTOGRAM EQUALIZATION
C
C    PARAMETERS:
C        FILNAM    NAME OF FILE
C        LUN1      LOGICAL UNIT NUMBER OF FILE
C        LUN2      LOGICAL UNIT NUMBER OF FILE
C        NSAM      NUMBER OF SAMPLES
C        NROW      NUMBER OF ROWS
C        NSLICE    NUMBER OF SLICES
C
C--*********************************************************************

      SUBROUTINE EHIST(FILNAM,LUN1,LUN2,NSAM,NROW,NSLICE)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      COMMON         BUF(1)
      COMMON /IOBUF/ REDBUF(NBUFSIZ)

      CHARACTER *(*) FILNAM
      CHARACTER      NULL,ANS
      LOGICAL        MAPOUT,HISOUT

      NULL   =  CHAR(0)

      MAPOUT = .FALSE.
      HISOUT = .FALSE.
      CALL RDPRMC(ANS,NC,.TRUE.,
     &   'MAPPING FUNCTION AND HISTOGRAM PRINTOUT? (Y/N)',NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (ANS .NE. 'N') THEN
        MAPOUT = .TRUE.
        HISOUT = .TRUE.
      ENDIF

C     FIND HISTOGRAM OF IMAGE, PLACE IT IN BUF
CCC   CALL HIST(LUN1,NSAM,NSLICE,NROW,0,0,.FALSE.) aug 98 al
      CALL HIST(LUN1,0,0,NSAM,NROW,NSLICE,HMIN,HMAX,HSIG,HMODE)

C     FIND CUMULATIVE DISTRIBUTION (HISTOGRAM MAPPING FUNCTION)
      DO K = 2,128
         BUF(K) = BUF(K) + BUF(K-1)
      ENDDO

      IF (MAPOUT) THEN
C       PRINT HISTOGRAM MAPPING FUNCTION IN RESULTS FILE
        CALL PDATES(FILNAM,0)
        WRITE(NDAT,*) ' HISTOGRAM MAPPING FUNCTION'
        CALL GRAPHS(NDAT,BUF,128,1,1,1.0,IRTFLG)
CCC     CALL GRAPHS(BUF(NSAM+1),128,1) aug 98 al
      ENDIF

C     HMIN IS MIN OF CUMULATIVE FUNCTION, HMAX IS MAX (TOTAL # OF PIXELS)
      HMIN = BUF(1)
      HMAX = BUF(1)
      DO  K = 1,128
         HB = BUF(K)
         IF (HB .LT. HMIN) HMIN = HB
         IF (HB .GT. HMAX) HMAX = HB
      ENDDO

      VAL = 2.0 / (HMAX-HMIN)
      DO  K = 1,128
        BUF(K) = (BUF(K) - HMIN) * VAL
      ENDDO

C     HINC IS INCREMENT OF ORIGINAL IMAGE RANGE CORRESPOINDING TO 1 HISTOGRAM UNIT
      HINC  = (FMAX-FMIN) / 127.0
      FAC   = (-FMIN/HINC)+ 1.5
      HINCR = 1.0 / HINC

C     CONVERT IMAGE VALUES TO MAPPED VALUES
      DO  I = 1,NROW*NSLICE
         CALL REDLIN(LUN1,REDBUF,NSAM,I)
         DO  K = 1,NSAM

C          MAP IMAGE VALUE TO RANGE 1...128
C**        MAP  = (REDBUF(K)-FMIN)/HINC +1.5
           MAP  = REDBUF(K) * HINCR + FAC

C          REPLACE IMAGE VALUE BY HISTOGRAM MAPPED VALUE
           REDBUF(K) = BUF(MAP)
	 ENDDO
         CALL WRTLIN(LUN2,REDBUF,NSAM,I)
      ENDDO

C     SET IMAGE HEADER FOR NEW IMAGE CONTENTS
      CALL SETPRM(LUN2,NSAM,NROW,0.0,0.0,0.0,'U')

C     PRINT HISTOGRAM IN RESULTS FILE IF DESIRED
      IF (HISOUT) CALL HIST(LUN2,0,0,NSAM,NROW,NSLICE,
     &                      HMIN,HMAX,HSIG,HMODE)

      RETURN
      END
