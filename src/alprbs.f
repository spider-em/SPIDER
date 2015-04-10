C ++********************************************************************
C                                                                      *
C   ALPRBS.F                                                           *
C                                                                      *
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
C  ALPRBS(NUMR,NRING,LCIRC,MODE)                                                                     *
C                                                                       
C  PURPOSE: APPEARS TO EXTRACT CIRCULAR RINGS AND POSTITION
C           IN A LINEAR ARRAY THAT HOLDS RINGS CONCATENATED TOGETHER.
C           OUTPUT IS DEPENDENT ON NUMBER OF RINGS 
C           FOR USE WITH SPIDER (e.g. NOT!! FFTW3) FFT.
C                                                                       
C  PARAMETERS: NUMR(1,I) - RING NUMBER                      (SENT)
C              NUMR(2,I) - BEGINNING IN CIRC                (RET.)
C              NUMR(3,I) - LENGTH IN CIRC                   (RET.)
C              NRING                                        (SENT)
C              LCIRC - TOTAL LENGTH OF CIRC.                (RET.)
C
C***********************************************************************

        SUBROUTINE  ALPRBS(NUMR,NRING,LCIRC,MODE)

        INCLUDE 'CMLIMIT.INC'

        INTEGER      NUMR(3,NRING)
        REAL*8       PI
        CHARACTER*1  MODE

C       PREPARATION OF PARAMETERS

        PI = 4.0*DATAN(1.0D0)
        IF (MODE.EQ.'F')  PI=2*PI

        LCIRC = 0

        DO I=1,NRING
           JP = PI*NUMR(1,I)
           IP = 2**LOG2(JP)

C          IF CLOSE TO BOUNDARY OF POWER OF TWO, GOTO NEXT POWER
           IF (I.LT.NRING .AND. JP.GT.IP+IP/2)  IP=MIN0(MAXFFT,2*IP)

C          LAST RING OVERSAMPLED TO ALLOW HIGHER ACCURACY
C          OF PEAK LOCATION (?).
           IF (I.EQ.NRING .AND. JP.GT.IP+IP/5)  IP=MIN0(MAXFFT,2*IP)

C          ALL THE RINGS ARE POWER-OF-TWO.
           NUMR(3,I) = IP

           IF (I .EQ. 1)  THEN
              NUMR(2,1) = 1
           ELSE
              NUMR(2,I) = NUMR(2,I-1)+NUMR(3,I-1)
           ENDIF

           LCIRC = LCIRC + IP
        ENDDO 

        END
