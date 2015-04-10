C++*********************************************************************
C
C ALPRBS_Q.F
C                 FOR FFTW3                      MAR. 2008 ARDEAN LEITH
C                 ALPRBS_Q_NEW                   JUN. 2008 ARDEAN LEITH
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
C ALPRBS_Q(NUMR,NRING,LCIRC,MODE)
C
C PURPOSE: CREATES NUMR ARRAY HOLDING THE SPECIFICATIONS FOR THE 
C          RADIAL RINGS, FOR USE WITH NON-SPIDER (e.g. FFTW3) FFT.
C
C PARAMETERS:   NRING       NUMBER OF RINGS                     SENT
C               NUMR(1,I)   RING NUMBER = RADIUS                SENT
C               NUMR(2,I)   BEGINNING IN CIRC                    RET.
C               NUMR(3,I)   LENGTH IN CIRC                       RET.
C               LCIRC       TOTAL LENGTH OF CIRC                 RET.
C               IRAY        RAY SKIP INCREMENT (1,2,4,8...)     SENT.
C   
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ALPRBS_Q(NUMR,NRING,LCIRC,MODE)

        INTEGER, INTENT(INOUT)   :: NUMR(3,NRING)
        CHARACTER*1, INTENT(IN)  :: MODE

        REAL*8                   :: PI
        INTEGER                  :: MAXFFTT

#ifndef SP_LIBFFTW3
C       FOR MAXFFT
        INCLUDE 'CMLIMIT.INC'

        MAXFFTT = MAXFFT
#else
        MAXFFTT = HUGE(MAXFFTT)
#endif

C       PREPARATION OF PARAMETERS
        PI = 4.0 * DATAN(1.0D0)
        IF (MODE .EQ. 'F')  PI = 2 * PI

        LCIRC = 0              ! TOTAL LENGTH OF CIRCLE ARRAYS

        DO I=1,NRING
           JP = PI * NUMR(1,I)     ! CIRCUMFERENCE OF CIRCLE
           IP = 2 ** LOG2(JP)      ! LENGTH IS POWER OF 2: 8,16,32,64....

C          IF CLOSE TO BOUNDARY OF POWER OF TWO, GOTO NEXT POWER
           IF (I .LT. NRING .AND. JP .GT. IP+IP/2)IP = MIN(MAXFFTT,2*IP)

C          LAST RING OVERSAMPLED TO IMPROVE ACCURACY OF PEAK LOCATION (?).
           IF (I.EQ. NRING .AND. JP .GT. IP+IP/5) IP = MIN(MAXFFTT,2*IP)

C          ALL THE RINGS ARE POWER-OF-TWO. ADD 2 TO LEN TO USE NEW FFT.
           IP        = IP + 2
           NUMR(3,I) = IP                ! LENGTH OF CIRCLE + PAD
           LCIRC     = LCIRC + IP        ! TOTAL LENGTH OF CIRCLE ARRAYS

           IF (I .EQ. 1) THEN
              NUMR(2,1) = 1
           ELSE
              NUMR(2,I) = NUMR(2,I-1) + NUMR(3,I-1)
           ENDIF

        ENDDO

        END

C       ----------- ALPRBS_Q_NEW ------------------------------------

        SUBROUTINE ALPRBS_Q_NEW(NUMR,NRING,LCIRC,MODE,IRAY)

        INTEGER, INTENT(INOUT)   :: NUMR(3,NRING)
        CHARACTER*1, INTENT(IN)  :: MODE

        REAL*8                   :: PI
        INTEGER                  :: MAXFFTT

#ifndef SP_LIBFFTW3
C       FOR MAXFFT
        INCLUDE 'CMLIMIT.INC'

        MAXFFTT = MAXFFT
#else
        MAXFFTT = HUGE(MAXFFTT)
#endif

C       PREPARATION OF PARAMETERS
        PI = 4.0 * DATAN(1.0D0)
        IF (MODE .EQ. 'F')  PI = 2 * PI

        LCIRC = 0              ! TOTAL LENGTH OF CIRCLE ARRAYS

        DO I=1,NRING
           JP = PI * NUMR(1,I)
           IP = 2 ** LOG2(JP)
           IF (I .LT. NRING .AND. JP .GT. IP+IP/2)IP = MIN(MAXFFTT,2*IP)
      
C          LAST RING OVERSAMPLED TO IMPROVE ACCURACY OF PEAK LOCATION (?).
           IF (I.EQ. NRING .AND. JP .GT. IP+IP/5) IP = MIN(MAXFFTT,2*IP)

           IF (IRAY .GT. 1) THEN
               !IPT = IP / IRAY 
               !write(6,*) ' ring,jp,ip --> ip:',i,jp,ip,ipt
               IP = IP / IRAY 
           ENDIF

C          ALL THE RINGS ARE POWER-OF-TWO. ADD 2 TO LEN TO USE NEW FFT.
           IP        = IP + 2
           NUMR(3,I) = IP                ! LENGTH OF CIRCLE + PAD
           LCIRC     = LCIRC + IP        ! TOTAL LENGTH OF CIRCLE ARRAYS

           IF (I .EQ. 1) THEN
              NUMR(2,1) = 1
           ELSE
              NUMR(2,I) = NUMR(2,I-1) + NUMR(3,I-1)
           ENDIF

        ENDDO

        END
