C++*********************************************************************
C
C CENT_PH.F
C        DOUBLE PRECSISON INCREASES REPEATABILITY    JUN 05 ArDean Leith
C        RENAMED FROM CENT_D.F                       MAR 12 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C CENT_PH(BUF,NX,NY,SNS,SNR)
C
C PURPOSE: PHASE APPROXIMATION TO CENTER OF GRAVITY
C
C PARAMETERS:  BUF                 DATA ARRAY
C              NX,NY               SIZE
C              SNS,SNR             CENTER RELATIVE TO IMAGE ORIGIN
C
C--*********************************************************************

         SUBROUTINE CENT_PH(BUF, NX,NY, SNS,SNR)

         IMPLICIT NONE
         REAL, INTENT(IN)   :: BUF(NX,NY)
         INTEGER            :: NX,NY
         REAL               :: SNS,SNR

         INTEGER            :: I,J,IM1,JM1
         REAL               :: RT
         INTEGER, PARAMETER :: R_8 = SELECTED_REAL_KIND(P=8)
         REAL(KIND=R_8)     :: C,S,P,T,FI

         C = 0.0
         S = 0.0
         P = 8 * DATAN(1.0D0) / NX

c$omp    parallel do private(i,im1,j,t),reduction(+:c,s)
         DO I=1,NX
            IM1 = (I-1)
            T = 0.0
            DO J=1,NY
               T = T + BUF(I,J)
            ENDDO

            C   = C + COS(P * IM1) * T
            S   = S + SIN(P * IM1) * T
         ENDDO

         FI = ATAN2(S,C)
         IF (FI < 0.0)  FI = FI + 8 * DATAN(1.0D0)
         SNS = FI / P  +1.0

         C = 0.0
         S = 0.0
         P = 8 * DATAN(1.0D0) / NY

c$omp    parallel do private(j,jm1,i,t),reduction(+:c,s)
         DO J=1,NY
            JM1 = (J-1)

            T = 0.0
            DO I=1,NX
               T = T + BUF(I,J)
            ENDDO

            C   = C + COS(P * JM1) * T
            S   = S + SIN(P * JM1) * T
         ENDDO

         FI = ATAN2(S,C)
         IF (FI < 0.0) FI = FI + 8 * DATAN(1.0D0)
         SNR = FI / P + 1.0

         END
