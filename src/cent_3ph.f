C++*********************************************************************
C
C CENT_3PH.F
C         REAL(KIND=R_8)     :: C,S,P,T,FI          July 2005 A. Leith
C         RENAMED FROM CENT_3                       Mar  2012 A. Leith
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
C  CENT_3PH(LUN,NX,NY,NZ,SNS,SNR,SNC)
C
C  PURPOSE: PHASE APPROXIMATION TO CENTER OF GRAVITY
C
C--*********************************************************************

         SUBROUTINE  CENT_3PH(LUN,NX,NY,NZ,SNS,SNR,SNC)

         REAL               ::  B(NX),X(NX),Y(NY),Z(NZ)

         INTEGER, PARAMETER :: R_8 = SELECTED_REAL_KIND(P=8)
         REAL(KIND=R_8)     :: C,S,P,T,FI

         X = 0.0
         Y = 0.0
         Z = 0.0
        
         DO K=1,NZ
            DO J=1,NY
               CALL REDLIN(LUN,B,NX,J+(K-1)*NY)
               T    = SUM(B)
               X    = X + B
               Y(J) = Y(J) + T
               Z(K) = Z(K) + T
            ENDDO
         ENDDO

         C = 0.0
         S = 0.0
         P = 8 * DATAN(1.0D0) / NX

         DO I=1,NX
            C = C + COS(P * (I-1)) * X(I)
            S = S + SIN(P * (I-1)) * X(I)
         ENDDO

         FI = ATAN2(S,C)
         IF (FI .LT. 0.0)  FI = FI + 8 * DATAN(1.0D0)
         SNS = FI / P + 1.0
C
         C = 0.0
         S = 0.0
         P = 8 * DATAN(1.0D0) / NY

         DO    J=1,NY
            C = C + COS(P * (J-1)) * Y(J)
            S = S + SIN(P * (J-1)) * Y(J)
         ENDDO

         FI = ATAN2(S,C)
         IF (FI < 0.0)  FI = FI + 8 * DATAN(1.0D0)
         SNR=FI/P+1.0

         C = 0.0
         S = 0.0
         P = 8 * DATAN(1.0D0) / NZ

         DO K=1,NZ
            C = C + COS(P * (K-1)) * Z(K)
            S = S + SIN(P * (K-1)) * Z(K)
         ENDDO

         FI = ATAN2(S,C)
         IF (FI < 0.0)  FI = FI + 8 * DATAN(1.0D0)

         SNC = FI / P + 1.0

         END
