C ++********************************************************************
C                                                                      *
C  FBS3      NEW                          OCT 2011 GREGORY KISHCHENKO *
C            NXP ADDED                    DEC 2011 ARDEAN LEITH       *
C            /NXLD BUG, WX DIM. BUG       DEC 2012 Gregory Kishchenko *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NROW 12204.    *
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
C                                                                      *
C   FUNCTION FBS3(X,Y,Z,NXLD,NX,NY,NZ,F1,NXP,XYZ,X1,Y1,Z1,XY2,XZ2,YZ2) *
C                                                                      *
C   PURPOSE: 3D FOURIER-BASED TRICUBIC SPLINE INTERPOLATION            *
C            BETWEEN VOXELS. ALGORITHM IS FAIRLY FAST AND              *
C            PRESERVES FINE DETAILS OF IMAGES                          *
C                                                                      *
C   PARAMETERS: X,Y,Z         LOCATION WITHIN F1                  SENT *
C               NXLD          FFTW PAD SIZE                       SENT *
C               NX,NY,NZ      SIZE                                SENT *
C               F1            PADDED IMAGE                        SENT *
C               NXP           F1 X DIM                            SENT *
C               XYZ,X1,Y1     DERIVATIVES                         SENT *
C               XY2,XZ2,YZ2   DERIVATIVES                         SENT *
C                                                                      *
C   NOTES:                                                             *
C   If the values of a function F(X) and its first derivatives F'(X)   *
C   are known at X=0 and X=1, the function can be interpolated on the  *
C   interval [0,1] as a third degree polynomial (cubic                 *
C   interpolation formula):                                            *
C                                                                      *
C      F(X)=A0 + A1*X + A2*X**2 + A3*X**3                              *
C                                                                      *
C       where A0, A1, A2, and A3 are given by                          *
C       A0 = F(0)                                                      *
C       A1 = F'(0)                                                     *
C       A2 = 3*(F(1)-F(0) - 2*F'(0) - F'(1)                            *
C       A3 = 2*(F(0)-F(1)) +  F'(0) + F'(1)                            *
C                                                                      *
C   In order to interpolate a three-dimensional grid [0,1]×[0,1]×[0,1] *
C   we used a particular sequential combination of three               *
C   one-dimensional cubic interpolations.                              *
C       The following derivatives were used for interpolation:         *
C   dF/dX, dF/dY, dF/dZ, d2F/dXdY, d2F/dXdZ, d2F/dYdZ, d3F/dXdYdZ.     *
C      The use of this set satisfies the isotropic requirement of      *
C   interpolation                                                      *
C      All mentioned above derivatives for 8 grid nodes were obtained  *
C   by formula for derivatives' calculation using 3D Fourier transform.*
C   This approach is computationally efficient and very convenient     *
C   for calculating the derivatives of a  function defined as a        *
C   discrete data set, and allows to calculate the derivative in any   *
C   local point without the finite difference approximation involving  *
C   the data from neighboring points.                                  *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       REAL FUNCTION FBS3(X,Y,Z, 
     &               NXLD, NX,NY,NZ,
     &               F1,NXP, XYZ,X1,Y1,Z1,
     &               XY2,XZ2,YZ2)

       IMPLICIT NONE

       REAL          :: X,  Y,  Z
       INTEGER       :: NXLD, NX, NY, NZ

       REAL          :: F1 (NXP,  NY,NZ)  ! DATA CUBE 
       INTEGER       :: NXP

       REAL          :: XYZ(NXLD, NY,NZ)  ! XYZ DERIVATIVE OF F1  
       REAL          :: X1 (NXLD, NY,NZ)  ! X   DERIVATIVE OF F1
       REAL          :: Y1 (NXLD, NY,NZ)  ! Y   DERIVATIVE OF F1
       REAL          :: Z1 (NXLD, NY,NZ)  ! Z   DERIVATIVE OF F1
       REAL          :: XY2(NXLD, NY,NZ)  ! XY  DERIVATIVE OF F1
       REAL          :: XZ2(NXLD, NY,NZ)  ! XZ  DERIVATIVE OF F1
       REAL          :: YZ2(NXLD, NY,NZ)  ! YZ  DERIVATIVE OF F1

       INTEGER       :: I, J, L
       INTEGER       :: I2,J2,L2
       INTEGER       :: I3,J3,L3

       REAL          :: DX, DY, DZ
       REAL          :: A0, A1, A2, A3
       REAL          :: ADX, BDX, DADX, DBDX
       REAL          :: C0, DC0, C1, DC1
       REAL          :: DXSQR,DYSQR,DXCUB,DYCUB


             I  = FLOOR(X)
             J  = FLOOR(Y)
             L  = FLOOR(Z)

             I2 = MODULO(I-1,NX) + 1
             J2 = MODULO(J-1,NY) + 1
             L2 = MODULO(L-1,NZ) + 1

             I3 = I2 + 1
             IF (I3 > NX) I3 = 1
             J3 = J2 + 1
             IF (J3 > NY) J3 = 1
             L3 = L2 + 1
             IF (L3 > NZ) L3 = 1

             DX    = X - I
             DY    = Y - J
             DZ    = Z - L
             DXSQR = DX **2
             DXCUB = DX **3
             DYSQR = DY **2
             DYCUB = DY **3

C            1ST Plane - Function
             A0  = F1(I2,J2,L2)
             A1  = X1(I2,J2,L2)
             A2  = 3*( F1(I3,J2,L2)-A0) - 2*A1-X1(I3,J2,L2)
             A3  = 2*(-F1(I3,J2,L2)+A0) +   A1+X1(I3,J2,L2)
             ADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0 = F1(I2,J3,L2)
             A1 = X1(I2,J3,L2)
             A2 = 3*( F1(I3,J3,L2)-A0) - 2*A1 - X1(I3,J3,L2)
             A3 = 2*(-F1(I3,J3,L2)+A0) +   A1 + X1(I3,J3,L2)
             BDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0  = Y1(I2,J2,L2)
             A1  = XY2(I2,J2,L2)
             A2  = 3*( Y1(I3,J2,L2)-A0) -2*A1 - XY2(I3,J2,L2)
             A3  = 2*(-Y1(I3,J2,L2)+A0) +  A1 + XY2(I3,J2,L2)
             DADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0 = Y1(I2,J3,L2)
             A1 = XY2(I2,J3,L2)
             A2 = 3*( Y1(I3,J3,L2)-A0) - 2*A1 - XY2(I3,J3,L2)
             A3 = 2*(-Y1(I3,J3,L2)+A0) +   A1 + XY2(I3,J3,L2)
             DBDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A2 = 3*(BDX - ADX) - 2*DADX - DBDX
             A3 = 2*(ADX - BDX) +   DADX + DBDX

             C0 = ADX + DADX * DY + A2 * DYSQR + A3 * DYCUB

C            1ST Plane - Z-Derivative
             A0  = Z1(I2,J2,L2)
             A1  = XZ2(I2,J2,L2)
             A2  = 3*( Z1(I3,J2,L2)-A0) - 2*A1-XZ2(I3,J2,L2)
             A3  = 2*(-Z1(I3,J2,L2)+A0) +   A1+XZ2(I3,J2,L2)
             ADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0 = Z1(I2,J3,L2)
             A1 = XZ2(I2,J3,L2)
             A2 = 3*( Z1(I3,J3,L2)-A0) - 2*A1 - XZ2(I3,J3,L2)
             A3 = 2*(-Z1(I3,J3,L2)+A0) +   A1 + XZ2(I3,J3,L2)
             BDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0  = YZ2(I2,J2,L2)
             A1  = XYZ(I2,J2,L2)
             A2  = 3*( YZ2(I3,J2,L2)-A0) -2*A1 - XYZ(I3,J2,L2)
             A3  = 2*(-YZ2(I3,J2,L2)+A0) +  A1 + XYZ(I3,J2,L2)
             DADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0 = YZ2(I2,J3,L2)
             A1 = XYZ(I2,J3,L2)
             A2 = 3*( YZ2(I3,J3,L2)-A0) - 2*A1 - XYZ(I3,J3,L2)
             A3 = 2*(-YZ2(I3,J3,L2)+A0) +   A1 + XYZ(I3,J3,L2)
             DBDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A2 = 3*(BDX - ADX) - 2*DADX - DBDX
             A3 = 2*(ADX - BDX) +   DADX + DBDX

             DC0 = ADX + DADX*DY + A2*DYSQR + A3*DYCUB

C            2nd Plane - Function
             A0  = F1(I2,J2,L3)
             A1  = X1(I2,J2,L3)
             A2  = 3*( F1(I3,J2,L3)-A0) - 2*A1-X1(I3,J2,L3)
             A3  = 2*(-F1(I3,J2,L3)+A0) +  A1+X1(I3,J2,L3)
             ADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0 = F1(I2,J3,L3)
             A1 = X1(I2,J3,L3)
             A2 = 3*( F1(I3,J3,L3)-A0) - 2*A1 - X1(I3,J3,L3)
             A3 = 2*(-F1(I3,J3,L3)+A0) +   A1 + X1(I3,J3,L3)
             BDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0  = Y1 (I2,J2,L3)
             A1  = XY2(I2,J2,L3)
             A2  = 3*( Y1(I3,J2,L3)-A0) - 2*A1 - XY2(I3,J2,L3)
             A3  = 2*(-Y1(I3,J2,L3)+A0) +   A1 + XY2(I3,J2,L3)
             DADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0 = Y1 (I2,J3,L3)
             A1 = XY2(I2,J3,L3)
             A2 = 3*( Y1(I3,J3,L3)-A0) - 2*A1 - XY2(I3,J3,L3)
             A3 = 2*(-Y1(I3,J3,L3)+A0) +   A1 + XY2(I3,J3,L3)
             DBDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A2 = 3*(BDX - ADX) - 2*DADX - DBDX
             A3 = 2*(ADX - BDX) +   DADX + DBDX

             C1 = ADX + DADX* DY + A2 * DYSQR + A3 * DYCUB

C            2nd Plane - Z-Derivative
             A0  = Z1 (I2,J2,L3)
             A1  = XZ2(I2,J2,L3)
             A2  = 3*( Z1(I3,J2,L3)-A0) - 2*A1-XZ2(I3,J2,L3)
             A3  = 2*(-Z1(I3,J2,L3)+A0) +   A1+XZ2(I3,J2,L3)
             ADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0 = Z1 (I2,J3,L3)
             A1 = XZ2(I2,J3,L3)
             A2 = 3*( Z1(I3,J3,L3)-A0) - 2*A1 - XZ2(I3,J3,L3)
             A3 = 2*(-Z1(I3,J3,L3)+A0) +   A1 + XZ2(I3,J3,L3)
             BDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0  = YZ2(I2,J2,L3)
             A1  = XYZ(I2,J2,L3)
             A2  = 3*( YZ2(I3,J2,L3)-A0) -2*A1 - XYZ(I3,J2,L3)
             A3  = 2*(-YZ2(I3,J2,L3)+A0) +  A1 + XYZ(I3,J2,L3)
             DADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A0 = YZ2(I2,J3,L3)
             A1 = XYZ(I2,J3,L3)
             A2 = 3*( YZ2(I3,J3,L3)-A0) - 2*A1 - XYZ(I3,J3,L3)
             A3 = 2*(-YZ2(I3,J3,L3)+A0) +   A1 + XYZ(I3,J3,L3)
             DBDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

             A2 = 3*(BDX - ADX) - 2*DADX - DBDX
             A3 = 2*(ADX - BDX) +   DADX + DBDX

             DC1 = ADX + DADX * DY + A2 * DYSQR + A3 * DYCUB

C            FINAL LINE INSIDE CUBE
C            F(Z)=A0 + A1*DZ + A2*DZ**2 + A3*DZ**3
C              where A2, and A3 are given by
C            A2 = 3*(F(1)-F(0)) - 2*F'(0) - F'(1)
C            A3 = 2*(F(0)-F(1)) + F'(0) + F'(1)

             A2 = 3*(C1-C0) - 2*DC0 - DC1
             A3 = 2*(C0-C1) +   DC0 + DC1

             FBS3 = C0 + DC0*DZ + A2*DZ**2 + A3*DZ**3

       END FUNCTION FBS3 

C            TRI-LINEAR (JUST AN EXAMPLE)
C            C0 = F1(I2,J2,L2)*(1-DX)*(1-DY)*(1-DZ)
C            C0 = C0 + F1(I3,J2,L2)  * DX*(1-DY)*(1-DZ)
C            C0 = C0 + F1(I2,J3,L2)  * (1-DX)*DY*(1-DZ)
C            C0 = C0 + F1(I3,J3,L2)  * DX*DY*(1-DZ)
C            C0 = C0 + F1(I2,J2,L3)  * (1-DX)*(1-DY)*DZ
C            C0 = C0 + F1(I3,J2,L3)  * DX*(1-DY)*DZ
C            C0 = C0 + F1(I2,J3,L3)  * (1-DX)*DY*DZ
C            C0 = C0 + F1(I3,J3,L3)  * DX*DY*DZ
C            FBS3 = C0/8.0


C ************************ FBS3_PREP ***********************************
C                                                                      *
C   FUNCTION FBS3_PREP(F0, NXLD,NX,NY,NZ,                              *
C                      X1,Y1,Z1,XY2,XZ2,YZ2)                           *
C                                                                      *
C   PURPOSE: 3D FOURIER-BASED TRICUBIC SPLINE INTERPOLATION            *
C            BETWEEN VOXELS. ALGORITHM IS FAIRLY FAST AND              *
C            PRESERVES FINE DETAILS OF IMAGES                          *
C                                                                      *
C   PARAMETERS: F0            IMAGE                               SENT *
C               NXLD          FFTW PAD SIZE                       SENT *
C               NX,NY,NZ      SIZE                                SENT *
C               F1            PADDED IMAGE                        SENT *
C               NXP           F1 X DIM                            SENT *
C               XYZ,X1,Y1     DERIVATIVES                          RET *
C               XXY2,XZ2,YZ2  DERIVATIVES                          RET *
C                                                                      *
C **********************************************************************

       SUBROUTINE FBS3_PREP(F0, NXLD,NX,NY,NZ,
     &                      X1,Y1,Z1,XY2,XZ2,YZ2)

       REAL              :: F0 (NXLD, NY,NZ) ! XYZ DERIV. RETURNED
       REAL              :: X1 (NXLD, NY,NZ)
       REAL              :: Y1 (NXLD, NY,NZ)
       REAL              :: Z1 (NXLD, NY,NZ)
       REAL              :: XY2(NXLD, NY,NZ)
       REAL              :: XZ2(NXLD, NY,NZ)
       REAL              :: YZ2(NXLD, NY,NZ)
C      XYZ3 = F0
       INTEGER           :: NXLD, NX, NY, NZ
       INTEGER           :: K, J, L
       REAL              :: WX(0:NXLD/2-1)
       REAL              :: WY(0:NY-1)
       REAL              :: WZ(0:NZ-1)
       REAL              :: A4

       REAL, PARAMETER   :: PI2 = 6.28318530717958647692

C      FORWARD FFT OF: F0
       INV = +1
       CALL  FMRS_3(F0,NX,NY,NZ,INV)

C      CALCULATE DERIVATIVES (dF0/dX, dF0/dY, dF0/dZ)

       A4 = PI2 / NX

       DO K=0,NXLD/2-1
          WX(K) = K * A4
       ENDDO

       A4 = PI2 / NY

       DO J=0,NY/2
          WY(J) = J * A4
       ENDDO

       WY(NY/2) = 0

       DO J=NY/2+1,NY-1
          WY(J) = (J-NY) * A4
       ENDDO

       A4 = PI2 / NZ

       DO L=0,NZ/2
          WZ(L) = L * A4
       ENDDO

       WZ(NZ/2) = 0

       DO L=NZ/2+1,NZ-1
          WZ(L) = (L-NZ) * A4
       ENDDO

       ! ARRAY INITIALIZATIONS (needed??)
       X1  = 0
       Y1  = 0
       Z1  = 0
       XY2 = 0
       XZ2 = 0
       YZ2 = 0

c$omp  parallel do private(l,j,k)
       DO L=0,NZ-1
          DO J=0,NY-1
             DO K=0,NXLD/2-1
                X1 (2*K+1,J+1,L+1) =  F0(2*K+2,J+1,L+1) * WX(K)
                X1 (2*K+2,J+1,L+1) = -F0(2*K+1,J+1,L+1) * WX(K)
                Y1 (2*K+1,J+1,L+1) =  F0(2*K+2,J+1,L+1) * WY(J)
                Y1 (2*K+2,J+1,L+1) = -F0(2*K+1,J+1,L+1) * WY(J)
                XY2(2*K+1,J+1,L+1) =  X1(2*K+2,J+1,L+1) * WY(J)
                XY2(2*K+2,J+1,L+1) = -X1(2*K+1,J+1,L+1) * WY(J)
             ENDDO
          ENDDO
       ENDDO

c$omp  parallel do private(l,j,k)
       DO L=0,NZ-1
          DO J=0,NY-1
             DO K=0,NXLD/2-1
                Z1 (2*K+1,J+1,L+1)  =  F0(2*K+2,J+1,L+1)  * WZ(L)
                Z1 (2*K+2,J+1,L+1)  = -F0(2*K+1,J+1,L+1)  * WZ(L)

                XZ2(2*K+1,J+1,L+1)  =  X1(2*K+1,J+1,L+1)  * WZ(L)
                XZ2(2*K+2,J+1,L+1)  = -X1(2*K+1,J+1,L+1)  * WZ(L)
                YZ2(2*K+1,J+1,L+1)  =  Y1(2*K+2,J+1,L+1)  * WZ(L)
                YZ2(2*K+2,J+1,L+1)  = -Y1(2*K+1,J+1,L+1)  * WZ(L)

                F0 (2*K+1,J+1,L+1)  =  XY2(2*K+2,J+1,L+1) * WZ(L)
                F0 (2*K+2,J+1,L+1)  = -XY2(2*K+1,J+1,L+1) * WZ(L)

              ENDDO
           ENDDO
        ENDDO

C      REVERSE FFT
       INV= -1
       CALL  FMRS_3(X1, NX,NY,NZ,INV)
       INV= -1
       CALL  FMRS_3(Y1, NX,NY,NZ,INV)
       INV= -1
       CALL  FMRS_3(Z1, NX,NY,NZ,INV)
       INV= -1
       CALL  FMRS_3(XY2,NX,NY,NZ,INV)
       INV= -1
       CALL  FMRS_3(XZ2,NX,NY,NZ,INV)
       INV= -1
       CALL  FMRS_3(YZ2,NX,NY,NZ,INV)
       INV= -1
       CALL  FMRS_3(F0, NX,NY,NZ,INV)

       END
