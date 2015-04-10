C **********************************************************************
C                                                                      *
C  FBS2    NEW                            JUL  2011 Gregory Kishchenko *
C          CHKBOUND                       DEC  2011 ArDean Leith       *
C          NXP                            DEC  2011 ArDean Leith       *
C          /NXLD, WX DIM. BUG             DEC  2012 Gregory Kishchenko *
C **********************************************************************
C=* Author: Gregory Kishchenko                                         *
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
C                                                                      *
C   REAL FUNCTION FBS2(X,Y, NXLD,NX,NY, FDATA,NXP X1,Y1,XY2, CHKBOUND) *
C                                                                      *
C   PURPOSE:    2D FOURIER-BASED BICUBIC SPLINE INTERPOLATION          *
C               BETWEEN PIXELS.                                        *
C               ALGORITHM IS FAIRLY FAST AND PRESERVES FINE DETAILS    *
C               OF  IMAGES                                             *
C                                                                      *
C   PARAMETERS: X,Y        LOCATION OF PIXEL                     SENT  *
C               NXLD       FFTW PAD FOR X DIM.                   SENT  *
C               NX,NY      IMAGE DIMENSIONS                      SENT  *
C               FDATA      IMAGE                                 SENT  *
C               NXP        IMAGE X DIMENSION (PAD OR NOT)        SENT  *
C               X1,Y1,XY2  DERIVATIVES                           SENT  *
C               CHKBOUND   FLAG TO CHECK VALUES OF LOCATION      SENT  *
C                                                                      *
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
C       A3 = 2*(F(0)-F(1)) + F'(0) + F'(1)                             *
C                                                                      *
C     In order to interpolate a two dimensional grid [0,1] × [0,1], we *
C   sequentially used 1D cubic interpolation formula. First,  we       *
C   interpolated the intensities and normal to boundaries first        *
C   derivatives at two horizontal boundary lines [0,0]-[1,0] and       *
C   [0,1]-[1,1].                                                       *
C     For intensities' interpolation the intensities and tangential    *
C   first derivatives in grid nodes were used.                         *
C     For normal derivatives' interpolation the normal derivatives and *
C   cross-derivatives in grid nodes were used.                         *
C     Thereafter we carried out the vertical cubic interpolation on    *
C   line between 2 horizontal line with a given value of X to obtain   *
C   the intensity on vertical coordinate Y inside the square cell.     *
C   This last procedure was done using previously interpolated values  *
C   of intensities and their first normal derivatives across cell      *
C   boundaries.                                                        *
C       The first derivatives dF/dX, dF/dY and a cross-derivative      *
C   d2F/dXdY  in grid nodes were obtained by calculating 2D            *
C   Fourier transform of image, and then calculating the inverse       *
C   Fourier transform of {ik*F(k,l)}, il*F(k,l)}, and {ikl*F(k,l)}.    *
C      This well-known formula is computationally efficient and very   *
C   convenient for calculating the derivatives of a  function defined  *
C   as a discrete data set, and allows to calculate the derivative in  *
C   any local point without the finite difference approximation        *
C   involving the data from neighboring points.                        *
C                                                                      *
C      The result of this interpolation is similar to those of         *
C   standard bicubic spline interpolation (the densities and their     *
C   both derivatives are continuous at boundaries), but has fewer      *
C   interpolation artifacts, preserving the fine details and sharp     *
C   boundaries from blurring, because our algorithm involves the       *
C   densities and three partial derivatives just from 4 pixels         *
C   surrounding the area for interpolation instead 16 pixels as in     *
C   standard algorithm with the finite difference approximation.       *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      REAL FUNCTION FBS2(X,Y, NXLD,NX,NY, FDATA,NXP, X1,Y1,XY2,CHKBOUND)

      IMPLICIT NONE

      REAL    :: X, Y
      INTEGER :: NXLD,NX,NY
      REAL    :: FDATA(NXP,  NY)    ! IMAGE ARRAY
      INTEGER :: NXP                ! IMAGE X DIM.                
      REAL    :: X1   (NXLD, NY)    ! X  DERIVATIVE OF IMAGE
      REAL    :: Y1   (NXLD, NY)    ! Y  DERIVATIVE OF IMAGE
      REAL    :: XY2  (NXLD, NY)    ! XY DERIVATIVE OF IMAGE
      LOGICAL :: CHKBOUND           ! X,Y MAY BE OUT OF IMAGE?

      INTEGER :: I, J, I2, J2, I3, J3 
      REAL    :: A0, A1, A2, A3, B1
      REAL    :: ADX, BDX, DADX, DBDX
      REAL    :: DX, DY
      REAL    :: DXSQR,DXCUB

      I     = FLOOR(X)
      J     = FLOOR(Y)

      DX    = X - I
      DY    = Y - J
      DXSQR = DX **2
      DXCUB = DX **3

      IF (CHKBOUND) THEN
         I2 = MODULO(I-1,NX) + 1
         J2 = MODULO(J-1,NY) + 1
         I3 = MODULO(I  ,NX) + 1
         J3 = MODULO(J  ,NY) + 1
      ELSE
         I2 = I
         J2 = J
         I3 = I + 1
         J3 = J + 1
      ENDIF

      A0   = FDATA(I2,J2)
      A1   = X1(I2,J2)
      A2   = 3*( FDATA(I3,J2)-A0) -2*A1 - X1(I3,J2)
      A3   = 2*(-FDATA(I3,J2)+A0)  + A1 + X1(I3,J2)
      ADX  = A0 + A1*DX + A2*DXSQR + A3*DXCUB

      A0   = FDATA(I2,J3)
      A1   = X1(I2,J3)
      A2   = 3*( FDATA(I3,J3)-A0) - 2*A1 - X1(I3,J3)
      A3   = 2*(-FDATA(I3,J3)+A0) +   A1 + X1(I3,J3)
      BDX  = A0 + A1*DX + A2*DXSQR + A3*DXCUB

      A0   = Y1 (I2,J2)
      A1   = XY2(I2,J2)
      A2   = 3*( Y1(I3,J2)-A0) -2*A1 - XY2(I3,J2)
      A3   = 2*(-Y1(I3,J2)+A0) +  A1 + XY2(I3,J2)
      DADX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

      A0   = Y1 (I2,J3)
      A1   = XY2(I2,J3)
      A2   = 3*( Y1(I3,J3)-A0) - 2*A1 - XY2(I3,J3)
      A3   = 2*(-Y1(I3,J3)+A0) +   A1 + XY2(I3,J3)
      DBDX = A0 + A1*DX + A2*DXSQR + A3*DXCUB

      A2   = 3*(BDX - ADX) - 2*DADX - DBDX
      A3   = 2*(ADX - BDX) +   DADX + DBDX

      FBS2 = ADX + DADX * DY + A2 * DY**2 + A3 * DY**3

      END FUNCTION FBS2 


C *********************** FBS2_PREP ************************************
C                                                                      *
C   REAL FUNCTION FBS2(X,Y, NXLD,NX,NY, FDATA,NXP X1,Y1,XY2, CHKBOUND) *
C                                                                      *
C   PURPOSE:    2D FOURIER-BASED BICUBIC SPLINE INTERPOLATION          *
C               BETWEEN PIXELS.                                        *
C               ALGORITHM IS FAIRLY FAST AND PRESERVES FINE DETAILS    *
C               OF  IMAGES                                             *
C                                                                      *
C   PARAMETERS: F0         IMAGE                                 SENT  *
C               X1,Y1,XY2  DERIVATIVES                           RET   *
C               NXLD       FFTW PAD FOR X DIM.                   SENT  *
C               NX,NY      IMAGE DIMENSIONS                      SENT  *
C                                                                      *
C                                                                      *
C **********************************************************************

       SUBROUTINE FBS2_PREP(F0, X1,Y1,XY2, NXLD,NX,NY, IRTFLG)

       IMPLICIT NONE

       REAL              :: F0 (0:NXLD-1,0:NY-1)
       REAL              :: X1 (0:NXLD-1,0:NY-1) ! X  DERIVATIVE OF F0
       REAL              :: Y1 (0:NXLD-1,0:NY-1) ! Y  DERIVATIVE OF F0
       REAL              :: XY2(0:NXLD-1,0:NY-1) ! XY DERIVATIVE OF F0

       INTEGER           :: NXLD, NX, NY, IRTFLG

       REAL              :: WX(0:NXLD/2-1)          
       REAL              :: WY(0:NY-1)
            
       INTEGER           :: INV,I,J
       REAL              :: A4

c                           PI  = 3.14159265358979323846
       REAL, PARAMETER   :: PI2 = 6.28318530717958647692

C      FORWARD FFT
       INV = 1
       CALL FMRS(F0, NX,NY,1, 0.0D0, .TRUE.,.TRUE., INV, IRTFLG)

C      CALCULATE DERIVATIVES (dFDATA/dX and dFDATA/dY)

       A4 = PI2 / NX

       DO I=0,NXLD/2-1
          WX(I) = I * A4
       ENDDO

       A4 = PI2 / NY

       DO J=0,NY/2
          WY(J) = J * A4
       ENDDO

       DO J=NY/2+1,NY-1
          WY(J) = (J-NY)*A4
       ENDDO

       DO J=0,NY-1
          DO I=0,NXLD/2-1
             X1 (2*I,  J) =  F0(2*I+1,J)  * WX(I)
             X1 (2*I+1,J) = -F0(2*I,  J)  * WX(I)
             Y1 (2*I,  J) =  F0(2*I+1,J) * WY(J)
             Y1 (2*I+1,J) = -F0(2*I,  J) * WY(J)
             XY2(2*I,  J) =  X1(2*I+1,J) * WY(J)
             XY2(2*I+1,J) = -X1(2*I,  J) * WY(J)
          ENDDO
       ENDDO

C      REVERSE FFT
       INV = -1
       CALL FMRS(X1,  NX,NY,1, 0.0D0, .TRUE.,.TRUE., INV, IRTFLG)

C      REVERSE FFT
       INV = -1
       CALL FMRS(Y1,  NX,NY,1, 0.0D0, .TRUE.,.TRUE., INV, IRTFLG)

C      REVERSE FFT
       INV = -1
       CALL FMRS(XY2, NX,NY,1, 0.0D0, .TRUE.,.TRUE., INV, IRTFLG)

       END
