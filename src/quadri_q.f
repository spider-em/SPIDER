C ++********************************************************************
C                                                                      *
C QUADRI_Q                                                                     *
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
C                                                                      *
C  FUNCTION QUADRI_Q(XX, YY, LSD, NXDATA, NYDATA, FDATA)               *
C                                                                      *
C  PURPOSE: QUADRATIC INTERPOLATION                                                           *
C                                                                      *
C  PARAMETERS:       XX,YY TREATED AS CIRCULARLY CLOSED.
C                    FDATA - IMAGE 1..NXDATA, 1..NYDATA
C
C                    F3    FC       F0, F1, F2, F3 are the values
C                     +             at the grid points.  X is the
C                     + X           point at which the function
C              F2++++F0++++F1       is to be estimated. (It need
C                     +             not be in the First quadrant).
C                     +             FC - the outer corner point
C                    F4             nearest X.
C
C                                   F0 is the value of the FDATA at
C                                   FDATA(I,J), it is the interior mesh
C                                   point nearest  X.
C                                   The coordinates of F0 are (X0,Y0),
C                                   The coordinates of F1 are (XB,Y0),
C                                   The coordinates of F2 are (XA,Y0),
C                                   The coordinates of F3 are (X0,YB),
C                                   The coordinates of F4 are (X0,YA),
C                                   The coordinates of FC are (XC,YC),
C
C                   O               HXA, HXB are the mesh spacings
C                   +               in the X-direction to the left
C                  HYB              and right of the center point.
C                   +
C            ++HXA++O++HXB++O       HYB, HYA are the mesh spacings
C                   +               in the Y-direction.
C                  HYA
C                   +               HXC equals either  HXB  or  HXA
C                   O               depending on where the corner
C                                   point is located.
c
C                                   Construct the interpolant
C                                   F = F0 + C1*(X-X0) +
C                                       C2*(X-X0)*(X-X1) +
C                                       C3*(Y-Y0) + C4*(Y-Y0)*(Y-Y1)
C                                       + C5*(X-X0)*(Y-Y0)
C
C  IMAGE_PROCESSING_ROUTINE 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      FUNCTION QUADRI_Q(XX, YY, LSD, NXDATA, NYDATA, FDATA)

      DIMENSION  FDATA(LSD,NYDATA)

      X = XX
      Y = YY

C     CIRCULAR CLOSURE
      IF (X.LT.1.0)               X = X+(1 - IFIX(X) / NXDATA) * NXDATA
      IF (X.GT.FLOAT(NXDATA)+0.5) X = AMOD(X-1.0,FLOAT(NXDATA)) + 1.0
      IF (Y.LT.1.0)               Y = Y+(1 - IFIX(Y) / NYDATA) * NYDATA
      IF (Y.GT.FLOAT(NYDATA)+0.5) Y = AMOD(Y-1.0,FLOAT(NYDATA)) + 1.0

      I   = IFIX(X)
      J   = IFIX(Y)

      DX0 = X - I
      DY0 = Y - J

      IP1 = I + 1
      IM1 = I - 1
      JP1 = J + 1
      JM1 = J - 1

      IF (IP1 .GT. NXDATA) IP1 = IP1 - NXDATA     
      IF (IM1 .LT. 1)      IM1 = IM1 + NXDATA
      IF (JP1 .GT. NYDATA) JP1 = JP1 - NYDATA        
      IF (JM1 .LT. 1)      JM1 = JM1 + NYDATA

      F0  = FDATA(I,J)
      C1  = FDATA(IP1,J) - F0
      C2  = (C1 - F0 + FDATA(IM1,J)) * 0.5
      C3  = FDATA(I,JP1) - F0 
      C4  = (C3 - F0 + FDATA(I,JM1)) * 0.5 

      DXB = (DX0 - 1)
      DYB = (DY0 - 1)

C     HXC & HYC ARE EITHER 1 OR -1
      HXC = INT(SIGN(1.0,DX0))
      HYC = INT(SIGN(1.0,DY0)) 
 
      IC  = I + HXC
      JC  = J + HYC

      IF (IC .GT .NXDATA) THEN
         IC = IC - NXDATA    
      ELSEIF (IC .LT. 1)  THEN
         IC = IC + NXDATA
      ENDIF

      IF (JC .GT. NYDATA)  THEN
         JC = JC - NYDATA
      ELSEIF (JC .LT. 1)  THEN
         JC = JC + NYDATA
      ENDIF

      C5  =  ((FDATA(IC,JC) - F0 - 
     &         HXC * C1 - 
     &        (HXC * (HXC - 1.0)) * C2 -
     &         HYC * C3 - 
     &        (HYC * (HYC - 1.0)) * C4) * 
     &        (HXC * HYC)) 

      QUADRI_Q = F0 + 
     &         DX0 * (C1 + DXB * C2 + DY0 * C5) + 
     &         DY0 * (C3 + DYB * C4)


      END
