C ++********************************************************************
C                                                                      *
C QUADRI                                                               *
C                  QUADRI_FAST WITHOUT TRAPS     JUN 2008 ARDEAN LEITH *
C                  QUADRI_PAD                    OCT 2012 ARDEAN LEITH *
C                                                                      *
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
C  FUNCTION QUADRI(XX, YY, NX, NY, FDATA)
C                                                                      
C  PURPOSE: QUADRATIC INTERPOLATION. QUADRI HAS CIRCULAR CLOSING STEPS 
C           WHICH   QUADRI_FAST LACKS                                  
C                                                                      
C  PARAMETERS:       XX,YY TREATED AS CIRCULARLY CLOSED.
C                    FDATA - IMAGE 1..NX, 1..NYC
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
C                                   Coordinates of F0 are (X0,Y0),
C                                   Coordinates of F1 are (XB,Y0),
C                                   Coordinates of F2 are (XA,Y0),
C                                   Coordinates of F3 are (X0,YB),
C                                   Coordinates of F4 are (X0,YA),
C                                   Coordinates of FC are (XC,YC),
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
C                                       C3*(Y-Y0) + C4*(Y-Y0)*(Y-Y1) +
C                                       C5*(X-X0)*(Y-Y0)
C
C NOTE: QUADRI-FAST CAN BE USED WHEN THERE IS NO CHANCE THAT THE 
C       INPUT DATA GOES OUTSIDE THE BOUNDS OF THE IMAGE.  THE MAIN
C       CASE THAT THIS IS TRUE IS WHEN CREATING POLAR IMAGES AND
C       YOU CAN BE SURE THAT THE RADIUS IS < BOUNDARY.
C
C       MAYBE QUADRI SHOULD BE INLINED?
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

C     ------------------- QUADRI -----------------------------------

      FUNCTION QUADRI(XX,YY, NX,NY, FDATA)

      DIMENSION  FDATA(NX,NY)

      X = XX
      Y = YY

C     CIRCULAR CLOSURE
      IF (X < 1.0)               X = X+(1 - IFIX(X) / NX) * NX
      IF (X > FLOAT(NX)+0.5) X = AMOD(X-1.0,FLOAT(NX)) + 1.0
      IF (Y < 1.0)               Y = Y+(1 - IFIX(Y) / NY) * NY
      IF (Y > FLOAT(NY)+0.5) Y = AMOD(Y-1.0,FLOAT(NY)) + 1.0

      I   = IFIX(X)
      J   = IFIX(Y)

      DX0 = X - I
      DY0 = Y - J

      IP1 = I + 1
      IM1 = I - 1
      JP1 = J + 1
      JM1 = J - 1

      IF (IP1 > NX) IP1 = IP1 - NX     
      IF (IM1 < 1)      IM1 = IM1 + NX
      IF (JP1 > NY) JP1 = JP1 - NY        
      IF (JM1 < 1)      JM1 = JM1 + NY
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

      IF (IC .GT .NX)   THEN
         IC = IC - NX    
      ELSEIF (IC  <  1) THEN
         IC = IC + NX
      ENDIF

      IF (JC  >  NY)    THEN
         JC = JC - NY      
      ELSEIF (JC  <  1) THEN
         JC = JC + NY      
      ENDIF

      C5  =  ((FDATA(IC,JC) - F0 - 
     &         HXC * C1 - 
     &        (HXC * (HXC - 1.0)) * C2 -
     &         HYC * C3 - 
     &        (HYC * (HYC - 1.0)) * C4) * 
     &        (HXC * HYC)) 

      QUADRI = F0 + 
     &         DX0 * (C1 + DXB * C2 + DY0 * C5) + 
     &         DY0 * (C3 + DYB * C4)


      END

C     ------------------- QUADRI_FAST --------------------------------

      FUNCTION QUADRI_FAST(X,Y, NX,NY, FDATA)

      REAL    :: X,Y
      REAL    :: FDATA(NX,NY)
      INTEGER :: NX, NY

C     SKIP CIRCULAR CLOSURE, IT IS SLOW, ENSURE IT NOT NEEDED IN CALLER

      I   = IFIX(X)
      J   = IFIX(Y)

      DX0 = X - I
      DY0 = Y - J

      IP1 = I + 1
      IM1 = I - 1
      JP1 = J + 1
      JM1 = J - 1
 
      F0  = FDATA(I,J)
      C1  = FDATA(IP1,J) - F0               ! DIFF. FROM CENTER
      C2  = (C1 - F0 + FDATA(IM1,J)) * 0.5  ! DIFF OF X+1 AND X-1
      C3  = FDATA(I,JP1) - F0               ! DIFF. FROM CENTER
      C4  = (C3 - F0 + FDATA(I,JM1)) * 0.5  ! DIFF oF Y+1 AND Y-1

      DXB = (DX0 - 1)
      DYB = (DY0 - 1)

C     HXC & HYC ARE EITHER 1 OR -1
      HXC = INT(SIGN(1.0,DX0))   ! X <> INT(X)
      HYC = INT(SIGN(1.0,DY0))   ! Y <> INT(Y)
 
      IC  = I + HXC
      JC  = J + HYC

      C5  =  ((FDATA(IC,JC) - F0 - 
     &         HXC * C1 - 
     &        (HXC * (HXC - 1.0)) * C2 -
     &         HYC * C3 - 
     &        (HYC * (HYC - 1.0)) * C4) * 
     &        (HXC * HYC)) 

      QUADRI_FAST = F0 + 
     &         DX0 * (C1 + DXB * C2 + DY0 * C5) + 
     &         DY0 * (C3 + DYB * C4)

      END


C     ------------------- QUADRI_PAD -------------------------------

      REAL FUNCTION QUADRI_PAD(XX, YY, NX,NY, NXP,NYP, FDATA)

      IMPLICIT NONE

      REAL    :: XX, YY
      REAL    :: FDATA(NXP,NYP)
      INTEGER :: NX,NY, NXP,NYP

      REAL    :: X, Y, DX0,DY0,F0,C1,C2,C3,DXB,DYB,c4,hxc,hyc,c5
      INTEGER :: I,J,IP1,IM1,JP1,JM1,ic,jc

      X = XX
      Y = YY

C     CIRCULAR CLOSURE
      IF (X < 1.0)           X = X+(1 - IFIX(X) / NX) * NX
      IF (X > FLOAT(NX)+0.5) X = AMOD(X-1.0,FLOAT(NX)) + 1.0
      IF (Y < 1.0)           Y = Y+(1 - IFIX(Y) / NY) * NY
      IF (Y > FLOAT(NY)+0.5) Y = AMOD(Y-1.0,FLOAT(NY)) + 1.0

      I   = IFIX(X)
      J   = IFIX(Y)

      DX0 = X - I
      DY0 = Y - J

      IP1 = I + 1
      IM1 = I - 1
      JP1 = J + 1
      JM1 = J - 1

      IF (IP1 > NX) IP1 = IP1 - NX     
      IF (IM1 < 1)  IM1 = IM1 + NX
      IF (JP1 > NY) JP1 = JP1 - NY        
      IF (JM1 < 1)  JM1 = JM1 + NY

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

      IF (IC > NX)      THEN
         IC = IC - NX    
      ELSEIF (IC  <  1) THEN
         IC = IC + NX
      ENDIF

      IF (JC  >  NY)    THEN
         JC = JC - NY      
      ELSEIF (JC  <  1) THEN
         JC = JC + NY 
      ENDIF

      C5  =  ((FDATA(IC,JC) - F0 - 
     &         HXC * C1 - 
     &        (HXC * (HXC - 1.0)) * C2 -
     &         HYC * C3 - 
     &        (HYC * (HYC - 1.0)) * C4) * 
     &        (HXC * HYC)) 

      QUADRI_PAD = F0 + 
     &         DX0 * (C1 + DXB * C2 + DY0 * C5) + 
     &         DY0 * (C3 + DYB * C4)


      END





#ifdef NEVER
C     ------------------- FLIN -----------------------------------
     FUNCTION FLIN(X, Y, NX, NY, FDATA)

      REAL    :: X,Y
      REAL    :: FDATA(NX,NY)
      INTEGER :: NX, NY

      I   = IFIX(X)
      J   = IFIX(Y)

      IP1 = I + 1
      JP1 = J + 1
 
      F0  = FDATA(I,J)
      F1  = FDATA(IP1,J)
      F3  = FDATA(I,JP1)

      C1  = FDATA(IP1,J) - F0               ! DIFF. FROM CENTER
      C2  = (C1 - F0 + FDATA(IM1,J)) * 0.5  ! DIFF of x+1 and x-1
      C3  = FDATA(I,JP1) - F0               ! DIFF. FROM CENTER
      C4  = (C3 - F0 + FDATA(I,JM1)) * 0.5  ! DIFF of y+1 and y-1


      END

C          GET VALUE OF Y1, BILINEAR INTERPRELATAION
           Y1(I) = (1.-(K(I)-INT(K(I)))) * BUF(INT(K(I))) +
     &                 (K(I)-INT(K(I)))  * BUF(INT(K(I))+1)

        RICENT = ICENT + SHY
        RKCENT = KCENT + SHX  

        JJ     = 0
        DO  I = 1,NROW
           JJ = JJ+1
           IF (IFLAG1 .EQ. 1) II = NROWS + (I-1) * NROWSK
           IF (IFLAG1 .EQ. 0) II = I
           Y    = I - RICENT
           YCOD = Y * COD + RICENT
           YSID = -Y * SID + RKCENT
           DO K = 1,NSAM
              RBUF(K) = BACK
              X       = K - RKCENT
              XOLD    = X * COD + YSID
              YOLD    = X * SID + YCOD
              IYOLD   = YOLD
              YDIF    = YOLD - IYOLD
              YREM    = 1.   - YDIF
              IXOLD   = XOLD
              IF ((IYOLD .GE. 1 .AND. IYOLD .LE. NROW-1) .AND.
     &            (IXOLD .GE. 1 .AND. IXOLD .LE. NSAM-1)) THEN
c                INSIDE BOUNDARIES OF OUTPUT IMAGE
                 XDIF    = XOLD - IXOLD
                 XREM    = 1. - XDIF
                 NADDR   = (IYOLD-1) * NSAM + IXOLD 
                 RBUF(K) = YDIF*(BUF(NADDR+NSAM)*XREM
     &                    +BUF(NADDR+NSAM+1)*XDIF)
     &                    +YREM*(BUF(NADDR)*XREM + BUF(NADDR+1)*XDIF)
              ENDIF
           ENDDO

           CALL WRTLIN(LUNO,RBUF,NSAM,II)
        ENDDO
#endif
