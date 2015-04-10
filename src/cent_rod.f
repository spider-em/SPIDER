C **********************************************************************
C                                                                      *
C CENT_ROD       NEW                      FEB 2012  GREGORY KISHCHENKO *
C                ROTATION BUG FIXED       AUG 2014  ARDEAN LEITH       *                                                               *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NROW 12204.    *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) aNROW later version.                  *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANROW WARRANTY; without even the implied warranty of   *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C   SUBROUTINE CENT_ROD(BUF1,BUF2,NXLD,NX,NY,THETA,SHX,SHY,IRTFLG)     *
C                                                                      *
C   PARAMETERS: BUF1    DATA ARRAY                              INPUT  *
C               BUF2    DATA ARRAY                              OUTPUT *
C               SHX,SHY SHIFTS                                  OUTPUT *
C               THETA   ROTATION                                OUTPUT *
C               IRTFLG  ERROR FLAG                              OUTPUT *
C                                                                      *
C   PURPOSE: SHIFT AND ROTATE FILAMENTOUS/ROD OBJECT                   *
C   MOVES FILAMENT/ROD TO CENTER OF IMAGE AND ALIGN IT VERTICALLY.     *
C   USES 2D FOURIER-BASED BICUBIC SPLINE INTERPOLATION BETWEEN PIXELS. *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE CENT_ROD(BUF1,BUF2,
     &                     NXLD,NX,NY,
     &                     THETA,SHX,SHY, IRTFLG)

       IMPLICIT NONE

       REAL              :: BUF1(NXLD,NY)
       REAL              :: BUF2(NX,  NY)
       INTEGER           :: NXLD, NX, NY
       REAL              :: THETA,SHX,SHY
       INTEGER           :: IRTFLG

       INTEGER           :: INV, MWANT
       INTEGER           :: CX, CY
       INTEGER           :: K,I,J,I2,J2,I3,J3
       INTEGER           :: I4,J4,I5,J5
       INTEGER           :: ANG, NMAX
       INTEGER           :: NX25,NX75,NY25,NY75

       INTEGER           :: ix,iy
       real              :: fx,fy,xold,yold,ycod,ysid
       real              :: fy0,fy1,fy2,fx0,fx1,fx2,shypny,shxpnx

       real              :: f
       REAL              :: COSTH, SINTH
       REAL              :: X,Y, FICX,FJCY
       REAL              :: DX,DY
       REAL              :: SHX2, SHY2
       REAL              :: COR, TEMP
       REAL              :: WEIGHT
       REAL              :: A,B,C,D

       REAL              :: fbs2

       REAL, ALLOCATABLE :: F0(:,:)
       REAL, ALLOCATABLE :: F1(:,:)
       REAL, ALLOCATABLE :: X1(:,:)
       REAL, ALLOCATABLE :: Y1(:,:)
       REAL, ALLOCATABLE :: XY2(:,:)

       REAL, PARAMETER   :: PI7 = 3.14159265358979323846 /720.0


       ! ARRAY ALLOCATION
       ALLOCATE (F0 (0:NXLD-1,0:NY-1),
     &           F1 (0:NXLD-1,0:NY-1),
     &           X1 (0:NXLD-1,0:NY-1),
     &           Y1 (0:NXLD-1,0:NY-1),
     &           XY2(0:NXLD-1,0:NY-1),
     &           STAT=IRTFLG)

       IF (IRTFLG .NE. 0) THEN
          MWANT = 4* NXLD*NY
          CALL ERRT(46,'CENT_ROD, F0...',MWANT)
          RETURN
       ENDIF

       F0 = BUF1   !ARRAY OPERATION

       CALL FBS2_PREP(F0, X1,Y1,XY2, NXLD, NX,NY, IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

C      COMPUTATION OF SQUARE OF A COMPLEX NUMBER AT EACH POINT
C      OF 2D FOURIER TRANSFORM OF IMAGE AFTER A 180-DEGREE ROTATION
       DO J=0,NY-1
          DO I=0,NXLD/2-1
             F1(2*I,J)   =     F0(2*I,J)**2 - F0(2*I+1,J)**2
             F1(2*I+1,J) = 2 * F0(2*I,J)    * F0(2*I+1,J)
          ENDDO
       ENDDO

C      COMPUTATION OF REVERSED 2D FOURIER TRANSFORM THAT IS
C      EQUAL TO 2D CROSS-CORRELATION BETWEEN
C      ORIGINAL IMAGE AND IMAGE AFTER A 180-DEGREE ROTATION
       INV = -1
       CALL FMRS(F1, NX,NY,1, 0.0D0, .TRUE.,.TRUE., INV, IRTFLG)

C      SEARCH FOR SHIFT WITH MAXIMUM CORRELATION
C      BETWEEN ORIGINAL IMAGE AND IMAGE AFTER A 180-DEGREE ROTATION
       COR = 0.0
       SHX = 0
       SHY = 0

       DO J=0,NY-1
          DO I=0,NX-1
            IF (F1(I,J) > COR) THEN
              COR = F1(I,J)
              SHX = I
              SHY = J
            ENDIF
          ENDDO
       ENDDO

       IF (SHX > NX*0.5) SHX = SHX - NX
       IF (SHY > NY*0.5) SHY = SHY - NY

       SHX = (SHX + MOD(NX,2)) * 0.5
       SHY = (SHY + MOD(NY,2)) * 0.5

       CX  = NX / 2 + 1
       CY  = NY / 2 + 1

C      GENERATION OF IMAGE WITH OBJECT CENTERED BY OPTIMAL SHIFT
c$omp  parallel do private(j,j2,i,i2,x,y)
       DO J=0, NY-1
          J2 = J  + 1
          DO I=0, NX-1
             I2 =  I + 1
             X  = I2  + SHX
             Y  = J2  + SHY
             BUF2(I2,J2) = FBS2(X,Y, NXLD,NX,NY, BUF1,NXLD,
     &                          X1,Y1,XY2, .TRUE.)
          ENDDO
        ENDDO

C      CORRECTION OF CENTER OF GRAVITY:
C      IF THE OBJECT'S CENTER OF QUASISYMMETRY IS NOT IN
C      CENTRAL QUARTER - MOVE IT BY HALF-SIZE OF IMAGE;
C      FOUR QUARTERS ARE TESTED.
       SHX2 = SHX
       SHY2 = SHY

       NX25 = INT(NX * 0.25)
       NX75 = INT(NX * 0.75)
       NY25 = INT(NY * 0.25)
       NY75 = INT(NY * 0.75)

       A = 0
       DO J=NY25, NY75
          DO I=NX25, NX75
             A = A + BUF2(I,J)
          ENDDO
       ENDDO
       WEIGHT = A

       B = 0
       DO J= 1,NY25
          DO I= NX25, NX75
             B = B + BUF2(I,J)
          ENDDO
       ENDDO
       DO J= NY75,NY
          DO I= NX25, NX75
             B = B + BUF2(I,J)
          ENDDO
       ENDDO

       IF (B > WEIGHT) THEN
          WEIGHT = B
          SHX2 = SHX
          SHY2 = SHY - NY*0.5
       ENDIF

       C = 0
       DO J=NY25, NY75
          DO I=1, NX25
             C = C + BUF2(I,J)
          ENDDO
       ENDDO
       DO J= NY25, NY75
          DO I= NX75, NX
             C = C + BUF2(I,J)
          ENDDO
       ENDDO

       IF (C > WEIGHT) THEN
          WEIGHT = C
          SHX2   = SHX - NX*0.5
          SHY2   = SHY
       ENDIF

       D = 0
       DO J= 1,NY25
          DO I=1, NX25
             D = D + BUF2(I,J)
          ENDDO
       ENDDO

       DO J= 1,NY25
          DO I= NX75, NX
             D = D + BUF2(I,J)
          ENDDO
       ENDDO

       DO J= NY75,NY
          DO I=1, NX25
             D = D + BUF2(I,J)
          ENDDO
       ENDDO

       DO J= NY75,NY
          DO I= NX75, NX
             D = D + BUF2(I,J)
          ENDDO
       ENDDO

       IF (D > WEIGHT) THEN
          SHX2 = SHX - NX*0.5
          SHY2 = SHY - NY*0.5
       ENDIF

C      SEARCH FOR ROTATION ANGLE WITH MAXIMUM CORRELATION BETWEEN
C      ORIGINAL IMAGE ROTATED BY ANGLE THETA
C      AND ORIGINAL IMAGE ROTATED BY ANGLE THETA
C      AND THEN REFLECTED ALONG VERTICAL AXIS
       F0 = BUF1

       CX   = NX / 2 + 1
       CY   = NY / 2 + 1
       COR  = 0.0
       ANG  = 0.0

       NMAX = INT (ATAN(FLOAT(NX)/FLOAT(NY)) / PI7)

       F0 = BUF1   ! ARRAY ASSIGNMENT

       DO K=-NMAX, NMAX
          THETA = K * PI7
          COSTH = COS(THETA)
          SINTH = SIN(THETA)

          DO J=0, NY-1
             J2   = J  + 1
             FJCY =  J2 - CY + SHY
             DO I=0, NX-1
               I2   = I + 1
               FICX = I2 - CX + SHX
               X    = COSTH*FICX - SINTH*FJCY + CX
               Y    = SINTH*FICX + COSTH*FJCY + CY
               I3   = FLOOR(X)
               J3   = FLOOR(Y)
               DX   = X - I3
               DY   = Y - J3
               I4   = MODULO(I3-1, NX) + 1
               J4   = MODULO(J3-1, NY) + 1
               I5   = MODULO(I3, NX) + 1
               J5   = MODULO(J3, NY) + 1
               A    = BUF1(I4,J4)
               B    = BUF1(I5,J4)-BUF1(I4,J4)
               C    = BUF1(I4,J5)-BUF1(I4,J4)
               D    = BUF1(I4,J4)-BUF1(I5,J4)
     &               -BUF1(I4,J5)+BUF1(I5,J5)

               BUF2(I2,J2) = A + B*DX + C*DY + D*DX*DY
             ENDDO
          ENDDO

          TEMP = 0

          DO J3=1, NY
             DO I3=1, NX
               TEMP = TEMP + BUF2(I3,J3) * BUF2(NX-I3+1,J3)
             ENDDO
          ENDDO
          IF (TEMP > COR) THEN
             COR = TEMP
             ANG = K
          ENDIF
       ENDDO

       THETA = ANG * PI7
       COSTH = COS(THETA)
       SINTH = SIN(THETA)

      write(6,*), 'Optimal Angle =',THETA
      write(6,*), 'SHX, SHY = ', SHX, SHY

c$omp  parallel do private(j,j2,fjcy,i,i2,ficx,x,y)
       DO J=0, NY-1
          J2   = J  + 1
          FJCY = J2 - CY + SHY2

          DO I=0, NX-1
             I2   = I + 1
             FICX = I2 - CX + SHX2
             X    = COSTH*FICX - SINTH*FJCY + CX
             Y    = SINTH*FICX + COSTH*FJCY + CY

             BUF2(I2,J2) = FBS2(X,Y, NXLD,NX,NY,F0,
     &		                NXLD,X1,Y1,XY2, .TRUE.)

          ENDDO
       ENDDO


       THETA = ANG * 180.0 / 720.0
 
999    IF (ALLOCATED(F0))  DEALLOCATE(F0)
       IF (ALLOCATED(F1))  DEALLOCATE(F1)
       IF (ALLOCATED(X1))  DEALLOCATE(X1)
       IF (ALLOCATED(Y1))  DEALLOCATE(Y1)
       IF (ALLOCATED(XY2)) DEALLOCATE(XY2)

      write(6,*), 'Optimal Angle =',THETA
      write(6,*), 'SHX2, SHY2 = ', SHX2, SHY2
 
       END

#ifdef NEVER
c$omp  parallel do private(j,j2,fjcy,i,i2,ficx,x,y)
       DO J=0, NY-1
          J2   = J  + 1
          FJCY = J2 - CY + SHY

          DO I=0, NX-1
             I2   = I + 1
             FICX = I2 - CX + SHX
             X    = COSTH*FICX - SINTH*FJCY + CX
             Y    = SINTH*FICX + COSTH*FJCY + CY

             BUF2(I2,J2) = FBS2(X,Y, NXLD,NX,NY,F0,
     &		                NXLD,X1,Y1,XY2, .TRUE.)

          ENDDO
       ENDDO
#endif

#ifdef NEVER

C      GENERATION OF IMAGE WITH OBJECT CENTERED BY OPTIMAL
C      SHIFT AND ROTATION

       FY0    = - SHY - CY
       FY1    = - SHY + NY - CY
       FY2    = - SHY - NY - CY

       FX0    = - SHX - CX
       FX1    = - SHX + NX - CX
       FX2    = - SHX - NX - CX

       SHYPNY = SHY + NY
       SHXPNX = SHX + NX

c$omp  parallel do private(iy,fy,ycod,ysid, ix,fx, xold,yold)
       DO IY=1, NY
          FY = IY + FY0 
          IF ((IY-1) <  SHY)    FY = IY + FY1 
          IF ((IY-1) >= SHYPNY) FY = IY + FY2 

          YCOD =  COSTH * FY + CY
          YSID = -SINTH * FY + CX

          DO IX=1, NX
            FX = IX + FX0 
            IF ((IX-1) <  SHX)    FX = IX + FX1
            IF ((IX-1) >= SHXPNX) FX = IX + FX2

            XOLD    = COSTH * FX + YSID
            YOLD    = SINTH * FX + YCOD

            BUF2(IX,IY) = FBS2(XOLD,YOLD, NXLD,NX,NY, F0,NXLD,
     &                         X1,Y1,XY2, .TRUE.)
          ENDDO
       ENDDO
#endif
