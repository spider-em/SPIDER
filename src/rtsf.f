C **********************************************************************
C                                                                      *
C RTSF       NEW                         MAY 2011  Gregory Kishchenko  * 
C            SPLINE                      JUL 2011  ArDean Leith        *
C            OMP                         OCT 2011  ArDean Leith        * 
C            USED FBS2                   MAY 2011  Gregory Kishchenko  *
C            SCALE BUG                   DEC 2011  Gregory Kishchenko  *
C            SPEEDUP                     JAN 2012  ArDean Leith        *
C            RTSF_PAD                    JAN 2012  ArDean Leith        *
C            RYE1 BUG                    MAY 2012  ArDean Leith        *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C   RTSF(BUF1,BUF2, NXLD,NX,NY,                                        *
C        THETAT,SCLI,SHX,SHY, IRTFLG)                                  *
C                                                                      *  
C   RTSF_PAD(XIMG,BUFOUT, NX,NY, NXP,NYP,                              *
C            THETA,SCLI,SHXI,SHYI, IRTFLG)                             *
C                                                                      *
C   PARAMETERS: BUF1         INPUT ARRAY                         INPUT *
C               BUF2         OUTPUT ARRAY                        OUTPUT*
C               THETAT       ANGLE                               INPUT *
C               SCLI         SCALE FACTOR                        INPUT *
C               SHXI,SHYI    SHIFTS                              INPUT *
C               IRTFLG       ERROR FLAG                          INPUT *
C                                                                      *
C   VARIABLES:  X1   - D/DX DERIVATIVE                                 *
C               Y1   - D/DY DERIVATIVE                                 *
C               XY2  - D2/DXDY DERIVATIVE                              *
C                                                                      *
C   PURPOSE:    Image rotation, shift, & scale                         *
C               Using 2D Fourier-based bicubic spline interpolation    *
c               between pixels.                                        *
C               Algorithm is fairly fast and preserves fine details    *
c               of  images                                             *
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

C      INPUT: FFT PADDED  OUTPUT: UNPADDED

       SUBROUTINE RTSF(BUF1, BUF2,
     &                 NXLD, NX,  NY,
     &                 THETAT,SCLI, SHXI,SHYI, IRTFLG)

       IMPLICIT NONE

       REAL              :: BUF1(NXLD,NY)
       REAL              :: BUF2(NX,  NY)
       INTEGER           :: NXLD, NX, NY
       REAL              :: THETAT,SCLI,SCL
       REAL              :: THETA, COSTH, SINTH
       REAL              :: SHXI,SHYI
       INTEGER           :: IRTFLG

       INTEGER           :: CX, CY, BNDI, BNDJ
       INTEGER           :: I,J,K,IX,IY
       INTEGER           :: I4,J4
       INTEGER           :: INV,MWANT
       REAL              :: FI,FJ
       REAL              :: XOLD,YOLD, FICX,FJCY
       REAL              :: SHX,SHY,YCOD,YSID
       REAL              :: FY,FX,FX0,FY0,FY1,FY2,FX1,FX2,SHYPNY,SHXPNX

       REAL              :: fbs2

       REAL, ALLOCATABLE :: F0(:,:)
       REAL, ALLOCATABLE :: X1(:,:)
       REAL, ALLOCATABLE :: Y1(:,:)
       REAL, ALLOCATABLE :: XY2(:,:)

       REAL, PARAMETER   :: PI = 3.14159265358979323846

       ALLOCATE (F0 (NXLD, NY),
     &           X1 (NXLD, NY),
     &           Y1 (NXLD, NY),
     &           XY2(NXLD, NY),
     &           STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN 
          MWANT = 4* NXLD*NY 
          CALL ERRT(46,'RTSF, F0...',MWANT)
          RETURN
       ENDIF 

       ! ARRAY ASSIGNMENT FOR PADDED ARRAY
       F0 = BUF1

       CALL FBS2_PREP(F0, X1,Y1,XY2, NXLD, NX,NY, IRTFLG)

C      CREATE TRANSFORMATION MATRIX
       THETA = THETAT * PI / 180
       COSTH = COS(THETA)
       SINTH = SIN(THETA)

C      SPIDER IMAGE CENTER
       CX    = NX / 2 + 1
       CY    = NY / 2 + 1

C      CONSTRAIN SHIFT TO CIRCULAR BOUNDARY
       SHX = MODULO(SHXI, FLOAT(NX))
       SHY = MODULO(SHYI, FLOAT(NY))

       IF (SCLI == 1.0) THEN !-----------------------------------------

         FY0    = - SHY - CY
         FY1    = - SHY + NY - CY
         FY2    = - SHY - NY - CY

         FX0    = - SHX - CX
         FX1    = - SHX + NX - CX
         FX2    = - SHX - NX - CX

         SHYPNY = SHY + NY
         SHXPNX = SHX + NX

c$omp    parallel do private(iy,fy,ycod,ysid, ix,fx, xold,yold)
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

              BUF2(IX,IY) = FBS2(XOLD,YOLD, NXLD,NX,NY, BUF1,NXLD,
     &                           X1,Y1,XY2, .TRUE.)
            ENDDO
          ENDDO

       ELSE !------------------------------------------- Needs speedup

         SCL   = 1 / SCLI

         IF (SCLI >= 1) THEN
            BNDI = INT(NX*(SCLI-1)/2)
            BNDJ = INT(NY*(SCLI-1)/2)
         ELSE
            BNDI = INT(NX*(SCLI-1)/2) - 1
            BNDJ = INT(NY*(SCLI-1)/2) - 1
         ENDIF

C        'RT SQ' OMP GAVE TIME: 3 ON 8 PROCS VS  6 ON 1
C        'RT SF' OMP GAVE TIME: 7 ON 8 PROCS VS 21 ON 1

c$omp    parallel do private(j,fj,fjcy,j4,i,fi,ficx,i4,xold,yold)
         DO J=BNDJ, NY*SCLI-BNDJ-1
            FJ   = MODULO(INT(J - BNDJ - SHY), NY) + BNDJ
            FJCY = FJ * SCL - CY + 1
            J4   = MODULO(J - BNDJ, NY) + 1

               DO I=BNDI, NX*SCLI-BNDI-1
                 FI   = MODULO(INT(I - BNDI- SHX), NX) + BNDI
                 FICX = FI * SCL - CX + 1
                 I4   = MODULO(I - BNDI, NX)+1

                 XOLD = COSTH*FICX - SINTH*FJCY + CX
                 YOLD = SINTH*FICX + COSTH*FJCY + CY

                 !if (i4 < 1 .or. i4 > nx .or. j4 < 1 .or. j4 > ny ) THEN
                 !   write(6,*) ' bad index:',i,j,i4,j4,x,y
                 !   stop
                 !endif

                 BUF2(I4,J4) = FBS2(XOLD,YOLD, NXLD,NX,NY, BUF1,NXLD,
     &                        X1,Y1,XY2, .TRUE.)
              ENDDO
          ENDDO
       ENDIF

       IF (ALLOCATED(F0))  DEALLOCATE(F0)
       IF (ALLOCATED(X1))  DEALLOCATE(X1)
       IF (ALLOCATED(Y1))  DEALLOCATE(Y1)
       IF (ALLOCATED(XY2)) DEALLOCATE(XY2)

       END

C      -------------------- RTSF_BACK ---------------------------------

C      INPUT: FFT PADDED  OUTPUT: UNPADDED

       SUBROUTINE RTSF_BACK(BUF1, BUF2,
     &                      NXLD, NX,  NY,
     &                      THETAT,SCLI,SHXI,SHYI, 
     &                      USEBACK,BACK,  IRTFLG)

       IMPLICIT NONE

       REAL              :: BUF1(NXLD,NY)
       REAL              :: BUF2(NX,  NY)
       INTEGER           :: NXLD, NX, NY
       REAL              :: THETAT,SCLI,SCL
       REAL              :: THETA,COSTH, SINTH
       REAL              :: SHXI,SHYI
       LOGICAL           :: USEBACK
       REAL              :: BACK
       INTEGER           :: IRTFLG

       INTEGER           :: CX, CY, BNDI, BNDJ
       INTEGER           :: I,J,K,IX,IY
       INTEGER           :: I4,J4
       INTEGER           :: INV,MWANT
       REAL              :: FI,FJ
       REAL              :: XOLD,YOLD
       REAL              :: FICX,FJCY
       REAL              :: SHX,SHY,YCOD,YSID
       REAL              :: FY,FX,FX0,FY0,FY1,FY2,FX1,FX2,SHYPNY,SHXPNX

       REAL              :: fbs2

       REAL, ALLOCATABLE :: F0 (:,:)
       REAL, ALLOCATABLE :: X1 (:,:)
       REAL, ALLOCATABLE :: Y1 (:,:)
       REAL, ALLOCATABLE :: XY2(:,:)

       REAL, PARAMETER   :: PI = 3.14159265358979323846

       ALLOCATE (F0 (NXLD, NY),
     &           X1 (NXLD, NY),
     &           Y1 (NXLD, NY),
     &           XY2(NXLD, NY),
     &           STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN 
          MWANT = 4* NXLD*NY 
          CALL ERRT(46,'RTSF_BACK, F0...',MWANT)
          RETURN
       ENDIF 

       ! ARRAY ASSIGNMENT
       F0  = BUF1

       CALL FBS2_PREP(F0, X1,Y1,XY2, NXLD, NX,NY, IRTFLG)

C      CREATE TRANSFORMATION MATRIX
       THETA = THETAT * PI / 180
       COSTH = COS(THETA)
       SINTH = SIN(THETA)

C      SPIDER IMAGE CENTER
       CX    = NX / 2 +1
       CY    = NY / 2 +1

C      CONSTRAIN SHIFT TO CIRCULAR BOUNDARY
       SHX = MODULO(SHXI, FLOAT(NX))
       SHY = MODULO(SHYI, FLOAT(NY))

       IF (SCLI == 1.0) THEN !-----------------------------------------

         FY0    = - SHY - CY
         FY1    = - SHY + NY - CY
         FY2    = - SHY - NY - CY

         FX0    = - SHX - CX
         FX1    = - SHX + NX - CX
         FX2    = - SHX - NX - CX

         SHYPNY = SHY + NY
         SHXPNX = SHX + NX

c$omp    parallel do private(iy,fy,ycod,ysid, ix,fx, xold,yold)
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

              IF  (USEBACK .AND. 
     &            (XOLD < 1 .OR. XOLD > NX .OR. 
     &             YOLD < 1 .OR. YOLD > NY)) THEN
C                CORNER LOCATION IN NEW IMAGE
                 BUF2(IX,IY) = BACK     
                 CYCLE
              ENDIF

              BUF2(IX,IY) = FBS2(XOLD,YOLD, NXLD,NX,NY, BUF1,NXLD,
     &                             X1,Y1,XY2, .TRUE.)
            ENDDO
          ENDDO
      
       ELSE !------------------------------------------ Needs speedup

          SCL   = 1 / SCLI

          IF (SCLI >= 1) THEN
             BNDI = INT(NX*(SCLI-1)/2)
             BNDJ = INT(NY*(SCLI-1)/2)
          ELSE
             BNDI = INT(NX*(SCLI-1)/2) -1
             BNDJ = INT(NY*(SCLI-1)/2) -1
          ENDIF

C         'RT SQ' OMP GAVE TIME:   ON 8 PROCS VS    ON 1
C         'RT SF' OMP GAVE TIME:   ON 8 PROCS VS    ON 1

c$omp     parallel do private(j,fj,fjcy,j4, i,fi,ficx,i4, xold,yold)
          DO J=BNDJ, NY*SCLI-BNDJ-1

             FJ   = MODULO(INT(J - BNDJ- SHY), NY) + BNDJ
             FJCY = FJ*SCL - CY + 1
             J4   = MODULO(J - BNDJ, NY)+1

              DO I=BNDI, NX*SCLI-BNDI-1
                 FI   = MODULO(INT(I - BNDI- SHX), NX) + BNDI
                 FICX = FI * SCL - CX + 1
                 I4   = MODULO(I - BNDI, NX) + 1

                 XOLD = COSTH*FICX - SINTH*FJCY + CX
                 YOLD = SINTH*FICX + COSTH*FJCY + CY

                 IF  (USEBACK .AND. 
     &               (XOLD < 1 .OR. XOLD > NX .OR. 
     &                YOLD < 1 .OR. YOLD > NY)) THEN
C                   CORNER LOCATION IN NEW IMAGE
                    BUF2(I4,J4) = BACK     
                    CYCLE
                 ENDIF

                 BUF2(I4,J4) = FBS2(XOLD,YOLD, NXLD,NX,NY, BUF1,NXLD,
     &                           X1,Y1,XY2, .TRUE.)
             ENDDO
          ENDDO
       ENDIF

       IF (ALLOCATED(F0))  DEALLOCATE(F0)
       IF (ALLOCATED(X1))  DEALLOCATE(X1)
       IF (ALLOCATED(Y1))  DEALLOCATE(Y1)
       IF (ALLOCATED(XY2)) DEALLOCATE(XY2)

       END


C      -------------------- RTSF_PADIN ---------------------------------

C      INPUT: FFT or 2xFFT PADDED,   OUTPUT: UNPADDED

       SUBROUTINE RTSF_PADIN(BUF1, BUF2,
     &                 NXLD, NX,  NY, NXP,
     &                 THETAT,SCLI, SHXI,SHYI, IRTFLG)

       IMPLICIT NONE

       REAL              :: BUF1(NXP,NY)
       REAL              :: BUF2(NX, NY)
       INTEGER           :: NXLD, NX, NY, NXP
       REAL              :: THETAT,SCLI,SCL
       REAL              :: THETA, COSTH, SINTH
       REAL              :: SHXI,SHYI
       INTEGER           :: IRTFLG

       INTEGER           :: CX, CY, BNDI, BNDJ
       INTEGER           :: I,J,K,IX,IY
       INTEGER           :: I4,J4
       INTEGER           :: INV,MWANT
       REAL              :: FI,FJ
       REAL              :: XOLD,YOLD, FICX,FJCY
       REAL              :: SHX,SHY,YCOD,YSID
       REAL              :: FY,FX,FX0,FY0,FY1,FY2,FX1,FX2,SHYPNY,SHXPNX

       REAL              :: fbs2

       REAL, ALLOCATABLE :: F0(:,:)
       REAL, ALLOCATABLE :: X1(:,:)
       REAL, ALLOCATABLE :: Y1(:,:)
       REAL, ALLOCATABLE :: XY2(:,:)

       REAL, PARAMETER   :: PI = 3.14159265358979323846

       ALLOCATE (F0 (NXLD, NY),
     &           X1 (NXLD, NY),
     &           Y1 (NXLD, NY),
     &           XY2(NXLD, NY),
     &           STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN 
          MWANT = 4* NXLD*NY 
          CALL ERRT(46,'RTSF, F0...',MWANT)
          RETURN
       ENDIF 

       ! ARRAY ASSIGNMENT FOR PADDED ARRAY
       F0 = BUF1(1:NXLD,1:NY)

       CALL FBS2_PREP(F0, X1,Y1,XY2, NXLD, NX,NY, IRTFLG)

C      CREATE TRANSFORMATION MATRIX
       THETA = THETAT * PI / 180
       COSTH = COS(THETA)
       SINTH = SIN(THETA)

C      SPIDER IMAGE CENTER
       CX    = NX / 2 + 1
       CY    = NY / 2 + 1

C      CONSTRAIN SHIFT TO CIRCULAR BOUNDARY
       SHX = MODULO(SHXI, FLOAT(NX))
       SHY = MODULO(SHYI, FLOAT(NY))

       IF (SCLI == 1.0) THEN !-----------------------------------------

         FY0    = - SHY - CY
         FY1    = - SHY + NY - CY
         FY2    = - SHY - NY - CY

         FX0    = - SHX - CX
         FX1    = - SHX + NX - CX
         FX2    = - SHX - NX - CX

         SHYPNY = SHY + NY
         SHXPNX = SHX + NX

c$omp    parallel do private(iy,fy,ycod,ysid, ix,fx, xold,yold)
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

              BUF2(IX,IY) = FBS2(XOLD,YOLD, NXLD,NX,NY, BUF1,NXP,
     &                           X1,Y1,XY2, .TRUE.)
            ENDDO
          ENDDO

       ELSE !------------------------------------------- Needs speedup

         SCL   = 1 / SCLI

         IF (SCLI >= 1) THEN
            BNDI = INT(NX*(SCLI-1)/2)
            BNDJ = INT(NY*(SCLI-1)/2)
         ELSE
            BNDI = INT(NX*(SCLI-1)/2) - 1
            BNDJ = INT(NY*(SCLI-1)/2) - 1
         ENDIF

C        'RT SQ' OMP GAVE TIME: 3 ON 8 PROCS VS  6 ON 1
C        'RT SF' OMP GAVE TIME: 7 ON 8 PROCS VS 21 ON 1

c$omp    parallel do private(j,fj,fjcy,j4,i,fi,ficx,i4,xold,yold)
         DO J=BNDJ, NY*SCLI-BNDJ-1
            FJ   = MODULO(INT(J - BNDJ - SHY), NY) + BNDJ
            FJCY = FJ * SCL - CY + 1
            J4   = MODULO(J - BNDJ, NY) + 1

               DO I=BNDI, NX*SCLI-BNDI-1
                 FI   = MODULO(INT(I - BNDI- SHX), NX) + BNDI
                 FICX = FI * SCL - CX + 1
                 I4   = MODULO(I - BNDI, NX)+1

                 XOLD = COSTH*FICX - SINTH*FJCY + CX
                 YOLD = SINTH*FICX + COSTH*FJCY + CY

                 !if (i4 < 1 .or. i4 > nx .or. j4 < 1 .or. j4 > ny ) THEN
                 !   write(6,*) ' bad index:',i,j,i4,j4,x,y
                 !   stop
                 !endif

                 BUF2(I4,J4) = FBS2(XOLD,YOLD, NXLD,NX,NY, BUF1,NXLD,
     &                        X1,Y1,XY2, .TRUE.)
              ENDDO
          ENDDO
       ENDIF

       IF (ALLOCATED(F0))  DEALLOCATE(F0)
       IF (ALLOCATED(X1))  DEALLOCATE(X1)
       IF (ALLOCATED(Y1))  DEALLOCATE(Y1)
       IF (ALLOCATED(XY2)) DEALLOCATE(XY2)

       END



C******************************** RTSF_PAD ****************************

C        INPUT: UNPADDED  OUTPUT: PADDED 

         SUBROUTINE RTSF_PAD(XIMG,BUFOUT, NX,NY, NXP,NYP,
     &                       THETA,SCLI,SHXI,SHYI, IRTFLG)

         IMPLICIT NONE
         REAL              :: XIMG(NX,NY)
         REAL              :: BUFOUT(NXP,NYP)
         INTEGER           :: NX,NY, NXP,NYP
         REAL              :: THETA,SCLI,SHXI,SHYI
         INTEGER           :: IRTFLG

         REAL              :: SHX,SHY,RY1,RX1,RY2,RX2,COD,SID,XI
         REAL              :: CODDSCLI,SIDDSCLI,FIXCENMSHX,FIYCENMSHY 
         REAL              :: RYE2,RYE1,RXE2,RXE1,YI
         REAL              :: YCOD,YSID,YOLD,XOLD
         INTEGER           :: IYCEN,IXCEN,IX,IY,MWANT,NXLD

         LOGICAL           :: CHKBOUND = .TRUE.

         REAL, ALLOCATABLE :: F0 (:,:)
         REAL, ALLOCATABLE :: X1 (:,:)
         REAL, ALLOCATABLE :: Y1 (:,:)
         REAL, ALLOCATABLE :: XY2(:,:)

	 REAL, PARAMETER   :: QUADPI = 3.14159265358979323846
	 REAL, PARAMETER   :: DGR_TO_RAD = (QUADPI/180)

         REAL              :: fbs2,quadri

         NXLD = NX + 2 - MOD(NX,2)       ! PAD FOR FFTW

         ALLOCATE (F0 (NXLD, NY),
     &             X1 (NXLD, NY),
     &             Y1 (NXLD, NY),
     &             XY2(NXLD, NY),
     &             STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            MWANT = 4* NXLD*NY 
            CALL ERRT(46,'RTSF_PAD, F0...',MWANT)
            RETURN
         ENDIF 

         ! ARRAY ASSIGNMENT
         F0(1:NX,1:NY) = XIMG(1:NX,1:NY)

         CALL FBS2_PREP(F0, X1,Y1,XY2, NXLD, NX,NY, IRTFLG)

C        SHIFT WITHIN IMAGE BOUNDARY
         SHX = MOD(SHXI,FLOAT(NX))
         SHY = MOD(SHYI,FLOAT(NY))

C        SPIDER IMAGE CENTER
         IXCEN = NX/2+1
         IYCEN = NY/2+1

C        IMAGE DIMENSIONS AROUND ORIGIN
         RX1   = -NX/2
         RX2   =  NX/2
         RY1   = -NY/2
         RY2   =  NY/2

         IF (MOD(NX,2) == 0) THEN
            RX2  =  RX2 - 1.0
            RXE1 = -NX
            RXE2 =  NX
         ELSE
            RXE1 = -NX - 1
            RXE2 =  NX + 1
         ENDIF

         IF (MOD(NY,2) == 0) THEN
            RY2  =  RY2 - 1.0
            RYE1 = -NY 
            RYE2 =  NY
         ELSE
            RYE1 = -NY - 1
            RYE2 =  NY + 1

         ENDIF

C        CREATE TRANSFORMATION MATRIX
         COD = COS(THETA * DGR_TO_RAD)
         SID = SIN(THETA * DGR_TO_RAD)

C        ADJUST FOR SCALING
         CODDSCLI = COD / SCLI
         SIDDSCLI = SID / SCLI

C        -(CENTER PLUS SHIFT)
         FIXCENMSHX = -IXCEN - SHX
         FIYCENMSHY = -IYCEN - SHY

c$omp    parallel do private(iy,yi,ycod,ysid, ix,xi,xold,yold)
         DO IY=1,NY
            YI = IY + FIYCENMSHY
            IF (YI < RY1) YI = MIN(YI+RYE2, RY2)
            IF (YI > RY2) YI = MAX(YI+RYE1, RY1)

            YCOD =  YI * CODDSCLI + IYCEN
            YSID = -YI * SIDDSCLI + IXCEN

            DO IX=1,NX
               XI = IX + FIXCENMSHX                           
               IF (XI  <  RX1) XI = MIN(XI+RXE2, RX2)   
               IF (XI  >  RX2) XI = MAX(XI+RXE1, RX1) 
 
               YOLD          = XI * SIDDSCLI + YCOD  
               XOLD          = XI * CODDSCLI + YSID 
 
               BUFOUT(IX,IY) = FBS2(XOLD,YOLD, NXLD,NX,NY, 
     &                              XIMG,NX,
     &                              X1,Y1,XY2, CHKBOUND)
            ENDDO
         ENDDO

         !call chkfile('jnk-rot-1',98,1,nx,ny,1, bufout,irtflg)

         IF (ALLOCATED(F0))  DEALLOCATE(F0)
         IF (ALLOCATED(X1))  DEALLOCATE(X1)
         IF (ALLOCATED(Y1))  DEALLOCATE(Y1)
         IF (ALLOCATED(XY2)) DEALLOCATE(XY2)

         IRTFLG = 0

         END



