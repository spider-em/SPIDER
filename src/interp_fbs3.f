C ++********************************************************************
C                                                                      *
C  INTERP_FBS3      NEW             SEPTEMBER 2011  GREGORY KISHCHENKO *                                                                           *
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
C   INTERP_FBS3(BUFIN,BUFOUT, NX, NY, NZ, NX2,NY2,NZ2, NXLD)           *                                                        *
C                                                                      *
C   PURPOSE: RESAMPLING OF 3D IMAGES BY FOUIER-BASED TRICUBIC SPLINE   *
C            INTERPOLATION BETWEEN VOXELS.                             *
C            ALGORITHM IS FAIRLY FAST AND PRESERVES FINE DETAILS       *
C            OF  IMAGES                                                *
C                                                                      *
C            SUBROUTINES FBS3_PREP and FUNCTION FBS3 are used          *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE INTERP_FBS3(BUFIN,BUFOUT,
     &                        NX,NY,NZ, NX2,NY2,NZ2, NXLD)

       IMPLICIT NONE

       INCLUDE 'CMBLOCK.INC'

       REAL              :: BUFIN (NXLD, NY,  NZ)
       REAL              :: BUFOUT(NX2,  NY2, NZ2)
       REAL, ALLOCATABLE :: XYZ(:,:,:)
       REAL, ALLOCATABLE :: X1 (:,:,:)
       REAL, ALLOCATABLE :: Y1 (:,:,:)
       REAL, ALLOCATABLE :: Z1 (:,:,:)
       REAL, ALLOCATABLE :: XY2(:,:,:)
       REAL, ALLOCATABLE :: XZ2(:,:,:)
       REAL, ALLOCATABLE :: YZ2(:,:,:)

       INTEGER           :: NXLD, NX, NY, NZ
       INTEGER           :: NX2,  NY2, NZ2
       INTEGER           :: I, J, L
       INTEGER           :: K1,K2,K3
       INTEGER           :: K11,K21,K31
       INTEGER           :: INV,MWANT,IRTFLG
       REAL              :: SCALEX,SCALEY,SCALEZ
       REAL              :: X, Y, Z

       REAL              :: fbs3

       WRITE(NOUT,*) ' Fourier based spline 3D interpolation'

       ALLOCATE (XYZ (NXLD, NY, NZ),
     &           X1  (NXLD, NY, NZ),
     &           Y1  (NXLD, NY, NZ),
     &           Z1  (NXLD, NY, NZ),
     &           XY2 (NXLD, NY, NZ),
     &           XZ2 (NXLD, NY, NZ),
     &           YZ2 (NXLD, NY, NZ), STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN 
           MWANT = 7*NXLD*NX*NY*NZ 
          CALL ERRT(46,'INTERP_FSB; XYZ...',MWANT)
          GOTO 9999
       ENDIF 

       SCALEX = FLOAT(NX) / FLOAT(NX2)
       SCALEY = FLOAT(NY) / FLOAT(NY2)
       SCALEZ = FLOAT(NZ) / FLOAT(NZ2)

       XYZ = BUFIN  ! ARRAY ASSIGNMENT

C      RETURNS XYZ = XYZ3

       CALL FBS3_PREP(XYZ, NXLD, NX, NY, NZ,
     &                X1, Y1, Z1, XY2, XZ2, YZ2)
       IF (IRTFLG .NE. 0) GOTO 9999

c$omp  parallel do private(k3,z,k31, k2,y,k21, k1,x,k11)
       DO K3 = 0,NZ2-1
          Z   = K3 * SCALEZ + 1
          K31 = K3 + 1

          DO K2 = 0,NY2-1
             Y   = K2 * SCALEY + 1
             K21 = K2 + 1

             DO K1 = 0,NX2-1
             X   = K1 * SCALEX + 1
             K11 = K1 + 1

             BUFOUT(K11,K21,K31) = FBS3(X,Y,Z,
     &              NXLD, NX, NY, NZ,
     &              BUFIN,NXLD, XYZ,
     &              X1, Y1, Z1,
     &              XY2,XZ2,YZ2)
             ENDDO
          ENDDO
       ENDDO
       IRTFLG = 0

9999   IF (ALLOCATED(XYZ)) DEALLOCATE(XYZ)
       IF (ALLOCATED(X1))  DEALLOCATE(X1)
       IF (ALLOCATED(Y1))  DEALLOCATE(Y1)
       IF (ALLOCATED(Z1))  DEALLOCATE(Z1)
       IF (ALLOCATED(XY2)) DEALLOCATE(XY2)
       IF (ALLOCATED(XZ2)) DEALLOCATE(XZ2)
       IF (ALLOCATED(YZ2)) DEALLOCATE(YZ2)

       END
