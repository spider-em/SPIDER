C++*********************************************************************
C                                                                      *
C WPRO_FBS.F  FROM: WPRO_N                      DEC 2011 G. KISHCHENKO *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C WPRO_FBS(PROJ, NXLD,NX,NY,NZ, CUBE,                                  *
C          X1,Y1,Z1,XY2,XZ2,YZ2,XYZ,  IPCUBE, NN,                      *
C          PHI,THETA,PSI, RI, LDPX,LDPY,LDPZ)                          * 
C                                                                      *
C PURPOSE:  COMPUTES PROJECTION(S) OF VOLUME ACCORDING                 *
C           THREE EULERIAN ANGLES.                                     *
C           DOES A WHOLE PROJECTION SERIES.                            *
C           USES 3D FBS INTERPOLATION.                                 *
C                                                                      *
C IPCUBE: A  LIST OF VOXELS ON EACH LINE IN THE              *
C         VOLUME  WHICH ARE WITHIN A SPECIFED RADIUS SPHERE IN VOL.    *
C                1 - BEGINNING VOXEL ON LINE                           *
C                2 - ENDING    VOXEL ON LINE                                     *
C                3 - IX FOR VOXEL                                      *
C                4 - IY FOR VOXEL                                      *
C                5 - IZ FOR VOXEL                                      *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE WPRO_FBS(PROJ, NXLD,NX,NY,NZ, CUBE,
     &                       X1,Y1,Z1,XY2,XZ2,YZ2,XYZ,
     &                       IPCUBE, NN,
     &                       PHI,THETA,PSI, RI, LDPX,LDPY,LDPZ)

         IMPLICIT NONE

         REAL                :: PROJ(NX,NY)
         INTEGER             :: NXLD,NX,NY,NZ
         REAL                :: CUBE(NX,NY,NZ)
         REAL                :: X1 (NXLD,NY,NZ),Y1 (NXLD,NY,NZ)
         REAL                :: Z1 (NXLD,NY,NZ)
         REAL                :: XY2(NXLD,NY,NZ),XZ2(NXLD,NY,NZ)
         REAL                :: YZ2(NXLD,NY,NZ),XYZ(NXLD,NY,NZ)
         INTEGER             :: IPCUBE(5,NN)
         INTEGER             :: NN
         REAL                :: PHI,THETA,PSI,RI
         INTEGER             :: LDPX,LDPY,LDPZ

         REAL                :: DM(9)
         DOUBLE PRECISION    :: CPHI,SPHI,CTHE,STHE
         DOUBLE PRECISION    :: CPSI,SPSI

         REAL                :: fbs3

         INTEGER             :: K,J,I,IRTFLG,IOX,IOY,IOZ
         REAL                :: DIM,XB,YB,ZB,DM1,DM2,DM3

         DOUBLE PRECISION,PARAMETER :: QUADPI = 3.1415926535897932384626
         DOUBLE PRECISION,PARAMETER :: DGR_TO_RAD = (QUADPI/180)

         CPHI = DCOS(DBLE(PHI)  *DGR_TO_RAD)
         SPHI = DSIN(DBLE(PHI)  *DGR_TO_RAD)

         CTHE = DCOS(DBLE(THETA)*DGR_TO_RAD)
         STHE = DSIN(DBLE(THETA)*DGR_TO_RAD)

         CPSI = DCOS(DBLE(PSI)  *DGR_TO_RAD)
         SPSI = DSIN(DBLE(PSI)  *DGR_TO_RAD)

         DM(1) =  CPHI*CTHE*CPSI - SPHI*SPSI
         DM(2) =  SPHI*CTHE*CPSI + CPHI*SPSI
         DM(3) = -STHE*CPSI

         DM(4) = -CPHI*CTHE*SPSI - SPHI*CPSI
         DM(5) = -SPHI*CTHE*SPSI + CPHI*CPSI
         DM(6) =  STHE*SPSI

         DM(7) =  STHE*CPHI
         DM(8) =  STHE*SPHI
         DM(9) =  CTHE

         DM1   = DM(1)
         DM2   = DM(2)
         DM3   = DM(3)

C        ZERO THE WHOLE PROJ ARRAY
         PROJ  = 0.0

         DIM  = MIN(NX,NY,NZ)

         IF ((2*(RI+1) + 1) .LE. DIM) THEN
C            NO NEED TO CHECK IQX & IQY BOUNDARIES

             DO I=1,NN

               K  = IPCUBE(4,I)+LDPY

               XB = IPCUBE(3,I)*DM(1) + IPCUBE(4,I)*DM(4) +
     &              IPCUBE(5,I)*DM(7) + LDPX

               YB = IPCUBE(3,I)*DM(2) + IPCUBE(4,I)*DM(5) +
     &              IPCUBE(5,I)*DM(8) + LDPY

               ZB = IPCUBE(3,I)*DM(3) + IPCUBE(4,I)*DM(6) +
     &              IPCUBE(5,I)*DM(9) + LDPZ

               DO J=IPCUBE(1,I),IPCUBE(2,I)  ! OVER VOXELS ON THIS LINE

                  PROJ(J,K) = PROJ(J,K) + 
     &                     FBS3(XB,YB,ZB,
     &                     NXLD,  NX, NY, NZ,
     &                     CUBE,NX, XYZ,X1,Y1,Z1,
     &                     XY2,XZ2,YZ2)

                  XB = XB + DM1
                  YB = YB + DM2
                  ZB = ZB + DM3
               ENDDO
            ENDDO
         ELSE
C            MUST CHECK IQX & IQY BOUNDARIES

             DO I=1,NN

               K  = IPCUBE(4,I)+LDPY

               XB = IPCUBE(3,I)*DM(1) + IPCUBE(4,I)*DM(4) +
     &              IPCUBE(5,I)*DM(7) + LDPX

               YB = IPCUBE(3,I)*DM(2) + IPCUBE(4,I)*DM(5) +
     &              IPCUBE(5,I)*DM(8) + LDPY

               ZB = IPCUBE(3,I)*DM(3) + IPCUBE(4,I)*DM(6) +
     &              IPCUBE(5,I)*DM(9) + LDPZ

               DO J=IPCUBE(1,I),IPCUBE(2,I)

C                 CHECK FOR PIXELS OUT OF BOUNDS 
                  IOX = IFIX(XB)
                  IOY = IFIX(YB)
                  IOZ = IFIX(ZB)

                  IF ((IOX > 0 .AND. IOX <= NX) .AND.
     &                (IOY > 0 .AND. IOY <= NY) .AND.
     &                (IOZ > 0 .AND. IOZ <= NZ)) THEN

                       PROJ(J,K) = PROJ(J,K) + 
     &                          FBS3(XB,YB,ZB,
     &                          NXLD, NX, NY, NZ,
     &                          CUBE,NX, XYZ,X1,Y1,Z1,
     &                          XY2,XZ2,YZ2)
                  ENDIF

                  XB = XB + DM1
                  YB = YB + DM2
                  ZB = ZB + DM3
               ENDDO
            ENDDO
          ENDIF

          END
