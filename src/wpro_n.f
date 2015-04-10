C++*********************************************************************
C                                                                      *
C WPRO_N.F    SPEEDED UP                         FEB 2000 ARDEAN LEITH *
C             LDPX,LDPY,LDPZ                     NOV 2011 ARDEAN LEITH *
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
C WPRO_N(B,NSAM,NROW,CUBE,LTB,IPCUBE,NN,                               *
C        PHI,THETA,PSI,RI,LDPX,LDPY,LDPZ)                              *
C                                                                      *
C PURPOSE: COMPUTES PROJECTION(S) OF VOLUME ACCORDING TO EULERIAN ANGLES. 
C                                                                      *
C IPCUBE: A RUN LENGTH LIST OF VOXELS ON EACH LINE IN THE              *
C         VOLUME  WHICH ARE WITHIN A SPECIFED RADIUS SPHERE IN VOL.    *
C                1 - BEGINNING VOXEL ON LINE                           *
C                2 - LENGTH OF RUN                                     *
C                3 - IX FOR VOXEL                                      *
C                4 - IY FOR VOXEL                                      *
C                5 - IZ FOR VOXEL                                      *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE WPRO_N(B,NSAM,NROW,NSLICE,CUBE,LTB,IPCUBE,
     &                     NN,PHI,THETA,PSI,RI,LDPX,LDPY,LDPZ)

         REAL,DIMENSION(NSAM,NROW)        :: B
         INTEGER                          :: NSAM,NROW,NSLICE
         REAL,DIMENSION(NSAM,NROW,NSLICE) :: CUBE
         INTEGER                          :: LTB
         INTEGER,DIMENSION(5,NN)          :: IPCUBE
         INTEGER                          :: NN
         REAL                             :: PHI,THETA,PSI,RI
         INTEGER                          :: LDPX,LDPY,LDPZ

         REAL,DIMENSION(9)                :: DM
         DOUBLE PRECISION                 :: CPHI,SPHI,CTHE,STHE
         DOUBLE PRECISION                 :: CPSI,SPSI

         DOUBLE PRECISION                 :: QUADPI,DGR_TO_RAD

         PARAMETER (QUADPI = 3.1415926535897932384626)
         PARAMETER (DGR_TO_RAD = (QUADPI/180))

         CPHI = DCOS(DBLE(PHI)  *DGR_TO_RAD)
         SPHI = DSIN(DBLE(PHI)  *DGR_TO_RAD)
         CTHE = DCOS(DBLE(THETA)*DGR_TO_RAD)
         STHE = DSIN(DBLE(THETA)*DGR_TO_RAD)
         CPSI = DCOS(DBLE(PSI)  *DGR_TO_RAD)
         SPSI = DSIN(DBLE(PSI)  *DGR_TO_RAD)

         DM(1) = CPHI*CTHE*CPSI-SPHI*SPSI
         DM(2) = SPHI*CTHE*CPSI+CPHI*SPSI
         DM(3) = -STHE*CPSI
         DM(4) = -CPHI*CTHE*SPSI-SPHI*CPSI
         DM(5) = -SPHI*CTHE*SPSI+CPHI*CPSI
         DM(6) = STHE*SPSI

         DM(7) = STHE*CPHI
         DM(8) = STHE*SPHI
         DM(9) = CTHE

         DM1   = DM(1)
         DM2   = DM(2)
         DM3   = DM(3)

C        ZERO THE WHOLE B ARRAY
         B  = 0.0

         DIM  = MIN(NSAM,NROW,NSLICE)

         IF ((2*(RI+1) + 1) .LE. DIM) THEN
C            NO NEED TO CHECK IQX & IQY BOUNDARIES
             DO I=1,NN

               K=IPCUBE(4,I)+LDPY

               XB = IPCUBE(3,I)*DM(1) + IPCUBE(4,I)*DM(4) +
     &              IPCUBE(5,I)*DM(7) + LDPX

               YB = IPCUBE(3,I)*DM(2) + IPCUBE(4,I)*DM(5) +
     &              IPCUBE(5,I)*DM(8) + LDPY

               ZB = IPCUBE(3,I)*DM(3) + IPCUBE(4,I)*DM(6) +
     &              IPCUBE(5,I)*DM(9) + LDPZ

               DO J=IPCUBE(1,I),IPCUBE(2,I)

                 IOX = IFIX(XB)
                 IOY = IFIX(YB)
                 IOZ = IFIX(ZB)
		 DX  = XB-IOX
		 DY  = YB-IOY
		 DZ  = ZB-IOZ

                 A1 = CUBE(IOX,IOY,IOZ)
                 A2 = CUBE(IOX+1,IOY,IOZ) - A1
                 A3 = CUBE(IOX,IOY+1,IOZ) - A1
                 A4 = CUBE(IOX,IOY,IOZ+1) - A1
                 A5 = -A2 - CUBE(IOX,IOY+1,IOZ) + CUBE(IOX+1,IOY+1,IOZ)
                 A61= - CUBE(IOX,IOY,IOZ+1) + CUBE(IOX+1,IOY,IOZ+1)
                 A6 = -A2 + A61
                 A7 = -A3 - CUBE(IOX,IOY,IOZ+1) + CUBE(IOX,IOY+1,IOZ+1)
                 A8 = -A5 - A61 - CUBE(IOX,IOY+1,IOZ+1) + 
     &                CUBE(IOX+1,IOY+1,IOZ+1)

                 B(J,K) = B(J,K) +
     &              A1 + DZ*(A4 + A6*DX + (A7 + A8*DX)*DY) + A3*DY +
     &              DX*(A2 + A5*DY)

                  XB = XB + DM1
                  YB = YB + DM2
                  ZB = ZB + DM3
               ENDDO
            ENDDO
         ELSE
C            MUST CHECK IQX & IQY BOUNDARIES

            DO    I=1,NN

               K=IPCUBE(4,I)+LDPY

               XB = IPCUBE(3,I)*DM(1) + IPCUBE(4,I)*DM(4) +
     &              IPCUBE(5,I)*DM(7) + LDPX

               YB = IPCUBE(3,I)*DM(2) + IPCUBE(4,I)*DM(5) +
     &              IPCUBE(5,I)*DM(8) + LDPY

               ZB = IPCUBE(3,I)*DM(3) + IPCUBE(4,I)*DM(6) +
     &              IPCUBE(5,I)*DM(9) + LDPZ

               DO J=IPCUBE(1,I),IPCUBE(2,I)

C                 CHECK FOR PIXELS OUT OF BOUNDS 
                  IOX = IFIX(XB)
                  IF (IOX.LT.1 .OR. IOX.GE.NSAM) GOTO 2
                  IOY = IFIX(YB)
                  IF (IOY.LT.1 .OR. IOY.GE.NROW) GOTO 2
                  IOZ    = IFIX(ZB)
                  IF (IOZ.LT.1 .OR. IOZ.GE.NSLICE) GOTO 2

                  IOX = IFIX(XB)
                  IOY = IFIX(YB)
                  IOZ = IFIX(ZB)
		  DX  = XB-IOX
		  DY  = YB-IOY
		  DZ  = ZB-IOZ

                  A1 = CUBE(IOX,IOY,IOZ)
                  A2 = CUBE(IOX+1,IOY,IOZ) - A1
                  A3 = CUBE(IOX,IOY+1,IOZ) - A1
                  A4 = CUBE(IOX,IOY,IOZ+1) - A1
                  A5 = -A2 - CUBE(IOX,IOY+1,IOZ) +CUBE(IOX+1,IOY+1,IOZ)
                  A61= -CUBE(IOX,IOY,IOZ+1) + CUBE(IOX+1,IOY,IOZ+1)
                  A6 = -A2 + A61
                  A7 = -A3 - CUBE(IOX,IOY,IOZ+1) +CUBE(IOX,IOY+1,IOZ+1)
                  A8 = -A5 - A61 - CUBE(IOX,IOY+1,IOZ+1) + 
     &                CUBE(IOX+1,IOY+1,IOZ+1) 

                  B(J,K) = B(J,K)+ A1 + DZ*(A4 + A6*DX +
     &                     (A7 + A8*DX)*DY) + A3*DY + DX*(A2 + A5*DY)

2                 XB = XB + DM1
                  YB = YB + DM2
                  ZB = ZB + DM3
               ENDDO
            ENDDO
          ENDIF

         END

#ifdef NEVER
         A1 = CUBE(IOX,IOY,IOZ)
         A2 = CUBE(IOX+1,IOY,IOZ) - CUBE(IOX,IOY,IOZ)
         A3 = CUBE(IOX,IOY+1,IOZ) - CUBE(IOX,IOY,IOZ)
         A4 = CUBE(IOX,IOY,IOZ+1) - CUBE(IOX,IOY,IOZ)
         A5 = CUBE(IOX,IOY,IOZ) - CUBE(IOX+1,IOY,IOZ) 
     &   - CUBE(IOX,IOY+1,IOZ)+ CUBE(IOX+1,IOY+1,IOZ)
         A6 = CUBE(IOX,IOY,IOZ) - CUBE(IOX+1,IOY,IOZ) 
     &   - CUBE(IOX,IOY,IOZ+1)+ CUBE(IOX+1,IOY,IOZ+1)
         A7 = CUBE(IOX,IOY,IOZ) - CUBE(IOX,IOY+1,IOZ) 
     &   - CUBE(IOX,IOY,IOZ+1)+ CUBE(IOX,IOY+1,IOZ+1)
         A8 = CUBE(IOX+1,IOY,IOZ) + CUBE(IOX,IOY+1,IOZ)
     &   + CUBE(IOX,IOY,IOZ+1)
     &   - CUBE(IOX,IOY,IOZ)- CUBE(IOX+1,IOY+1,IOZ) 
     &   - CUBE(IOX+1,IOY,IOZ+1)
     &   - CUBE(IOX,IOY+1,IOZ+1) + CUBE(IOX+1,IOY+1,IOZ+1)
         B(J,K)=B(J,K)+
     &    A1 + DZ*(A4 + A6*DX + (A7 + A8*DX)*DY) + A3*DY
     &   + DX*(A2 + A5*DY)
#endif
