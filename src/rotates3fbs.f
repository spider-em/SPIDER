C++*********************************************************************
C                                                                      *
C ROTATES3FBS.F   NEW                         Oct  4 2011 ArDean Leith *
C                 OMP                         Oct 31 2011 ArDean Leith *
C                 KLY,KNY BUG                 APR 27 2012 ArDean Leith *                                                                     *
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
C                                                                      *
C  ROTATES3FBS(LUN2,Q1,KLX,KNX,KLY,KNY,KLZ,KNZ,DRM,USEBACK,BACK)       *
C                                                                      *
C  PURPOSE:     3D ROTATION USING MATRIX DRM.                          *
C               AND FOUIER-BASED TRICUBIC SPLINE INTERPOLATION.        *
C               VOLUME TRUNCATED AT ORIGINAL BORDERS.                  *
C                                                                      *
C               SUBROUTINES FBS3_PREP and FUNCTION FBS3 are used       *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

       SUBROUTINE ROTATES3FBS(LUN2,Q1 ,KLX,KNX,  KLY,KNY, KLZ,KNZ,
     &                          DRM, USEBACK,BACK, IRTFLG)

       IMPLICIT NONE

       INTEGER            :: LUN2,KLX,KNX,KLY,KNY,KLZ,KNZ
       REAL               :: Q1(KLX:KNX, KLY:KNY, KLZ:KNZ)
       DOUBLE PRECISION   :: DRM(3,3)
       LOGICAL            :: USEBACK
       REAL               :: BACK
       INTEGER            :: IRTFLG

       INTEGER            :: I, J, K
       INTEGER            :: NXLD, NX, NY, NZ,ISLICE
       INTEGER            :: IZ, IY, IX
       INTEGER            :: IOX, IOY, IOZ

       REAL               :: Q2(KLX:KNX,KLY:KNY)  ! AUTOMATIC ARRAY
       REAL               :: QR(3)
       REAL               :: QR1,QR2,QR3

       REAL               :: FBS3

       REAL, ALLOCATABLE  :: F1  (:,:,:)

       REAL, ALLOCATABLE  :: XYZ (:,:,:)
       REAL, ALLOCATABLE  :: X1  (:,:,:)
       REAL, ALLOCATABLE  :: Y1  (:,:,:)
       REAL, ALLOCATABLE  :: Z1  (:,:,:)
       REAL, ALLOCATABLE  :: XY2 (:,:,:)
       REAL, ALLOCATABLE  :: XZ2 (:,:,:)
       REAL, ALLOCATABLE  :: YZ2 (:,:,:)

       NX   = KNX - KLX + 1
       NY   = KNY - KLY + 1
       NZ   = KNZ - KLZ + 1

       NXLD = NX + 2 - MOD(NX,2)

       ALLOCATE (F1  (NXLD, NY, NZ),
     &           XYZ (NXLD, NY, NZ),
     &           X1  (NXLD, NY, NZ),
     &           Y1  (NXLD, NY, NZ),
     &           Z1  (NXLD, NY, NZ),
     &           XY2 (NXLD, NY, NZ),
     &           XZ2 (NXLD, NY, NZ),
     &           YZ2 (NXLD, NY, NZ),
     &                  STAT=IRTFLG)

        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'ROTATES3FBS; XYZ...',7*NXLD*NY*NZ)
           RETURN
        ENDIF

C       PAD: Q1 INTO: F1 WITH PADDING
        DO K = 1,NZ
           DO J = 1,NY
              DO I = NX + 1,NXLD
                 F1(I,J,K) = 0
              ENDDO
           ENDDO
        ENDDO

        DO K = 1,NZ
           DO J = 1,NY
              DO I = 1,NX
                 F1(I,J,K) = Q1(I+KLX-1,J+KLY-1,K+KLZ-1)
              ENDDO
           ENDDO
        ENDDO

        XYZ = F1     ! ARRAY ASSIGNMENT

        CALL FBS3_PREP(XYZ, NXLD, NX, NY, NZ,
     &                  X1, Y1, Z1, XY2, XZ2, YZ2)

        ISLICE = 0

        DO IZ=KLZ,KNZ

C          OMP GAVE TIME: 56 ON 8 PROCS VS 80 ON 1
c$omp      parallel do private(iy,qr,ix,iox,ioy,ioz,qr1,qr2,qr3)
           DO IY=KLY,KNY

              QR(1) = DRM(1,1)*KLX + DRM(2,1)*IY + DRM(3,1)*IZ
              QR(2) = DRM(1,2)*KLX + DRM(2,2)*IY + DRM(3,2)*IZ
              QR(3) = DRM(1,3)*KLX + DRM(2,3)*IY + DRM(3,3)*IZ

              DO IX=KLX,KNX

C                IOX..  INTEGER LOCATION IN -NSAM/2...NSAM/2 ARRAY
                 IOX = FLOOR(QR(1))   ! QR CHANGES IN LOOP!
                 IOY = FLOOR(QR(2))   
                 IOZ = FLOOR(QR(3))   

                 IF ((IOX > KLX .AND. IOX < KNX) .AND.
     &               (IOY > KLY .AND. IOY < KNY) .AND.
     &               (IOZ > KLZ .AND. IOZ < KNZ)) THEN
C                  ROTATED POSITION IS INSIDE OF VOLUME

C                  EVALUATE INTENSITY AT: QR(..)

                   QR1 = QR(1) - KLX + 1
                   QR2 = QR(2) - KLY + 1
                   QR3 = QR(3) - KLZ + 1

                   Q2(IX,IY) = FBS3(QR1,QR2,QR3,
     &                            NXLD, NX, NY, NZ,
     &                            F1,NXLD, XYZ,
     &                            X1, Y1, Z1,
     &                            XY2,XZ2,YZ2)
                 ELSE
C                   ROTATED POSITION IS OUTSIDE VOLUME
                    IF (USEBACK) THEN
                       Q2(IX,IY) = BACK  
                    ELSE
                       Q2(IX,IY) = Q1(IX,IY,IZ)
                    ENDIF
                 ENDIF

                 QR(1) = QR(1) + DRM(1,1)
                 QR(2) = QR(2) + DRM(1,2)
                 QR(3) = QR(3) + DRM(1,3)
              ENDDO

           ENDDO

           ISLICE = ISLICE + 1
           CALL WRTVOL(LUN2,NX,NY,ISLICE,ISLICE,Q2,IRTFLG)
        ENDDO

        IF (ALLOCATED(F1))    DEALLOCATE(F1)
        IF (ALLOCATED(XYZ))   DEALLOCATE(XYZ)
        IF (ALLOCATED(X1))    DEALLOCATE(X1)
        IF (ALLOCATED(Y1))    DEALLOCATE(Y1)
        IF (ALLOCATED(Z1))    DEALLOCATE(Z1)
        IF (ALLOCATED(XY2))   DEALLOCATE(XY2)
        IF (ALLOCATED(XZ2))   DEALLOCATE(XZ2)
        IF (ALLOCATED(YZ2))   DEALLOCATE(YZ2)

        IRTFLG = 0

        END
