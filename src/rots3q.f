C++*********************************************************************
C
C ROTS3Q.F                                 NEW APRIL 4 2002 ArDean Leith
C
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
C
C  ROTS3Q(LUN2,Q1,KLX,KNX,KLY,KNY,KLZ,KNZ,PSI,THETA,PHI)
C 
C  PURPOSE:        3D ROTATION USING EULER ANGLES OF VOLUME. 
C                  TRI-QUADRATIC INTERPOLATION
C                  VOLUME TRUNCATED AT ORIGINAL BORDERS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE ROTS3Q(LUN2,Q1,KLX,KNX,KLY,KNY,KLZ,KNZ,
     &                     PSI,THETA,PHI)

         DIMENSION  Q1(KLX:KNX,KLY:KNY,KLZ:KNZ), Q2(KLX:KNX)
         DOUBLE     PRECISION  QR(3),RM(3,3),DX,DY,DZ

         INTEGER, PARAMETER :: NSIZE = 27

         INTEGER, DIMENSION(27) :: X,Y,Z
         REAL, DIMENSION(27) :: F

C        SET THE KNOWN COORDINATE GRID
C  Replaced by loops below, data does not agree with openmp.
c         DATA X/  
c     &          -1, 0, 1, -1, 0, 1, -1, 0, 1, 
c     &          -1, 0, 1, -1, 0, 1, -1, 0, 1, 
c     &          -1, 0, 1, -1, 0, 1, -1, 0, 1/ 

c         DATA Y/ 
c     &          -1,-1,-1,  0, 0, 0,  1, 1, 1, 
c     &          -1,-1,-1,  0, 0, 0,  1, 1, 1, 
c     &          -1,-1,-1,  0, 0, 0,  1, 1, 1/ 
 
c         DATA Z/  
c     &          -1,-1,-1, -1,-1,-1, -1,-1,-1, 
c     &           0, 0, 0,  0, 0, 0,  0, 0, 0,
c     &           1, 1, 1,  1, 1, 1,  1, 1, 1/
 

C        SET THE KNOWN COORDINATE GRID
	 DO  L=1,NSIZE,3
	  X(L)=-1
	  X(L+1)=0
	  X(L+2)=1
	  Y(L)=MOD(L/3,3)-1
	  Y(L+1)=MOD(L/3,3)-1
	  Y(L+2)=MOD(L/3,3)-1
	 ENDDO
	 DO  L=1,NSIZE
	  Z(L)=(L-1)/9-1
	 ENDDO
C
         CALL BLDR(RM,PSI,THETA,PHI)
         NSAM  = KNX - KLX + 1
         IREC  = 0

         DO IZ=KLZ,KNZ
           DO IY=KLY,KNY

              QR(1) = RM(1,1)*KLX+RM(2,1)*IY+RM(3,1)*IZ
              QR(2) = RM(1,2)*KLX+RM(2,2)*IY+RM(3,2)*IZ
              QR(3) = RM(1,3)*KLX+RM(2,3)*IY+RM(3,3)*IZ

              DO IX=KLX,KNX

C                IOX..  INTEGER LOCATION IN -NSAM/2...NSAM/2 ARRAY
                 IOX = FLOOR(QR(1))   
                 IOY = FLOOR(QR(2))   
                 IOZ = FLOOR(QR(3))   

C                DX.. OFFSET FROM INTEGER ARRAY
                 DX  = QR(1) - IOX
                 DY  = QR(2) - IOY
                 DZ  = QR(3) - IOZ

                 IF ((IOX.GT.KLX .AND. IOX.LT.(KNX)) .AND.
     &               (IOY.GT.KLY .AND. IOY.LT.(KNY)) .AND.
     &               (IOZ.GT.KLZ .AND. IOZ.LT.(KNZ))) THEN
C                   ROTATED POSITION IS INSIDE OF VOLUME

C                   FIND INTENSITIES ON 3x3x3 COORDINATE GRID
                    DO L = 1,NSIZE
                       I    = IOX + X(L)
                       J    = IOY + Y(L)
                       K    = IOZ + Z(L)
                       F(L) = Q1(I,J,K)
                    ENDDO

c             write(6,*) 'iox,ioy,ioz:',iox,ioy,ioz
c             write(6,*) 'qr:',qr
c             write(6,*) 'dx,dy,dz:',dx,dy,dz
c            write(6,*) 'dx,dy,dz:',dx,dy,dz
c             write(6,*) 'iox,ioy,ioz:',iox,ioy,ioz
c             stop

C                   EVALUATE INTENSITY AT PX,PY,PZ
                    Q2(IX) = TRIQUAD(DX,DY,DZ,F)

                 ELSE
C                   ROTATED POSITION IS OUTSIDE VOLUME
                    Q2(IX) = Q1(IX,IY,IZ)
                 ENDIF

                 QR(1) = QR(1) + RM(1,1)
                 QR(2) = QR(2) + RM(1,2)
                 QR(3) = QR(3) + RM(1,3)
              ENDDO

              IREC = IREC + 1
              CALL WRTLIN(LUN2,Q2,NSAM,IREC)
           ENDDO
        ENDDO
        END

