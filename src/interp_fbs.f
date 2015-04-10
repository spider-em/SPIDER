C ++********************************************************************
C                                                                      *
C  INTERP_FBS     NEW                     JUL 2011  GREGORY KISHCHENKO *                                                                           *
C                 OMP,FBS2_PREP           OCT 2011  ArDean Leith       *
C                 FBS2                    DEC 2011  GREGORY KISHCHENKO *
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
C   INTERP_FBS(BUF1, BUF2, NXLD,NX,  NY, NX2, NY2, IRTFLG)             *
C                                                                      *
C   PURPOSE: RESAMPLING OF 2D IMAGES BY FOURIER-BASED BICUBIC SPLINE   *
C            INTERPOLATION BETWEEN PIXELS.                             *
C            ALGORITHM IS FAIRLY FAST AND PRESERVES FINE DETAILS       *
C            OF  IMAGES                                                *
C                                                                      *
C            SUBROUTINES FBS2_PREP and FUNCTION FBS2 are used          *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE INTERP_FBS(BUF1, BUF2, 
     &                       NXLD,NX,  NY,
     &                            NX2, NY2, IRTFLG)

       IMPLICIT NONE

       INCLUDE 'CMBLOCK.INC'

       REAL              :: BUF1(NXLD,NY)
       REAL              :: BUF2(NX2,NY2)
       INTEGER           :: NXLD, NX, NY
       INTEGER           :: NX2, NY2
       INTEGER           :: IRTFLG

       INTEGER           :: I,J,K1,K2,MWANT
       INTEGER           :: M1,M2
       INTEGER           :: INV
       REAL              :: X,Y
       REAL              :: SCALEX, SCALEY

       REAL, ALLOCATABLE :: F0(:,:)
       REAL, ALLOCATABLE :: X1(:,:)
       REAL, ALLOCATABLE :: Y1(:,:)
       REAL, ALLOCATABLE :: XY2(:,:)

       REAL              :: fbs2

       WRITE(NOUT,*) ' Fourier based spline 2D interpolation'

       ALLOCATE (F0 (NXLD,NY),
     &           X1 (NXLD,NY),
     &           Y1 (NXLD,NY),
     &           XY2(NXLD,NY),
     &           STAT=IRTFLG)

       IF (IRTFLG .NE. 0) THEN 
          MWANT = 4* NXLD*NY 
          CALL ERRT(46,'INTERP_FBS, F0...',MWANT)
          RETURN
       ENDIF

       SCALEX = FLOAT(NX) / FLOAT(NX2)
       SCALEY = FLOAT(NY) / FLOAT(NY2)

       !WRITE(6,*) 'NX, NXLD, NX2= ',NX, NXLD, NX2
       !WRITE(6,*) 'NY, NY2=       ',NY, NY2
       !WRITE(6,*) 'SCALE X & Y=   ',SCALEX, SCALEY

       F0  = BUF1  ! ARRAY ASSIGNMENT

       CALL FBS2_PREP(F0, X1,Y1, XY2, NXLD, NX,NY, IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999

C      Timing 1400x1400 To: 500x500 
C     'IP   ' OMP GAVE TIME:  <1 ON 8 PROCS VS  <1 ON 1
C     'IP SF' OMP GAVE TIME:   3 ON 8 PROCS VS   4 ON 1

c$omp  parallel do private(k2,y,m2, k1,x,m1)
       DO K2 = 0,NY2-1
          Y  = K2 * SCALEY + 1
          M2 = K2 + 1

          DO K1 = 0,NX2-1
             X  = K1 * SCALEX + 1
             M1 = K1 + 1

             BUF2(M1,M2) = FBS2(X,Y, NXLD,NX,NY, BUF1,NXLD,
     &                          X1,Y1,XY2, .TRUE.)
          ENDDO
       ENDDO

       IRTFLG = 0

9999   IF (ALLOCATED(F0))  DEALLOCATE(F0)
       IF (ALLOCATED(X1))  DEALLOCATE(X1)
       IF (ALLOCATED(Y1))  DEALLOCATE(Y1)
       IF (ALLOCATED(XY2)) DEALLOCATE(XY2)

       END
