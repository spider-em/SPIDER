C **********************************************************************
C                                                                      *
C CENT_SYM       NEW                      FEB 2012  GREGORY KISHCHENKO *
C                SHIFTIT                  MAR 2012  ARDEAN LEITH       *                                                                   *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C                                                                      
C   SUBROUTINE CENT_SYM(BUF1,BUF2, SHIFTIT, NXLD,NX,NY, 
C                          SHX,SHY,IRTFLG)     
C                                                                      
C   PARAMETERS: BUF1    DATA ARRAY                              INPUT  
C               BUF2    DATA ARRAY                              OUTPUT 
C               SHX,SHY SHIFTS                                  OUTPUT 
C               SHIFTIT SHIFT FLAG                              INPUT 
C               IRTFLG  ERROR FLAG                              OUTPUT 
C                                                                      
C   PURPOSE: SHIFT OBJECT IN ORDER TO MOVE IT TO CENTER OF IMAGE       
C            USES 2D FOURIER-BASED BICUBIC SPLINE INTERPOLATION        
C            BETWEEN PIXELS.                                           
C                                                                      
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE CENT_SYM(BUF1,BUF2, SHIFTIT,
     &                     NXLD,NX,NY,
     &                     SHX,SHY, IRTFLG)

       IMPLICIT NONE

       REAL              :: BUF1(NXLD,NY)
       REAL              :: BUF2(NX,  NY)
       INTEGER           :: NXLD, NX, NY
       REAL              :: SHX, SHY
       LOGICAL           :: SHIFTIT
       INTEGER           :: IRTFLG

       INTEGER           :: INV, MWANT
       INTEGER           :: CX, CY
       INTEGER           :: I,J
       INTEGER           :: NX25,NX75,NY25,NY75

       REAL              :: X,Y
       REAL              :: SHX2, SHY2
       REAL              :: COR
       REAL              :: WEIGHT
       REAL              :: A,B,C,D

       REAL              :: fbs2

       REAL, ALLOCATABLE :: F0(:,:)
       REAL, ALLOCATABLE :: F1(:,:)
       REAL, ALLOCATABLE :: X1(:,:)
       REAL, ALLOCATABLE :: Y1(:,:)
       REAL, ALLOCATABLE :: XY2(:,:)

       ALLOCATE (F0 (0:NXLD-1,0:NY-1),
     &           F1 (0:NXLD-1,0:NY-1),
     &           X1 (0:NXLD-1,0:NY-1),
     &           Y1 (0:NXLD-1,0:NY-1),
     &           XY2(0:NXLD-1,0:NY-1),
     &           STAT=IRTFLG)

       IF (IRTFLG .NE. 0) THEN
          MWANT = 4* NXLD*NY
          CALL ERRT(46,'CENT_SYM, F0...',MWANT)
          RETURN
       ENDIF

       F0 = BUF1  ! ARRAY ASSIGN

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

       F0  = BUF1  ! ARRAY ASSIGNMENT

       CX  = NX / 2 + 1
       CY  = NY / 2 + 1

C      GENERATION OF IMAGE WITH OBJECT BY OPTIMAL SHIFT
c$omp  parallel do private(j,i,x,y)
       DO J=1, NY
          DO I=1, NX
            X = I  + SHX
            Y = J  + SHY
            BUF2(I,J) = FBS2(X,Y, NXLD,NX,NY,F0,
     &	                     NXLD,X1,Y1,XY2, .TRUE.)
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

       IF (SHIFTIT) THEN
C         GENERATION OF IMAGE WITH OBJECT BY CORRECTED OPTIMAL SHIFT

c$omp     parallel do private(j,i,x,y)
          DO J=1,NY 
             DO I=1,NX 
                X         = I + SHX2
                Y         = J + SHY2
                BUF2(I,J) = FBS2(X,Y, NXLD,NX,NY,F0,
     &	   	                 NXLD,X1,Y1,XY2, .TRUE.)
             ENDDO
          ENDDO
       ENDIF

       IF (ALLOCATED(F0))  DEALLOCATE(F0)
       IF (ALLOCATED(F1))  DEALLOCATE(F1)
       IF (ALLOCATED(X1))  DEALLOCATE(X1)
       IF (ALLOCATED(Y1))  DEALLOCATE(Y1)
       IF (ALLOCATED(XY2)) DEALLOCATE(XY2)

C       write(6,*), 'A,B,C,D =',A,B,C,D
C       write(6,*), 'SHX, SHY = ', SHX2, SHY2

       END
