C ++********************************************************************
C                                                                    
C SPHDECON     NEW                          FEB 12  GREGORY KISHCHENKO  
C                                                                    
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
C
C SPHDECON()
C
C PURPOSE:  SPHERICAL DECONVOLUTION OF VOLUME
C
C--*********************************************************************

        SUBROUTINE SPHDECON()

        IMPLICIT NONE

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FILNAM1,FILNAM2

        REAL, ALLOCATABLE     :: X1F(:,:,:), X2F(:,:,:)  
        REAL, ALLOCATABLE     :: X1F_90(:,:,:)  
        REAL, ALLOCATABLE     :: X2F_90(:,:,:)  

        REAL, PARAMETER       :: PI  = 3.14159265358979323846
        REAL                  :: PIRAD

        REAL                  :: SIGMA
        REAL                  :: X,  Y,  Z, RAD, TET

        INTEGER               :: MAXIM1,MAXIM2,IFORM1,IFORM2
        INTEGER               :: NX,NY,NZ
        INTEGER               :: IRTFLG,NOT_USED
        INTEGER               :: NXLD, NYLD, MWANT
        INTEGER               :: I, J, K
        INTEGER               :: CX, CY, CZ
        INTEGER               :: NRAD

	INTEGER, PARAMETER    :: LUN1   = 21
	INTEGER, PARAMETER    :: LUN2   = 22

C       OPEN INPUT FILE
        MAXIM1  = 0
        CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',IFORM1,NX,NY,NZ,
     &                MAXIM1,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL RDPRM1S(SIGMA,NOT_USED,'SIGMA',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C	GET OUTPUT FILE NAME
        IFORM2 = 3
        CALL OPFILEC(0,.TRUE.,FILNAM2,LUN2,'U',IFORM2,NX,NY,NZ,
     &                MAXIM2,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999


        CX     = NX / 2 + 1
        CY     = NY / 2 + 1
        CZ     = NZ / 2 + 1
        NXLD   = NX + 2 - MOD(NX,2)
        NYLD   = NY + 2 - MOD(NY,2)

        NRAD   = MIN(NX,NY,NZ) / 2 - 2

        ALLOCATE (X1F   (NXLD, NY,NZ),
     &            X2F   (NX,   NY,NZ),
     &            X1F_90(NYLD, NX,NZ),
     &            X2F_90(NY,   NX,NZ), STAT=IRTFLG)

        IF (IRTFLG .NE. 0) THEN 
            MWANT = NXLD*NY*NZ + 2*NX*NY*NZ + NYLD*NX*NZ
            CALL ERRT(46,'SPHDECON, X1F & X2F',MWANT)
            GOTO 9999
        ENDIF 

C       LOAD VOLUME WITH FOURIER PAD
        CALL READV(LUN1,X1F, NXLD,NY, NX,NY,NZ) ! PAD FOR FFT

        X1F_90 = 0

        DO K=1,NZ
           DO J=1,NY
              DO I=1,NX
                 X1F_90(J, I, K) = X1F(I, J, K)
              ENDDO
           ENDDO
        ENDDO

        CALL SPHDECON2(X1F,   X2F,   NXLD,NX,NY,NZ,SIGMA,IRTFLG)
        CALL SPHDECON2(X1F_90,X2F_90,NYLD,NY,NX,NZ,SIGMA,IRTFLG)

        PIRAD = 2 * PI / REAL(NRAD)

       DO K=1,NZ
          Z = (K-CZ)
          DO J=1,NY
             Y = (J-CY)
             DO I=1,NX
                X    = (I-CX)
                RAD  = SQRT(X**2 + Y**2 + Z**2)

                IF (RAD == 0) THEN
                   TET = 1.0
                ELSE
                   TET = ABS(X/RAD)
                ENDIF

                IF (TET >= 0.707) THEN
                         TET = (TET-0.707) / (1-0.707)
                X2F(I, J, K) = (1-TET) * X2F(I, J, K) +
     &                                TET * X2F_90(J,I,K)
                ENDIF

             ENDDO
          ENDDO
       ENDDO

C       X2F IS NOT FOURIER PADDED!
        CALL WRITEV(LUN2,X2F, NX,NY, NX,NY,NZ)

9999    CLOSE(LUN1)
        CLOSE(LUN2)

        IF (ALLOCATED(X1F))    DEALLOCATE(X1F)
        IF (ALLOCATED(X2F))    DEALLOCATE(X2F)
        IF (ALLOCATED(X1F_90)) DEALLOCATE(X1F_90)
        IF (ALLOCATED(X2F_90)) DEALLOCATE(X2F_90)

        END



C      --------------------- SPHDECON2 ------------------------------

       SUBROUTINE SPHDECON2(BUFIN,BUFOUT,
     &                      NXLD, NNX,NNY,NNZ, SIGMAT,IRTFLG )

       IMPLICIT NONE

       REAL              :: BUFIN (NXLD, NNY, NNZ)
       REAL              :: BUFOUT(NNX,   NNY, NNZ)
       INTEGER           :: NXLD, NNX, NNY, NNZ
       REAL              :: SIGMAT
       INTEGER           :: IRTFLG

       INTEGER           :: CX, CY, CZ
       INTEGER           :: NRAD, NRING, NRLD
       INTEGER           :: I, J, K, IRAD
       INTEGER           :: J2, K2, J3, K3
       INTEGER           :: IX,  IY,  IZ
       INTEGER           :: IX1, IY1, IZ1
       INTEGER           :: I2

       REAL              :: SIGMA, TEMP
       REAL              :: A1, A2, A3, A4 
       REAL              :: A5, A6, A61, A7, A8 
       REAL              :: X,  Y,  Z, DX, DY, DZ
       REAL              :: RAD, DELTA, FI, TET
       REAL              :: PIRING, RINGPI

       INTEGER           :: INV,NUMTH,MWANT

       REAL, ALLOCATABLE :: FLTR(:)
       REAL, ALLOCATABLE :: MER(:), FMER(:)
       REAL, ALLOCATABLE :: EQU(:,:,:), EQU_2(:,:,:), FEQU (:)

       LOGICAL, PARAMETER:: SPIDER_SIGN  = .FALSE.
       LOGICAL, PARAMETER:: SPIDER_SCALE = .FALSE.

       INTEGER, PARAMETER:: INV_F =  1
       INTEGER, PARAMETER:: INV_R = -1

       REAL, PARAMETER   :: PI  = 3.14159265358979323846
       REAL, PARAMETER   :: PI2 = 6.28318530717958647692

C      PLAN.. ARE POINTERS TO A STRUCTURE 
       INTEGER *8        :: PLAN_MERF = 0
       INTEGER *8        :: PLAN_MERR = 0
       INTEGER *8        :: PLAN_EQUF = 0
       INTEGER *8        :: PLAN_EQUR = 0

       NRAD   = MIN(NNX,NNY,NNZ) / 2 - 2

       X      = PI2 * NRAD
       X      = LOG(X) / LOG(2.0)
       NRING  = INT(X) + 1
       NRING  = 2**NRING

       NRLD   = NRING + 2 - MOD(NRING,2)
       PIRING = PI2 / REAL(NRING)
       RINGPI = REAL(NRING) / PI2

       CX     = NNX / 2 + 1
       CY     = NNY / 2 + 1
       CZ     = NNZ / 2 + 1

       SIGMA  = SIGMAT / 360
c       SIGMA = 0.5*(1.0/REAL(NRAD) + 0.0063739 * SIGMA)

       ALLOCATE(EQU ((NRING/2+1),NRING,0:NRAD+1),
     &          EQU_2(NRING,(NRING/2+1),0:NRAD+1), 
     &          FLTR(0:(NRLD/2 -1)),
     &          MER (NRING),
     &          FMER(NRLD),
     &          FEQU(NRLD),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
          MWANT = NRING * NRING * (NRING + 1) + NRING + 2*NRLD
          CALL ERRT(46,'SPHDECON2; EQU....',MWANT)
          RETURN
       ENDIF

       BUFOUT = 0  ! ARRAY ZERO
       MER    = 0  ! ARRAY ZERO
       FMER   = 0  ! ARRAY ZERO
       EQU    = 0  ! ARRAY ZERO
       EQU_2  = 0  ! ARRAY ZERO
       FEQU   = 0  ! ARRAY ZERO

       !write(6,'(A,I5,A,I5)') ' NRAD:',NRAD,'  NRING:',NRING

C      CREATE FFTW PLANS AND CACHE THEM
       NUMTH = 1    ! SINCE USED INSIDE PARALLEL REGION

       CALL FFTW3_MAKEPLAN(NRING,1,1,NUMTH, PLAN_MERF, 1,IRTFLG)
       CALL FFTW3_MAKEPLAN(NRING,1,1,NUMTH, PLAN_MERR,-1,IRTFLG)

       CALL FFTW3_MAKEPLAN(NRING,1,1,NUMTH, PLAN_EQUF, 1,IRTFLG)
       CALL FFTW3_MAKEPLAN(NRING,1,1,NUMTH, PLAN_EQUR,-1,IRTFLG)

C      ************* CALCULATE FOURIER FILTER *************
       TEMP = 0.65 / SIGMA
       I2   = INT(TEMP)

       IF (I2 >= (NRLD/2 -1)) THEN
          DO I=0,NRLD/2 -1
             FLTR(I) = SIGMA * I
             FLTR(I) = 1 - FLTR(I)
          ENDDO
       ELSE
          DO I=0,I2
             FLTR(I) = SIGMA * I
             FLTR(I) = 1 - FLTR(I)
          ENDDO
          DO I=I2+1,NRLD/2 -1
             FLTR(I) = SIGMA  * (I - TEMP)
             FLTR(I) = 0.35 * EXP(-FLTR(I )/ 0.35)
          ENDDO
       ENDIF

       DO I=0,NRLD/2 -1
          FLTR(I) = FLTR(I) / (0.95 * FLTR(I)**2 + 0.05)
       ENDDO

C      ************* SPHERICAL DECONVOLUTION *************

C      CYCLE ON RADIUS

       DO IRAD = 0, NRAD
           !write(6,*) 'IRAD  =',IRAD

C          CYCLE ON DELTA-INCLINATION (DIFFERENT "EQUATORS")
           DO K=1,NRING/2

             DELTA = (K-1) * PIRING

C            CYCLE ON FI
             DO J=1,NRING/2

               MER = 0      ! ARRAY(NRING) ZERO

C              CYCLE ON THETA

c$omp  parallel do private(i,tet,fi,x,y,z,
c$omp&                  ix,iy,iz,dx,dy,dz,
c$omp&                  a1,a2,a3,a4,a5,a61,a6,a7,a8)

             DO I=1,NRING

                  IF (I <= (NRING/4)) THEN
                     TET = 0.5*PI - REAL(I-1) * PIRING
                     FI  =          REAL(J-1) * PIRING
                  ENDIF

                  IF(I > (NRING/4) .AND. I <= (3*NRING/4)) THEN
                     TET = -0.5*PI + REAL(I-1) * PIRING
                     FI  = PI      + REAL(J-1) * PIRING
                  ENDIF

                  IF ((I) > (3*NRING/4)) THEN
                     TET = 2.5*PI - REAL(I-1) * PIRING
                     FI  =          REAL(J-1) * PIRING
                  ENDIF

                  X   = IRAD * SIN(TET) * COS(FI)  + CX
                  Y   = IRAD * (-SIN(DELTA) * COS(TET) +
     &                  COS(DELTA) * SIN(TET) * SIN(FI)) + CY
                  Z   = IRAD * (COS(DELTA)*COS(TET) +
     &                  SIN(DELTA)*SIN(TET)*SIN(FI)) + CZ

                  IX  = INT(X)
                  IY  = INT(Y)
                  IZ  = INT(Z)

                  DX  = X - IX
                  DY  = Y - IY
                  DZ  = Z - IZ

                  A1  =       BUFIN(IX,  IY,  IZ)
                  A2  =       BUFIN(IX+1,IY,  IZ)   - A1
                  A3  =       BUFIN(IX,  IY+1,IZ)   - A1
                  A4  =       BUFIN(IX,  IY,  IZ+1) - A1
                  A5  = -A2 - BUFIN(IX,  IY+1,IZ)   + 
     &                        BUFIN(IX+1,IY+1,IZ)
                  A61 = -     BUFIN(IX,  IY,  IZ+1) + 
     &                        BUFIN(IX+1,IY,  IZ+1)
                  A6  = -A2 + A61
                  A7  = -A3 - BUFIN(IX,  IY,  IZ+1) + 
     &                        BUFIN(IX,  IY+1,IZ+1)
                  A8  = -A5 - A61 -
     &                        BUFIN(IX,  IY+1,IZ+1) + 
     &                        BUFIN(IX+1,IY+1,IZ+1)

                  MER(I) = A1 + 
     &                   DZ*(A4 + A6*DX + 
     &                  (A7 + A8*DX)*DY) + A3*DY +
     &                   DX*(A2 + A5*DY)

               ENDDO ! END OF: DO I=1,NRING
c$omp end parallel do 

               DO I=1,NRING
                  FMER(I) = MER(I)
               ENDDO
               DO I=NRING+1,NRLD
                  FMER(I) = 0
               ENDDO

C              ************* Meridianal Deconvolution *************

C              1D FOURIER TRANSFORM OF MERIDIAN
               CALL FMRS(FMER, NRING,1,1, PLAN_MERF,
     &                  SPIDER_SIGN,SPIDER_SCALE, INV_F,IRTFLG)

               DO I=0,NRLD/2 -1
                  FMER(2*I+1) = FMER(2*I+1) * FLTR(I)
                  FMER(2*I+2) = FMER(2*I+2) * FLTR(I)
               ENDDO

C              REVERSE 1D FOURIER TRANSFORM OF MERIDIAN
               CALL FMRS(FMER, NRING,1,1, PLAN_MERR,
     &                  SPIDER_SIGN,SPIDER_SCALE, INV_R,IRTFLG)

               EQU(K,J,IRAD)  = FMER(1)

               J2             = J + NRING/2
               J2             = MOD(J2-1,NRING) +1
               EQU(K,J2,IRAD) = FMER(1+ NRING/2)

             ENDDO   ! END OF:  DO J=1,NRING/2 +1
           ENDDO     ! END OF:  DO K=1,NRING/2
        ENDDO        ! END OF:  DO IRAD = 0, NRAD

C       ************* Equatorial Deconvolution *************

        DO IRAD = 0, NRAD
           DO K=1,NRING/2

             DO J=1,NRING
                FEQU(J) = EQU(K,J,IRAD)
             ENDDO
             DO J=NRING+1,NRLD
                   FEQU(J) = 0
             ENDDO

             CALL FMRS(FEQU, NRING,1,1, PLAN_EQUF,
     &                  SPIDER_SIGN,SPIDER_SCALE, INV_F,IRTFLG)

               DO I=0,NRLD/2 - 1
                  FEQU(2*I+1) = FEQU(2*I+1) * FLTR(I)
                  FEQU(2*I+2) = FEQU(2*I+2) * FLTR(I)
               ENDDO

             CALL FMRS(FEQU, NRING,1,1, PLAN_EQUR,
     &                  SPIDER_SIGN,SPIDER_SCALE, INV_R,IRTFLG)

             DO J=1,NRING
                EQU(K,J,IRAD) =  FEQU(J)
             ENDDO

           ENDDO
        ENDDO

C       ****** Arrange Data in Regular Spherical Coordinates ******

        DO IRAD = 0, NRAD
           DO K=1,NRING/2

             DO J = 1, NRING/2 + 1
                EQU_2(K,J,IRAD) = EQU (K,J,IRAD)
             ENDDO

             DO J = NRING/2 +2, NRING
                K2                = MOD((K+NRING/2-1),NRING) + 1
                J2                = NRING - J + 2
                EQU_2(K2,J2,IRAD) = EQU (K,J,IRAD)
             ENDDO

           ENDDO
        ENDDO

C      ************* Back to Cartesian *************

       DO K=1,NNZ
          Z = (K-CZ)

          DO J=1,NNY
             Y = (J-CY)

             DO I=1,NNX
                X    = (I-CX)

                RAD  = SQRT(X**2 + Y**2 + Z**2)
                IRAD = INT(RAD)

             IF (IRAD > NRAD) THEN
                BUFOUT(I, J, K)= 0
             ELSE

                   FI  = ATAN2(Z,Y)*RINGPI

                   IF (Y == 0 .AND. Z == 0) THEN
                      TET = (NRING/4) * (1-(X+0.1) / ABS(X+0.1))
                   ELSE
                      TET = ACOS(X/RAD)*RINGPI
                   ENDIF

                   FI    = FI  + NRING
                   K2    = INT(FI)
                   DX    = FI - K2

                   K2    = K2 + 1
                   K2    = MOD(K2-1,NRING)+1
                   K3    = MOD(K2,NRING) + 1

                   J2    = INT(TET)
                   DY    = TET - J2
                   J2    = J2 + 1
                   J3    = J2 + 1

                   DZ    = RAD-IRAD

                   IF (J3 > NRING/2 + 1) THEN
                      J3 = NRING - J3 + 2
                      K3 = MOD(K3+NRING/2 -2,NRING) + 1
                   ENDIF

                   A1 = EQU_2(K2,J2,IRAD)
                   A2 = EQU_2(K3,J2,IRAD)   - A1
                   A3 = EQU_2(K2,J3,IRAD)   - A1
                   A4 = EQU_2(K2,J2,IRAD+1) - A1
                   A5 = -A2 - EQU_2(K2,J3,IRAD) +
     &                        EQU_2(K3,J3,IRAD)
                   A61= - EQU_2(K2,J2,IRAD+1) +
     &                    EQU_2(K3,J2,IRAD+1)
                   A6 = -A2 + A61
                   A7 = -A3 - EQU_2(K2,J2,IRAD+1) +
     &                        EQU_2(K2,J3,IRAD+1)
                   A8 = -A5 - A61 -
     &                  EQU_2(K2,J3,IRAD+1) +
     &                  EQU_2(K3,J3,IRAD+1)

                   BUFOUT(I, J, K) = A1 + 
     &                DZ*(A4 + A6*DX + 
     &                (A7 + A8*DX)*DY) + A3*DY +
     &                DX*(A2 + A5*DY)

                   BUFOUT(I, J, K)= BUFOUT(I, J, K)/(NRING**2)

                ENDIF
             ENDDO
          ENDDO
       ENDDO

       CALL FFTW3_KILLPLAN(PLAN_MERF,IRTFLG)
       CALL FFTW3_KILLPLAN(PLAN_MERR,IRTFLG)
       CALL FFTW3_KILLPLAN(PLAN_EQUF,IRTFLG)
       CALL FFTW3_KILLPLAN(PLAN_EQUR,IRTFLG)

9999   IF (ALLOCATED(MER))    DEALLOCATE(MER)
       IF (ALLOCATED(FMER))   DEALLOCATE(FMER)
       IF (ALLOCATED(EQU))    DEALLOCATE(EQU)
       IF (ALLOCATED(FEQU))   DEALLOCATE(FEQU)
       IF (ALLOCATED(EQU_2))  DEALLOCATE(EQU_2)
       IF (ALLOCATED(FLTR))   DEALLOCATE(FLTR)

       END

