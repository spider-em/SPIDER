
C++*********************************************************************
C
C    NORM3.F
C           SIG FOR BLANK IMAGE FLOAT SLOP         JAN 06 ArDean Leith
C           FLOAT                                  APR 07 ArDean Leith
C           ERROR MESSAGES                         NOV 10 ArDean Leith
C           NORMVALSP ADDED                        AUG 11 ArDean Leith
C           NORMVALSP ADDED                        AUG 11 ArDean Leith
C           USED SETPRMB, IMPLICIT                 JUL 19 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C    NORM3(LUN,NX,NY,NZ,FMAX,FMIN,AV)
C
C    PURPOSE:  DETERMINES MAX, MIN, AV, AND SIG FOR A SPIDER IMAGE
C
C    PARAMETERS:
C        LUN          LOGICAL UNIT NUMBER OF IMAGE             (SENT)
C        NX,NY,NZ     DIMENSIONS OF VOLUME                     (SENT)
C        FMAX         MAXIMUM OF VOLUME                        (RET.)
C        FMIN         MINIMUM OF VOLUME                        (RET.)
C        AV           AVERAGE OF VOLUME                        (RET.)
C
C    NOTE:    THIS USES UNLABELED COMMON!!
C             SIG RETURNED IN COMMON  (LEGACY ISSUE)
C
C--*********************************************************************

      SUBROUTINE NORM3(LUN,NX,NY,NZ,FMAX,FMIN,AV)

      IMPLICIT NONE

      INTEGER          :: LUN,NX,NY,NZ
      REAL             :: FMAX,FMIN,AV

      REAL             :: BUF        
      COMMON BUF(1)

      INTEGER          :: NSAMC,NROWC,IREC,NLABEL,IFORM,IMAMI,IHIST
      REAL             :: FMAXC,FMINC,AVC,SIG
      COMMON /MASTER/NSAMC,NROWC,IREC,NLABEL,IFORM,IMAMI,FMAXC,FMINC,
     &                    AVC,SIG,IHIST

      INTEGER          :: LUNT,NIN,NOUT
      COMMON /UNITS/LUNT,NIN,NOUT

      DOUBLE PRECISION :: DAV,DAV2,DTOP,FNALL,DTEMP
      REAL             :: B,DIFF 
      INTEGER          :: NE,IRECT,K

      FNALL = FLOAT(NX) * FLOAT(NY) * FLOAT(NZ)
      IF (FNALL <= 0) THEN
         WRITE(NOUT,*) '*** NX, NY, NZ:',NX,NY,NZ
         CALL ERRT(101,'NORM3, BAD DIMENSIONS',NE)
         RETURN
      ENDIF

      DAV  = 0.0
      DAV2 = 0.0

      FMIN = HUGE(FMIN)
      FMAX = -FMIN

      DO IRECT = 1,NY*NZ
         CALL REDLIN(LUN,BUF,NX,IRECT)
         DO K = 1,NX
            B    = BUF(K)
            FMAX = MAX(B,FMAX)
            FMIN = MIN(B,FMIN)
            DAV  = DAV  + B
            DAV2 = DAV2 + B * DBLE(B)
          ENDDO
      ENDDO

      AV    = DAV / FNALL
      AVC   = AV
      FMAXC = FMAX
      FMINC = FMIN

      DTOP  = DAV2 - DAV * DAV / FNALL

      DIFF  = FMAX - FMIN
      IF (DIFF <= TINY(DIFF)) THEN
C        BLANK IMAGE SOMETIMES LEADS TO SQRT NEG. NUMBER

      ELSEIF (DTOP < 0.0D0) THEN
C        SQRT OF NEGATIVE NUMBER
         WRITE(NOUT,*) '*** SQRT(',DTOP,') IMPOSSIBLE. ',
     &                 'ASSUMING THIS IS A BLANK IMAGE' 
         SIG = 0.0
 
      ELSEIF (FNALL == 1.0) THEN
C        DIVISION BY ZERO
         CALL ERRT(101,'SINGLE PIXEL --> NO STANDARD DEVIATION',NE) 
         RETURN

      ELSE
C        CAN CALCULATE SIG
         SIG = DSQRT( DTOP / DBLE(FNALL - 1.0))
      ENDIF


      !write(6,*)' In norm3, before setprmb, fmin,fmax: ',fmin,fmax

C     VALUES ARE SET IN COMMON BY NAME ABOVE

C     WRITE(3,90) FMIN,FMAX,AV,SIG
90    FORMAT('  FMIN: ', 1PG10.3,
     &       '  FMAX: ', 1PG10.3,
     &       '  AV: ',   1PG12.5,
     &       '  SIG: ',  1PG12.5)

C     SET VALUES IN FILE HEADER
      CALL SETPRMB(LUN,FMAX,FMIN,AV,SIG)

      END


C****************************** NORMVALS ******************************
C
C NORMVALS.F         NEW                         FEB 2011 ARDEAN LEITH
C
C++*********************************************************************
C
C NORMVALS(X,NS,NR,DAV,DSIG,USE_OMP)
C
C PURPOSE: DETERMINE NORMALIZATION PARAMETERS: AVG & VARIANCE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE NORMVALS(X, NX, NY, NZ,
     &                         DAV,DSIG, USE_OMP)

        IMPLICIT NONE

        REAL, INTENT(IN)              :: X(NX,NY,NZ) 
        INTEGER, INTENT(IN)           :: NX, NY, NZ    ! SIZE 
        DOUBLE PRECISION, INTENT(OUT) :: DAV,DSIG 
        LOGICAL, INTENT(IN)           :: USE_OMP 

        DOUBLE PRECISION              :: DN,DVR,DTEMP,DSUM
        INTEGER                       :: K,J,I       

        DSUM  = 0.0
        DVR   = 0.0
        DN    = FLOAT(NX*NY*NZ)

        IF (USE_OMP) THEN
c$omp      parallel do private(k,j,i),reduction(+:dsum,dvr)
           DO K=1,NZ
              DO J=1,NY
                 DO I=1,NX
                   DSUM = DSUM + X(I,J,K)
                   DVR  = DVR  + X(I,J,K) * DBLE(X(I,J,K))
                 ENDDO
              ENDDO
           ENDDO
c$omp      end parallel do

        ELSE
           DO K=1,NZ
              DO J=1,NY
                 DO I=1,NX
                   DSUM = DSUM + X(I,J,K)
                   DVR  = DVR  + X(I,J,K) * DBLE(X(I,J,K))
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

        DAV   = DSUM / DN

        DTEMP = (DVR - DN * DAV * DAV)

        IF (DTEMP .GT. 0) THEN
           DSIG   = DSQRT(DTEMP / (DN-1))
        ELSE
C          TRAP FOR BLANK IMAGE AREA IMPRECISION
           DSIG = 0
        ENDIF

        END



C***********************************************************************
C
C NORMVALS_LMASKED.F         NEW                  SEP 2012 ARDEAN LEITH
C
C***********************************************************************
C
C NORMVALS_LMASKED()
C
C PURPOSE: DETERMINE AVG & STD. DEV NORMALIZATION PARAMETERS, 
C          INSIDE / OUTSIDE OF MASK 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE NORMVALS_LMASKED(X, LMASK,NX,NY,NZ, 
     &                            USE_OMP,
     &                            NI,DAVI,DSIGI, NO,DAVO,DSIGO)

        IMPLICIT NONE

        REAL, INTENT(IN)              :: X(NX,NY,NZ) 
        LOGICAL, INTENT(IN)           :: LMASK(NX,NY,NZ) 
        INTEGER, INTENT(IN)           :: NX, NY, NZ    ! SIZE 
        LOGICAL, INTENT(IN)           :: USE_OMP 
        INTEGER, INTENT(OUT)          :: NI,NO      
        DOUBLE PRECISION, INTENT(OUT) :: DAVI,DSIGI, DAVO,DSIGO 

        DOUBLE PRECISION              :: DNI,DNO,DTEMP
        DOUBLE PRECISION              :: DVRI,DSUMI
        DOUBLE PRECISION              :: DVRO,DSUMO
        INTEGER                       :: IX,IY,IZ      

        DSUMI  = 0.0
        DVRI   = 0.0
        DSUMO  = 0.0
        DVRO   = 0.0

        NI     = 0
        NO     = 0

        IF (USE_OMP) THEN
c$omp      parallel do private(iz,iy,ix),
c$omp&                      reduction(+:ni,dsumi,dvri,no,dsumo,dvro)
           DO IZ=1,NZ
              DO IY=1,NY
                 DO IX=1,NX

                   IF (LMASK(IX,IY,IZ)) THEN 
                      DSUMI = DSUMI +      X(IX,IY,IZ)
                      DVRI  = DVRI  + DBLE(X(IX,IY,IZ)**2)
                      NI    = NI + 1
                   ELSE 
                      DSUMO = DSUMO +      X(IX,IY,IZ)
                      DVRO  = DVRO  + DBLE(X(IX,IY,IZ)**2)
                      NO    = NO + 1
                   ENDIF
                 ENDDO
              ENDDO
           ENDDO
c$omp      end parallel do

        ELSE
           DO IZ=1,NZ
              DO IY=1,NY
                 DO IX=1,NX

                   IF (LMASK(IX,IY,IZ)) THEN 
                      DSUMI = DSUMI +      X(IX,IY,IZ)
                      DVRI  = DVRI  + DBLE(X(IX,IY,IZ)**2)
                      NI    = NI + 1
                   ELSE 
                      DSUMO = DSUMO +      X(IX,IY,IZ)
                      DVRO  = DVRO  + DBLE(X(IX,IY,IZ)**2)
                      NO    = NO + 1
                   ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

        DNI   = FLOAT(NI)
        DAVI  = DSUMI / DNI
        DTEMP = (DVRI - DNI * (DAVI**2) )

C       TRAP FOR BLANK IMAGE AREA 
        DSIGI = 0
        IF (DTEMP > 0) DSIGI = DSQRT(DTEMP / (DNI-1))

        DNO   = FLOAT(NO)
        DAVO  = DSUMO / DNO
        DTEMP = (DVRO - DNO * (DAVO**2) )

C       TRAP FOR BLANK IMAGE AREA 
        DSIGO = 0
        IF (DTEMP > 0) DSIGO = DSQRT(DTEMP / (DNO-1))

        END

C***********************************************************************
C
C NORMVALS_CMASKED.F         NEW                  SEP 2012 ARDEAN LEITH
C
C****************************** NORMVALS *******************************
C
C NORMVALS_CMASKED()
C
C PURPOSE: DETERMINE AVG & STD. DEV NORMALIZATION PARAMETERS, 
C          INSIDE / OUTSIDE OF CIRCULAR MASK OF GIVEN RADIUS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE NORMVALS_CMASKED(X, NX,NY,NZ, IXCEN,IYCEN,IZCEN,
     &                      IRADI, USE_OMP,
     &                      NI,DAVI,DSIGI, NO,DAVO,DSIGO)

        IMPLICIT NONE

        REAL, INTENT(IN)              :: X(NX,NY,NZ) 
        INTEGER, INTENT(IN)           :: NX, NY, NZ    ! SIZE 
        INTEGER, INTENT(IN)           :: IXCEN,IYCEN,IZCEN 
        INTEGER, INTENT(IN)           :: IRADI
        LOGICAL, INTENT(IN)           :: USE_OMP 
        INTEGER, INTENT(OUT)          :: NI,NO      
        DOUBLE PRECISION, INTENT(OUT) :: DAVI,DSIGI, DAVO,DSIGO 

        DOUBLE PRECISION              :: DNI,DNO,DTEMP
        DOUBLE PRECISION              :: DVRI,DSUMI
        DOUBLE PRECISION              :: DVRO,DSUMO
        INTEGER                       :: IX,IY,IZ      
        REAL                          :: RADISQ,RADZ,RADYZ,RADT 

        RADISQ = IRADI**2

        DSUMI  = 0.0
        DVRI   = 0.0
        DSUMO  = 0.0
        DVRO   = 0.0

        NI     = 0
        NO     = 0

        IF (USE_OMP) THEN
c$omp      parallel do private(iz,iy,ix,radz,radyz,radt),
c$omp&                      reduction(+:ni,dsumi,dvri,no,dsumo,dvro)
           DO IZ=1,NZ
              RADZ = (IZ - IZCEN) **2

              DO IY=1,NY
                 RADYZ = (IY - IYCEN) **2 + RADZ

                 DO IX=1,NX
                   RADT = RADYZ + (IX - IXCEN) **2 

                   IF (RADT <= RADISQ) THEN 
                      DSUMI = DSUMI +      X(IX,IY,IZ)
                      DVRI  = DVRI  + DBLE(X(IX,IY,IZ)**2)
                      NI    = NI + 1
                   ELSE 
                      DSUMO = DSUMO +      X(IX,IY,IZ)
                      DVRO  = DVRO  + DBLE(X(IX,IY,IZ)**2)
                      NO    = NO + 1
                   ENDIF
                 ENDDO
              ENDDO
           ENDDO
c$omp      end parallel do

        ELSE
           DO IZ=1,NZ
              RADZ = (IZ - IZCEN) **2

              DO IY=1,NY
                 RADYZ = (IY - IYCEN) **2 + RADZ

                 DO IX=1,NX
                   RADT = RADYZ + (IX - IXCEN) **2 

                   IF (RADT <= RADISQ) THEN 
                      DSUMI = DSUMI +      X(IX,IY,IZ)
                      DVRI  = DVRI  + DBLE(X(IX,IY,IZ)**2)
                      NI    = NI + 1
                   ELSE 
                      DSUMO = DSUMO +      X(IX,IY,IZ)
                      DVRO  = DVRO  + DBLE(X(IX,IY,IZ)**2)
                      NO    = NO + 1
                   ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

        DNI   = FLOAT(NI)
        DAVI  = DSUMI / DNI
        DTEMP = (DVRI - DNI * (DAVI**2) )

C       TRAP FOR BLANK IMAGE AREA 
        DSIGI = 0
        IF (DTEMP > 0) DSIGI = DSQRT(DTEMP / (DNI-1))

        DNO   = FLOAT(NO)
        DAVO  = DSUMO / DNO
        DTEMP = (DVRO - DNO * (DAVO**2) )

C       TRAP FOR BLANK IMAGE AREA 
        DSIGO = 0
        IF (DTEMP > 0) DSIGO = DSQRT(DTEMP / (DNO-1))

        END





C***********************************************************************
C
C NORMVALSP.F         NEW                  AUG 2011 ARDEAN LEITH
C
C***********************************************************************
C
C NORMVALSP()
C
C PURPOSE: FIND STATISTICS ON MASKED IMAGE/VOLUME
C
C--*********************************************************************

        SUBROUTINE NORMVALSP(X, NX, NY, NZ,
     &                          NXP,NYP,NZP, 
     &                          DAV,DSIG, USE_OMP)

C       FOR PADDED BUFFER

        IMPLICIT NONE

        REAL, INTENT(IN)              :: X(NXP,NYP,NZP) 
        INTEGER, INTENT(IN)           :: NX, NY, NZ     ! ANALYZED SIZE 
        INTEGER, INTENT(IN)           :: NXP,NYP,NZP    ! PADDED  SIZE
        DOUBLE PRECISION, INTENT(OUT) :: DAV,DSIG 
        LOGICAL, INTENT(IN)           :: USE_OMP        ! UNUSED 

        DOUBLE PRECISION              :: DN,DVR,DTEMP,DSUM
        INTEGER                       :: K,J,I       

        DSUM  = 0.0
        DVR   = 0.0
        DN    = FLOAT(NX*NY*NZ)

c$omp   parallel do private(k,j,i),reduction(+:dsum,dvr)
        DO K=1,NZ
           DO J=1,NY
              DO I=1,NX
                DSUM = DSUM + X(I,J,K)
                DVR  = DVR  + X(I,J,K) * DBLE(X(I,J,K))
              ENDDO
           ENDDO
        ENDDO
c$omp   end parallel do

        DAV   = DSUM / DN 

        DTEMP = (DVR - DN * DAV * DAV)

        IF (DTEMP > 0) THEN
           DSIG   = DSQRT(DTEMP / (DN-1))
        ELSE
C          TRAP FOR BLANK IMAGE AREA IMPRECISION
           DSIG = 0
        ENDIF

        END






C****************************** NORMVALSP_NOOMP ***********************
!! UNUSED!!

        SUBROUTINE NORMVALSP_NOOMP(X, NX, NY, NZ,
     &                          NXP,NYP,NZP, 
     &                          DAV,DSIG, USE_OMP)

C       COMPILER FAILED IF USE_OMP WAS PASSED, GAVE ODD COMPILATION

        IMPLICIT NONE

        REAL, INTENT(IN)              :: X(NXP,NYP,NZP) 
        INTEGER, INTENT(IN)           :: NX, NY, NZ     ! ANALYZED SIZE 
        INTEGER, INTENT(IN)           :: NXP,NYP,NZP    ! PADDED  SIZE
        DOUBLE PRECISION, INTENT(OUT) :: DAV,DSIG 
        LOGICAL, INTENT(IN)           :: USE_OMP        ! UNUSED 

        DOUBLE PRECISION              :: DN,DVR,DTEMP,DSUM
        INTEGER                       :: K,J,I       

        DSUM  = 0.0
        DVR   = 0.0
        DN    = FLOAT(NX*NY*NZ)

           DO K=1,NZ
              DO J=1,NY
                 DO I=1,NX
                   DSUM = DSUM + X(I,J,K)
                   DVR  = DVR  + X(I,J,K) * DBLE(X(I,J,K))
                 ENDDO
              ENDDO
           ENDDO

        DAV   = DSUM / DN 

        DTEMP = (DVR - DN * DAV * DAV)

        IF (DTEMP .GT. 0) THEN
           DSIG   = DSQRT(DTEMP / (DN-1))
        ELSE
C          TRAP FOR BLANK IMAGE AREA IMPRECISION
           DSIG = 0
        ENDIF

        END

