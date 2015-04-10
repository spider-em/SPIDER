C++*********************************************************************
C
C  CROSRNG.F 
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
C  CROSRNG(CIRC1,CIRC2, LCIRC,NRING,T,Q,MAXRIN,
C            JACUP,NUMR,QN,TOT,MODE)
C
C PURPOSE: CROSS CORRELATION OF RADIAL RINGS FOR USE IN ROTATIONAL
C          ALIGNMENT. CHECKS ONLY UN-MIRRORED POSITION
C
C PARAMETERS:
C    CIRC1  - FT OF RINGS MULTIPLIED BY WEIGHTS           (SENT)
C    CIRC2  - FT OF RINGS MULTIPLIED BY WEIGHTS           (SENT)
C    LCIRC  - SIZE OF CIRCS ARRAYS                        (SENT)
C    NRING  - NUMBER OF RINGS                             (SENT)
C    MAXRIN - LONGEST RING                                (SENT)
C    NUMR   - RING LOCATION POINTERS                      (SENT)
C    Q      - CC ARRAY                                    (RETURNED)
C    QMAX   - CC MAX                                      (RETURNED)
C    POS_MAX - POSITION OF CC MAX                         (RETURNED)
C    TT     - USED FOR SGI FFT (UNUSED NOW)               (SENT)
C    NEG    - FLAG FOR CONJUGATE (MIRROR) OF 1'ST RING    (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE CROSRNG(CIRC1,CIRC2, LCIRC,NRING, T,Q, MAXRIN,
     &                     JACUP,NUMR,QN,TOT,MODE)

        INTEGER  NUMR(3,NRING),MAXRIN
        DIMENSION  CIRC1(LCIRC),CIRC2(LCIRC)
        DOUBLE PRECISION  T(MAXRIN),Q(MAXRIN)
        DOUBLE PRECISION  QN,QT,PI,T7(-3:3)
        CHARACTER*1  MODE

        PI = 4.0*DATAN(1.0D0)
        IF (MODE .EQ. 'F')  PI=2*PI
        IP = -LOG2(MAXRIN)

        Q = 0.0D0

        DO I=1,NRING
           WR = REAL(NUMR(1,I)) * PI/REAL(NUMR(3,I)) * 
     &          REAL(MAXRIN) / REAL(NUMR(3,I))

           T(1) = DBLE(CIRC1(NUMR(2,I))) * CIRC2(NUMR(2,I))
           IF (NUMR(3,I) .EQ. MAXRIN) THEN
              T(2) = DBLE(CIRC1(NUMR(2,I)+1)) * CIRC2(NUMR(2,I)+1)
              DO J=3,MAXRIN,2
                 JC = J + NUMR(2,I)-1
                 T(J)   =  DBLE(CIRC1(JC))   * CIRC2(JC)   + 
     &                     DBLE(CIRC1(JC+1)) * CIRC2(JC+1)
                 T(J+1) = -DBLE(CIRC1(JC))   * CIRC2(JC+1)  + 
     &                     DBLE(CIRC1(JC+1)) * CIRC2(JC)
              ENDDO

              Q = Q + T * WR
           ELSE
              T(NUMR(3,I)+1) = DBLE(CIRC1(NUMR(2,I)+1)) * 
     &                              CIRC2(NUMR(2,I)+1)/2.0
              T(2) = 0.0
              DO J=3,NUMR(3,I),2
                 JC=J+NUMR(2,I)-1
                 T(J)   =  DBLE(CIRC1(JC)) * CIRC2(JC)   + 
     &                     DBLE(CIRC1(JC+1))*CIRC2(JC+1)
                 T(J+1) = -DBLE(CIRC1(JC)) * CIRC2(JC+1) + 
     &                     DBLE(CIRC1(JC+1))*CIRC2(JC)
              ENDDO

              Q(1:NUMR(3,I)+1) = Q(1:NUMR(3,I)+1) + 
     &                           T(1:NUMR(3,I)+1) * WR
           ENDIF
        ENDDO

        CALL  FFTR_D(Q,IP)

        QN = -1.0D20
        DO J=1,MAXRIN
           IF (Q(J) .GE. QN)  THEN
              QN   = Q(J)
              JTOT = J
           ENDIF
        ENDDO

        IF (JACUP .EQ. 0) THEN
           TOT = JTOT
        ELSE
           DO K=-3,3
              J     = MOD(JTOT+K+MAXRIN-1,MAXRIN)+1
              T7(K) = Q(J)
           ENDDO
           CALL PRB1D(T7,7,POS)

           K = INT(100.0 / REAL(JACUP+1))
           TOT = REAL(JTOT) + REAL(IFIX(POS)) +
     &           REAL(K) * REAL(INT(POS*100.0 / REAL(K))) / 100.0
        ENDIF

        END


C       --------------- CROSRNG_NEW --------------------------------

        SUBROUTINE CROSRNG_NEW(CIRCR,CIRCE,LCIRC, 
     &                         NRING,MAXRIN,NUMR,
     &                         FFTW3PLAN, USE_MIR,
     &                         Q, QMAX,POS_MAX, MAXL)

C       USES NUMR TABLE FOR MAPPING INTO Q ARRAY 
C       USES SIMPLIFIED LOGIC FOR BOUNDARY VALUES, FLOATING PT. ARITH.

        IMPLICIT NONE

        REAL,            INTENT(IN)   :: CIRCR(LCIRC), CIRCE(LCIRC)
        INTEGER,         INTENT(IN)   :: LCIRC,NRING,MAXRIN
        INTEGER,         INTENT(IN)   :: NUMR(3,NRING)
        INTEGER *8                    :: FFTW3PLAN(*)
        REAL,            INTENT(OUT)  :: Q(MAXRIN+2)
        REAL,            INTENT(OUT)  :: QMAX
        INTEGER,         INTENT(OUT)  :: MAXL
        REAL,            INTENT(OUT)  :: POS_MAX
        LOGICAL,         INTENT(IN)   :: USE_MIR

        INTEGER                       :: I,IGO,IGOM1,NVAL,J1,J,JC

C       ZERO WHOLE Q ARRAY,  STRAIGHT  = CIRCR * CONJG(CIRCE)
        Q = 0.0
     
        DO I=1,NRING
           IGO   = NUMR(2,I)
           IGOM1 = IGO - 1
           NVAL  = NUMR(3,I)       
           J1    = 1

           IF (USE_MIR) THEN
C             FIRST RING SET IS CONJUGATED (MIRRORED)
              DO J=J1,NVAL,2
                 JC     = J + IGOM1

                 Q(J)   = Q(J)   + (CIRCR(JC))   * CIRCE(JC)   -
     &                             (CIRCR(JC+1)) * CIRCE(JC+1)
                 Q(J+1) = Q(J+1) - (CIRCR(JC))   * CIRCE(JC+1) -
     &                             (CIRCR(JC+1)) * CIRCE(JC)
              ENDDO
           ELSE
C             FIRST RING SET IS NON-CONJUGATED (NOT MIRRORED)
              DO J=J1,NVAL,2
                 JC     = J + IGOM1

                 Q(J)   = Q(J)   + (CIRCR(JC))   * CIRCE(JC)   +
     &                             (CIRCR(JC+1)) * CIRCE(JC+1)
                 Q(J+1) = Q(J+1) - (CIRCR(JC))   * CIRCE(JC+1) +
     &                             (CIRCR(JC+1)) * CIRCE(JC)
              ENDDO
           ENDIF
        ENDDO

        CALL CROSRNG_COM_R(Q,MAXRIN,FFTW3PLAN,QMAX,POS_MAX,MAXL)
 
        END

