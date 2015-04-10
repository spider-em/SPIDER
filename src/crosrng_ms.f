
C++*********************************************************************
C
C  CROSRNG_MS.F   
C              TRAP FOR COMPILER ERROR ON ALTIX   FEB 2005 ARDEAN LEITH
C              REWRITE USING CROSRNG_COM          FEB 2008 ARDEAN LEITH
C              REWRITE USING FFTW3 RINGS          FEB 2008 ARDEAN LEITH
C              REMOVED TT                         JUN 2010 ARDEAN LEITH
C              REMOVED UNUSED SUBROUTINES         OCT 2010 ARDEAN LEITH
C
C **********************************************************************
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
C  CROSRNG_MS(CIRC1,CIRC2,LCIRC,NRING,MAXRIN,NUMR,
C             QMAX,POS_MAX, QMAXM,POS_MAXM, TT)
C
C PURPOSE: CROSS CORRELATION OF RADIAL RINGS FOR USE IN ROTATIONAL
C          ALIGNMENT.  CHECKS BOTH STRAIGHT & MIRRORED POSITIONS
C          USES NUMR TABLE FOR MAPPING INTO Q ARRAY 
C          USES SIMPLIFIED LOGIC FOR BOUNDARY VALUES, FLOATING PT. ARITH.
C
C PARAMETERS:
C    CIRC1     - FT OF RINGS MULTIPLIED BY WEIGHTS           (SENT)
C    CIRC2     - FT OF RINGS MULTIPLIED BY WEIGHTS           (SENT)
C    LCIRC     - SIZE OF CIRCS ARRAYS                        (SENT)
C    NRING     - NUMBER OF RINGS                             (SENT)
C    MAXRIN    - LONGEST RING                                (SENT)
C    NUMR      - RING LOCATION POINTERS                      (SENT)
C    QMAX      - CC MAX                                      (RETURNED)
C    POS_MAX   - POSITION OF CC MAX                          (RETURNED)
C    QMAXM     - CC MAX MIRRORED                             (RETURNED)
C    POS_MAXM  - POSITION OF CC MAX MIRRORED                 (RETURNED)
C    FFTW3PLAN - REVERSE PLAN FOR RING FFT                   (SENT)
C
C  NOTES: AUG 04 ATTEMPTED SPEEDUP USING 
C       PREMULTIPLY  ARRAYS ie( CIRC12 = CIRC1 * CIRC2) much slower
C       VARIOUS  OTHER ATTEMPTS  FAILED TO YIELD IMPROVEMENT
C       THIS IS A VERY IMPORTANT COMPUTE DEMAND IN ALIGNMENT & REFINE.
C       OPTIONAL LIMIT ON ANGULAR SEARCH SHOULD BE ADDED.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

C       USED IN: oracfmskm.f:              

        SUBROUTINE CROSRNG_MS(CIRC1,CIRC2,LCIRC,NRING,MAXRIN,NUMR,
     &                        QMAX,POS_MAX, QMAXM,POS_MAXM, TT)

C       USES NUMR TABLE FOR MAPPING INTO Q ARRAY 
C       USES SIMPLIFIED LOGIC FOR BOUNDARY VALUES, FLOATING PT. ARITH.
C       USES SPIDER FFT,  HAS TT REMAINS OF SGI FFT 

        INTEGER, INTENT(IN) :: NUMR(3,NRING)
        REAL, INTENT(IN)    :: CIRC1(LCIRC), CIRC2(LCIRC)
        DOUBLE PRECISION    :: QMAX,QMAXM
        REAL, INTENT(OUT)   :: POS_MAX, POS_MAXM
        DOUBLE PRECISION    :: TT(*)

C       AUTOMATIC ARRAYS
        DOUBLE PRECISION    :: T(MAXRIN+2) 
        DOUBLE PRECISION    :: Q(MAXRIN+2) 

C       C - STRAIGHT  = CIRC1 * CONJG(CIRC2), ZERO Q ARRAY
        Q = 0.0D0

C       T - MIRRORED  = CONJG(CIRC1) * CONJG(CIRC2), ZERO T ARRAY
        T = 0.0D0

C       PREMULTIPLY  ARRAYS ie( CIRC12 = CIRC1 * CIRC2) much slower

        DO I=1,NRING
           IGO   = NUMR(2,I)
           NVAL  = NUMR(3,I)
           IGOM1 = IGO - 1
           J1    = 1

           T1   = CIRC1(IGO) * CIRC2(IGO)
           Q(1) = Q(1) + T1
           T(1) = T(1) + T1

           IF (NVAL .NE. MAXRIN) THEN
              T1        = CIRC1(IGO+1) * CIRC2(IGO+1)
              Q(NVAL+1) = Q(NVAL+1) + T1
              T(NVAL+1) = T(NVAL+1) + T1
C             T(NVAL+1) VALUE WAS NOT CHANGED EVEN IN EARLIEST VERSION!!

           ELSE
              T1   = CIRC1(IGO+1) * CIRC2(IGO+1)
              Q(2) = Q(2) + T1
              T(2) = T(2) + T1
           ENDIF

           DO J=J1,NVAL,2
              JC     = J + IGOM1

              C1     = CIRC1(JC)
              C2     = CIRC1(JC+1)
              D1     = CIRC2(JC)
              D2     = CIRC2(JC+1)

              T1     = C1 * D1
              T3     = C1 * D2
              T2     = C2 * D2
              T4     = C2 * D1

              Q(J)   = Q(J)   + T1 + T2
              Q(J+1) = Q(J+1) - T3 + T4
              T(J)   = T(J)   + T1 - T2
              T(J+1) = T(J+1) - T3 - T4
           ENDDO
        ENDDO

C       FOR UN-MIRRORED CASE
        CALL CROSRNG_COM(Q,LCIRC,MAXRIN,QMAX,POS_MAX,TT)

C       FOR MIRRORED CASE
        CALL CROSRNG_COM(T,LCIRC,MAXRIN,QMAXM,POS_MAXM,TT)

        END




