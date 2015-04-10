C++*********************************************************************
C
C FRNGS.F    ADDED FMRS FOR SPEED              FEB 2008 ArDean Leith
C            MODIFED FOR USING FFTW3           MAR 2008 ARDEAN LEITH
C            REMOVED FRNGS_NEW (UNUSED)        MAY 2010 ARDEAN LEITH
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
C  FRNGS(CIRC,LCIRC,NUMR,NRING)
C
C  PURPOSE:  FOURIER TRANSFORM A RADIAL RINGS ARRAY
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE FRNGS(CIRC,LCIRC,NUMR,NRING)

        REAL     :: CIRC(LCIRC)
        INTEGER  :: NUMR(3,NRING)
        INTEGER  :: LCIRC,NRING

        INTEGER  :: LOG2

        DO I=1,NRING
           L = LOG2(NUMR(3,I))            ! LENGTH OF RING
           CALL FFTR_Q(CIRC(NUMR(2,I)),L)
        ENDDO

        END

C       ------------------------- FRNGS_NEWT -------------------------

        SUBROUTINE FRNGS_NEWT(CIRC,LCIRC, NUMR,NRING, 
     &                        SPIDER_SIGN, FFTW_PLANS, OMP_FRNG)

        IMPLICIT NONE

        REAL, INTENT(INOUT)   :: CIRC(LCIRC)
        INTEGER, INTENT(IN)   :: NUMR(3,NRING)
        INTEGER, INTENT(IN)   :: LCIRC,NRING
        LOGICAL, INTENT(IN)   :: SPIDER_SIGN
        INTEGER*8, INTENT(IN) :: FFTW_PLANS(*)!POINTERS TO STRUCTURES
        LOGICAL, INTENT(IN)   :: OMP_FRNG

        INTEGER               :: LOG2
        INTEGER               :: INV,I,LOC,LEN,INDX,IRTFLG
        LOGICAL               :: SPIDER_SCALE = .FALSE.

        INV = +1

C       FOURIER TRANSFORM ON CIRC USING FFTW3 

        IF (OMP_FRNG) THEN
c$omp      parallel do private(i,loc,len,indx)
           DO I=1,NRING
              LOC  = NUMR(2,I)           ! START OF RING
              LEN  = NUMR(3,I) - 2       ! LENGTH OF RING - FOURIER PAD
              INDX = LOG2(LEN) - 1       ! INDEX FOR PLAN
             !write(6,*) ' frngs; i,loc,len,indx: ',i,loc,len,indx

              CALL FMRS(CIRC(LOC), LEN,1,1, FFTW_PLANS(INDX), 
     &                  SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
	   ENDDO
c$omp      end parallel do 

        ELSE
           DO I=1,NRING
              LOC  = NUMR(2,I)           ! START OF RING
              LEN  = NUMR(3,I) - 2       ! LENGTH OF RING - FOURIER PAD
              INDX = LOG2(LEN) - 1       ! INDEX FOR PLAN
              !write(6,*) ' frngs; i,loc,len,indx: ',i,loc,len,indx

              CALL FMRS(CIRC(LOC), LEN,1,1, FFTW_PLANS(INDX), 
     &                  SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
	   ENDDO
        ENDIF


        END



